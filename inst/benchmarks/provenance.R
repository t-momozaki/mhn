# provenance.R --- A one-row record of the software and machine behind a run.
#
# The benchmark and audit scripts of mhn source this file and append its
# columns to their provenance or diagnostics CSV, so that every result travels
# with a record of what produced it:
#
#   mhn_version, mhn_git_commit, mhn_built   the package that was measured
#   r_version, platform, os                  R and the operating system
#   cpu_model, cpu_cores_logical, cpu_cores_physical, cpu_cores_detail,
#   machine_model, memory_gb                 the hardware
#   cxx, cxxflags, cxx17, cxx17flags, cxx_version, makevars_user
#                                            the C++ compiler R uses; mhn
#                                            itself is compiled as C++17
#   blas, lapack                             the linear-algebra libraries
#   start_time, end_time, elapsed_min        the wall-clock span of the run
#   loadavg_start, loadavg_end               the system load around it
#
# An installed package carries no git metadata, so the commit is read from the
# environment variable MHN_GIT_COMMIT (set it to `git rev-parse HEAD` of the
# source tree that was installed), else from a file GIT_COMMIT in the installed
# package, written there by whoever installed it from a git checkout, else from
# the RemoteSha field that remotes and pak write when they install from
# GitHub; otherwise it is recorded as NA.
#
# The compiler entries are those of `R CMD config`, which includes a personal
# Makevars file (~/.R/Makevars or R_MAKEVARS_USER).  They describe the
# toolchain at the time of the run, which is the one that compiled mhn when
# the package was installed on the same machine just before.  makevars_user
# names the personal Makevars file in effect, if any.
#
# Usage:
#   run <- provenance_start()          # at the start of the run
#   ...
#   df  <- with_provenance(df, run)    # appends the columns df lacks
#
# Every probe is wrapped so that a missing tool yields NA, never an error.

provenance_start <- function() {
  list(time = Sys.time(), loadavg = .prov_loadavg())
}

run_provenance <- function(run, end = Sys.time()) {
  start <- if (is.null(run$time)) end else run$time
  mhn <- .prov_mhn()
  commit <- Sys.getenv("MHN_GIT_COMMIT", unset = "")
  if (!nzchar(commit)) commit <- mhn$git_commit
  if (is.na(commit) || !nzchar(commit)) commit <- mhn$remote_sha
  if (is.na(commit) || !nzchar(commit)) commit <- NA_character_
  cxx17 <- .prov_rconfig("CXX17")
  si <- utils::sessionInfo()
  data.frame(
    mhn_version        = mhn$version,
    mhn_git_commit     = commit,
    mhn_built          = mhn$built,
    r_version          = paste(R.version$major, R.version$minor, sep = "."),
    platform           = R.version$platform,
    os                 = if (is.null(si$running)) R.version$os else si$running,
    cpu_model          = .prov_cpu_model(),
    cpu_cores_logical  = .prov_cores(TRUE),
    cpu_cores_physical = .prov_cores(FALSE),
    cpu_cores_detail   = .prov_cores_detail(),
    machine_model      = .prov_machine_model(),
    memory_gb          = .prov_memory_gb(),
    cxx                = .prov_rconfig("CXX"),
    cxxflags           = .prov_rconfig("CXXFLAGS"),
    cxx17              = cxx17,
    cxx17flags         = .prov_rconfig("CXX17FLAGS"),
    cxx_version        = .prov_cxx_version(cxx17),
    makevars_user      = .prov_makevars_user(),
    blas               = if (is.null(si$BLAS)) NA_character_ else si$BLAS,
    lapack             = tryCatch(La_library(), error = function(e) NA_character_),
    start_time         = format(start, "%Y-%m-%dT%H:%M:%S%z"),
    end_time           = format(end, "%Y-%m-%dT%H:%M:%S%z"),
    elapsed_min        = round(as.numeric(difftime(end, start, units = "mins")), 2),
    loadavg_start      = if (is.null(run$loadavg)) NA_character_ else run$loadavg,
    loadavg_end        = .prov_loadavg(),
    stringsAsFactors   = FALSE)
}

# Appends to the one-row data frame `df` every provenance column it does not
# already have, so that the existing columns keep their names and meaning.
with_provenance <- function(df, run, end = Sys.time()) {
  p <- run_provenance(run, end)
  cbind(df, p[setdiff(names(p), names(df))], stringsAsFactors = FALSE)
}

# ---- probes ----------------------------------------------------------------

.prov_lines <- function(cmd, args = character()) {
  out <- tryCatch(
    suppressWarnings(system2(cmd, args, stdout = TRUE, stderr = FALSE)),
    error = function(e) character())
  status <- attr(out, "status")
  if (!is.null(status) && status != 0L) return(character())
  out <- trimws(out)
  out[nzchar(out)]
}

.prov_first <- function(cmd, args = character()) {
  out <- .prov_lines(cmd, args)
  if (length(out)) out[1L] else NA_character_
}

.prov_sysname <- function() Sys.info()[["sysname"]]

.prov_sysctl <- function(key) {
  if (.prov_sysname() != "Darwin") return(NA_character_)
  .prov_first("sysctl", c("-n", key))
}

.prov_read <- function(path) {
  if (!file.exists(path)) return(character())
  tryCatch(readLines(path, warn = FALSE), error = function(e) character())
}

.prov_cpu_model <- function() {
  sys <- .prov_sysname()
  if (sys == "Darwin") return(.prov_sysctl("machdep.cpu.brand_string"))
  if (sys == "Linux") {
    x <- .prov_read("/proc/cpuinfo")
    m <- grep("^(model name|Model|cpu model|Hardware)[[:space:]]*:", x, value = TRUE)
    if (length(m)) return(trimws(sub("^[^:]*:", "", m[1L])))
    m <- grep("^Model name:", .prov_lines("lscpu"), value = TRUE)
    if (length(m)) return(trimws(sub("^[^:]*:", "", m[1L])))
  }
  if (sys == "Windows") {
    v <- Sys.getenv("PROCESSOR_IDENTIFIER")
    if (nzchar(v)) return(v)
  }
  NA_character_
}

.prov_cores <- function(logical) {
  v <- tryCatch(parallel::detectCores(logical = logical),
                error = function(e) NA_integer_)
  if (length(v) != 1L) NA_integer_ else as.integer(v)
}

# Apple silicon splits its cores into performance levels: "16P+4E" is sixteen
# performance and four efficiency cores.  NA elsewhere.
.prov_cores_detail <- function() {
  n <- suppressWarnings(as.integer(.prov_sysctl("hw.nperflevels")))
  if (is.na(n) || n < 1L) return(NA_character_)
  parts <- vapply(seq_len(n), function(i) {
    k <- .prov_sysctl(sprintf("hw.perflevel%d.physicalcpu", i - 1L))
    paste0(k, if (i <= 2L) c("P", "E")[i] else paste0("L", i - 1L))
  }, character(1))
  paste(parts, collapse = "+")
}

.prov_machine_model <- function() {
  sys <- .prov_sysname()
  if (sys == "Darwin") return(.prov_sysctl("hw.model"))
  if (sys == "Linux") {
    x <- .prov_read("/sys/devices/virtual/dmi/id/product_name")
    if (length(x)) return(trimws(x[1L]))
  }
  NA_character_
}

.prov_memory_gb <- function() {
  sys <- .prov_sysname()
  bytes <- NA_real_
  if (sys == "Darwin") {
    bytes <- suppressWarnings(as.numeric(.prov_sysctl("hw.memsize")))
  } else if (sys == "Linux") {
    m <- grep("^MemTotal:", .prov_read("/proc/meminfo"), value = TRUE)
    if (length(m))
      bytes <- 1024 * suppressWarnings(as.numeric(gsub("[^0-9]", "", m[1L])))
  }
  if (is.na(bytes)) NA_real_ else round(bytes / 2^30, 1)
}

.prov_rconfig <- function(var) {
  .prov_first(file.path(R.home("bin"), "R"), c("CMD", "config", var))
}

.prov_cxx_version <- function(cxx) {
  if (is.na(cxx) || !nzchar(cxx)) return(NA_character_)
  .prov_first(strsplit(cxx, "[[:space:]]+")[[1L]][1L], "--version")
}

# The personal Makevars file R reads, abbreviated to ~ so that no user name
# reaches a published CSV.
.prov_makevars_user <- function() {
  env <- Sys.getenv("R_MAKEVARS_USER", unset = "")
  cand <- if (nzchar(env)) env else
    file.path("~", ".R", c(paste0("Makevars-", R.version$platform), "Makevars"))
  cand <- cand[file.exists(path.expand(cand))]
  if (!length(cand)) return(NA_character_)
  p <- path.expand(cand[1L])
  home <- path.expand("~")
  if (startsWith(p, home)) paste0("~", substring(p, nchar(home) + 1L)) else p
}

.prov_loadavg <- function() {
  sys <- .prov_sysname()
  if (sys == "Darwin") {
    v <- .prov_sysctl("vm.loadavg")          # "{ 1.23 2.34 3.45 }"
    if (is.na(v)) return(NA_character_)
    return(trimws(gsub("[{}]", "", v)))
  }
  if (sys == "Linux") {
    x <- .prov_read("/proc/loadavg")
    if (length(x))
      return(paste(strsplit(x[1L], "[[:space:]]+")[[1L]][1:3], collapse = " "))
  }
  NA_character_
}

.prov_mhn <- function() {
  d <- suppressWarnings(tryCatch(utils::packageDescription("mhn"),
                                 error = function(e) NA))
  if (!inherits(d, "packageDescription"))
    return(list(version = NA_character_, built = NA_character_,
                remote_sha = NA_character_, git_commit = NA_character_))
  f <- system.file("GIT_COMMIT", package = "mhn")
  git <- if (nzchar(f)) {
    tryCatch(trimws(readLines(f, n = 1L, warn = FALSE)),
             error = function(e) NA_character_)
  } else {
    NA_character_
  }
  if (!length(git)) git <- NA_character_
  list(version = if (is.null(d$Version)) NA_character_ else d$Version,
       built = if (is.null(d$Built)) NA_character_ else d$Built,
       remote_sha = if (is.null(d$RemoteSha)) NA_character_ else d$RemoteSha,
       git_commit = git)
}
