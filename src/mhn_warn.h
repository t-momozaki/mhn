// mhn_warn.h
//
// Raising an R warning from C++ so that the C++ frames are unwound when the
// warning does not return.

#ifndef MHN_WARN_H
#define MHN_WARN_H

#include <Rcpp.h>
#include <string>
#include <utility>

namespace mhn {

// Raises an R warning whose message is formatted as by Rcpp::warning().
//
// Rcpp::warning() calls Rf_warning().  When the warning is caught by an
// exiting handler, as in tryCatch(..., warning = ), or turned into an error by
// options(warn = 2), R leaves Rf_warning() by a long jump, and the destructors
// of the C++ objects alive in the calling frames do not run: an Rcpp vector
// among them stays protected for the rest of the session, and memory owned by
// a C++ container is not released.  Here the warning is raised by evaluating
// R's warning() instead, which Rcpp does under unwind protection (its default
// from Rcpp 1.0.10): the jump becomes a C++ exception, the frames unwind, and
// the jump is resumed where the call returns to R.  With an earlier Rcpp a
// warning caught by an exiting handler still leaves by a long jump.
//
// The warning carries no call, the entry point being an internal function.
template <typename... Args>
inline void warn(const char* fmt, Args&&... args) {
  const std::string msg = tfm::format(fmt, std::forward<Args>(args)...);
  const Rcpp::Function r_warning =
      Rcpp::Environment::base_namespace()["warning"];
  r_warning(msg, Rcpp::Named("call.") = false);
}

}  // namespace mhn

#endif  // MHN_WARN_H
