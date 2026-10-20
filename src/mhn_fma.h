// mhn_fma.h
//
// x * y + z rounded once.
//
// The package forms the exact residual of a product, and a few differences
// of nearly equal terms, with a fused multiply-add, and relies on its single
// rounding.  std::fma gives it where the C library computes it in hardware or
// correctly in software.  On x86_64 Windows the fma that Rtools45 links
// statically, mingw-w64 v11's, does not: it adds the four partial products of
// 26-bit halves of x and y to z, rounding at each step.  On Windows the
// package therefore uses fma_soft, on every architecture, the others not
// having been checked; it forms x * y + z exactly in integers and rounds it
// once, to nearest with ties to even.

#ifndef MHN_FMA_H
#define MHN_FMA_H

#include <cmath>

namespace mhn {

double fma_soft(double x, double y, double z);

// MHN_FMA_SOFT selects fma_soft on any platform, to test the Windows path
// where the C library's fma is correct.
inline double fma(double x, double y, double z) {
#if defined(_WIN32) || defined(MHN_FMA_SOFT)
  return fma_soft(x, y, z);
#else
  return std::fma(x, y, z);
#endif
}

}  // namespace mhn

#endif  // MHN_FMA_H
