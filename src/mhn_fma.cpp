// x * y + z rounded once, in integer arithmetic (see mhn_fma.h).
//
// Each finite nonzero operand is m 2^e with an integer m < 2^53.  The product
// mx my has at most 106 bits and is held exactly in 128; it and mz are each
// shifted so that their leading bit is bit 127, the smaller in magnitude is
// shifted right to the exponent of the larger, and the bits shifted out are
// kept as a sticky bit in its lowest place.  The sum or difference is then
// exact down to that place, more than 70 bits below where it is rounded, so
// the sticky bit decides the rounding as the bits it stands for would: the
// operand shifted by more than one place cannot cancel more than one leading
// bit of the other, and the one shifted by at most one place loses nothing.

#include "mhn_fma.h"

#include <Rcpp.h>

#include <cmath>
#include <cstdint>
#include <cstring>

namespace {

typedef std::uint64_t u64;

struct U128 {
  u64 hi;
  u64 lo;
};

u64 bits_of(double x) {
  u64 b;
  std::memcpy(&b, &x, sizeof b);
  return b;
}

double double_of(u64 b) {
  double x;
  std::memcpy(&x, &b, sizeof x);
  return x;
}

int clz64(u64 x) {
  if (x == 0) return 64;
  int n = 0;
  if (!(x >> 32)) { n += 32; x <<= 32; }
  if (!(x >> 48)) { n += 16; x <<= 16; }
  if (!(x >> 56)) { n += 8;  x <<= 8; }
  if (!(x >> 60)) { n += 4;  x <<= 4; }
  if (!(x >> 62)) { n += 2;  x <<= 2; }
  if (!(x >> 63)) { n += 1; }
  return n;
}

int clz128(const U128& v) {
  return v.hi ? clz64(v.hi) : 64 + clz64(v.lo);
}

U128 shift_left(const U128& v, int s) {
  if (s == 0) return v;
  if (s >= 64) return U128{v.lo << (s - 64), 0};
  return U128{(v.hi << s) | (v.lo >> (64 - s)), v.lo << s};
}

// v shifted right by s, with the bits shifted out ORed into the lowest bit.
U128 shift_right_sticky(const U128& v, int s) {
  if (s == 0) return v;
  if (s >= 128) return U128{0, (v.hi | v.lo) ? u64(1) : u64(0)};
  U128 r;
  u64 lost;
  if (s >= 64) {
    const int t = s - 64;
    r.hi = 0;
    r.lo = t ? (v.hi >> t) : v.hi;
    lost = v.lo | (t ? (v.hi << (64 - t)) : 0);
  } else {
    r.hi = v.hi >> s;
    r.lo = (v.lo >> s) | (v.hi << (64 - s));
    lost = v.lo << (64 - s);
  }
  if (lost) r.lo |= 1;
  return r;
}

// The 128-bit product of two 64-bit integers, from 32-bit halves.
U128 multiply(u64 a, u64 b) {
  const u64 a0 = a & 0xffffffffu, a1 = a >> 32;
  const u64 b0 = b & 0xffffffffu, b1 = b >> 32;
  const u64 p00 = a0 * b0, p01 = a0 * b1, p10 = a1 * b0, p11 = a1 * b1;
  const u64 mid = (p00 >> 32) + (p01 & 0xffffffffu) + (p10 & 0xffffffffu);
  return U128{p11 + (p01 >> 32) + (p10 >> 32) + (mid >> 32),
              (mid << 32) | (p00 & 0xffffffffu)};
}

bool less(const U128& a, const U128& b) {
  return a.hi != b.hi ? a.hi < b.hi : a.lo < b.lo;
}

// A finite nonzero double as m 2^e, m < 2^53.
void decompose(double x, u64* m, int* e) {
  const u64 b = bits_of(x);
  const int field = static_cast<int>((b >> 52) & 0x7ff);
  const u64 frac = b & ((u64(1) << 52) - 1);
  if (field == 0) {
    *m = frac;
    *e = -1074;
  } else {
    *m = frac | (u64(1) << 52);
    *e = field - 1075;
  }
}

}  // namespace

namespace mhn {

double fma_soft(double x, double y, double z) {
  // Infinities and NaN: a finite product added to an infinite z is z; any
  // other case is what the IEEE operations give, x * y being exact or
  // infinite there.  A zero factor makes the product an exact signed zero.
  if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)) {
    if (std::isfinite(x) && std::isfinite(y) && std::isinf(z)) return z;
    const double p = x * y;
    return p + z;
  }
  if (x == 0.0 || y == 0.0) {
    const double p = x * y;
    return p + z;
  }
  // A zero addend leaves the product, rounded once by the multiplication;
  // its sign is that of the product even where it rounds to zero.
  if (z == 0.0) return x * y;

  const u64 sign_p = (bits_of(x) ^ bits_of(y)) >> 63;
  const u64 sign_z = bits_of(z) >> 63;
  u64 mx, my, mz;
  int ex, ey, ez;
  decompose(x, &mx, &ex);
  decompose(y, &my, &ey);
  decompose(z, &mz, &ez);

  // Each term as M 2^E with the leading bit of M at bit 127.
  U128 mp = multiply(mx, my);
  const int lp = clz128(mp);
  mp = shift_left(mp, lp);
  int ep = ex + ey - lp;
  U128 mzz = U128{0, mz};
  const int lz = clz128(mzz);
  mzz = shift_left(mzz, lz);
  int ezz = ez - lz;

  // a is the term of the larger magnitude.
  bool p_larger = (ep > ezz) || (ep == ezz && !less(mp, mzz));
  U128 a = p_larger ? mp : mzz;
  U128 b = p_larger ? mzz : mp;
  int e = p_larger ? ep : ezz;
  const int eb = p_larger ? ezz : ep;
  const u64 sign = p_larger ? sign_p : sign_z;
  b = shift_right_sticky(b, e - eb);

  U128 s;
  if (sign_p == sign_z) {
    s.lo = a.lo + b.lo;
    const u64 carry_lo = s.lo < a.lo;
    const u64 hi = a.hi + b.hi;
    const u64 carry_hi = hi < a.hi;
    s.hi = hi + carry_lo;
    if (carry_hi || s.hi < hi) {
      // The sum reached 2^128: halve it, keeping the bit shifted out.
      const u64 out = s.lo & 1;
      s.lo = (s.lo >> 1) | (s.hi << 63) | out;
      s.hi = (s.hi >> 1) | (u64(1) << 63);
      e += 1;
    }
  } else {
    s.lo = a.lo - b.lo;
    s.hi = a.hi - b.hi - (a.lo < b.lo ? 1 : 0);
    if (s.hi == 0 && s.lo == 0) return 0.0;  // an exact zero is +0
    const int l = clz128(s);
    s = shift_left(s, l);
    e -= l;
  }

  // The leading bit of s stands for 2^k.  A normal result keeps 53 bits, so
  // 75 bits are dropped; below 2^-1022 the result is subnormal and keeps
  // fewer.
  const int k = e + 127;
  int drop = 75;
  if (k < -1022) drop += -1022 - k;

  u64 m = 0;
  bool round_bit = false, sticky = false;
  if (drop <= 128) {
    // t = s >> (drop - 1), with the bits below it gathered into sticky.
    // drop is at least 75, so t lies in the high word, shifted by 10 to 63.
    const int r = drop - 1;
    const u64 t = s.hi >> (r - 64);
    sticky = s.lo != 0 || (s.hi << (128 - r)) != 0;
    round_bit = (t & 1) != 0;
    m = t >> 1;
  }
  // drop > 128: the value is below half the smallest subnormal and rounds
  // to zero.
  if (round_bit && (sticky || (m & 1))) m += 1;

  if (k >= -1022) {
    int biased = k + 1023;
    if (m == (u64(1) << 53)) {
      m >>= 1;
      biased += 1;
    }
    if (biased >= 2047) return sign ? -INFINITY : INFINITY;
    return double_of((sign << 63) | (u64(biased) << 52) |
                     (m & ((u64(1) << 52) - 1)));
  }
  // Subnormal: m < 2^52, or 2^52 where it has rounded up to the smallest
  // normal double, whose bit pattern it then is.
  return double_of((sign << 63) | m);
}

}  // namespace mhn

// Test hook: fma_soft element by element, so that the tests can check it on
// every platform against x * y + z rounded once in Rmpfr.
// [[Rcpp::export(.fma_soft, rng = false)]]
Rcpp::NumericVector fma_soft_R(Rcpp::NumericVector x, Rcpp::NumericVector y,
                               Rcpp::NumericVector z) {
  const R_xlen_t n = x.size();
  if (y.size() != n || z.size() != n) {
    Rcpp::stop("x, y and z must have the same length");
  }
  Rcpp::NumericVector out(n);
  for (R_xlen_t i = 0; i < n; ++i) out[i] = mhn::fma_soft(x[i], y[i], z[i]);
  return out;
}
