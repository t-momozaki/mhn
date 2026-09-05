#ifndef MHN_PARAM_SLOTS_H
#define MHN_PARAM_SLOTS_H

// One cache slot per distinct parameter triple in a recycled call.
//
// The vectorised evaluators recycle alpha, beta and gamma independently, so the
// triple at element i depends only on i modulo L = lcm(na, nb, ng): the same L
// triples repeat for the whole call.  A single most-recently-used cache hits
// only when consecutive elements agree, which they do not under recycling --
// dmhn(x, alpha = c(a1, a2)) alternates, so every element missed and every
// element recomputed the Fox-Wright normalising constant, the most expensive
// quantity in the package.  Keeping one slot per residue turns that into L
// computations for the whole call however long it is.
//
// The stored triple is compared before a slot is reused.  Under recycling that
// comparison always succeeds, but it costs three doubles and makes the reuse
// correct by inspection rather than by an argument about periods.
//
// L is capped: it is a least common multiple and can run to millions for
// coprime lengths, while the whole point is to hold a small table.  Past the cap
// the class keeps a single slot and behaves exactly as the previous code did.

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>

namespace mhn {

inline R_xlen_t recycling_period(R_xlen_t na, R_xlen_t nb, R_xlen_t ng,
                                 R_xlen_t cap) {
  auto lcm_capped = [cap](R_xlen_t a, R_xlen_t b) -> R_xlen_t {
    if (a <= 0 || b <= 0) return cap + 1;
    R_xlen_t x = a, y = b;
    while (y != 0) { const R_xlen_t t = x % y; x = y; y = t; }  // gcd
    const R_xlen_t g = x;
    if (a / g > cap / std::max<R_xlen_t>(b, 1)) return cap + 1;  // would overflow the cap
    return (a / g) * b;
  };
  const R_xlen_t l1 = lcm_capped(na, nb);
  if (l1 > cap) return cap + 1;
  return lcm_capped(l1, ng);
}

template <class Cache>
class ParamSlots {
 public:
  // n is the length of the call; there is no point holding more slots than
  // elements.  max_slots bounds the table for coprime parameter lengths.
  ParamSlots(R_xlen_t na, R_xlen_t nb, R_xlen_t ng, R_xlen_t n,
             R_xlen_t max_slots = 256) {
    const R_xlen_t period = recycling_period(na, nb, ng, max_slots);
    period_ = (period <= max_slots && period <= n) ? period : 1;
    slots_.resize(static_cast<std::size_t>(period_));
    filled_.assign(static_cast<std::size_t>(period_), false);
    a_.assign(static_cast<std::size_t>(period_), 0.0);
    b_.assign(static_cast<std::size_t>(period_), 0.0);
    g_.assign(static_cast<std::size_t>(period_), 0.0);
  }

  // The cache for element i, recomputed only if this slot does not already
  // hold this triple.
  Cache& at(R_xlen_t i, double a, double b, double g) {
    const std::size_t k = static_cast<std::size_t>(i % period_);
    if (!filled_[k] || a_[k] != a || b_[k] != b || g_[k] != g) {
      slots_[k].recompute(a, b, g);
      a_[k] = a; b_[k] = b; g_[k] = g;
      filled_[k] = true;
    }
    return slots_[k];
  }

 private:
  R_xlen_t period_ = 1;
  std::vector<Cache> slots_;
  std::vector<bool> filled_;
  std::vector<double> a_, b_, g_;
};

}  // namespace mhn

#endif  // MHN_PARAM_SLOTS_H
