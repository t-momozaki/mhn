// Rejection samplers from Sun et al. (2023): Algorithm 1 (gamma > 0,
// alpha > 1; Normal or sqrt-Gamma proposal) and Algorithm 3
// (gamma <= 0; AM-GM proposal).  Each setup is built once per
// parameter triple via build_sun_algoX(), then reused for repeated
// sample_sun_algoX() draws.

#ifndef MHN_SUN_H
#define MHN_SUN_H

#include <Rcpp.h>

namespace mhn {

// Setup for Sun et al. (2023) Algorithm 1 (gamma > 0, alpha > 1).
// Theorems 1 and 2 of that paper establish the acceptance bounds.
struct SunAlgo1Setup {
  enum Proposal { NORMAL, SQRT_GAMMA } chosen;
  double alpha, beta, gamma;
  double mu_opt;          // optimal Normal-proposal mean
  double delta_opt;       // optimal sqrt-Gamma-proposal rate
  double sigma;           // = 1 / sqrt(2*beta), Normal-proposal SD
  double log_K1_minus_K2; // diagnostic: <0 -> NORMAL, >=0 -> SQRT_GAMMA
};

// Setup for Sun et al. (2023) Algorithm 3 (gamma <= 0, alpha > 0).
// Theorems 3 and 4 of that paper establish the acceptance bounds.
struct SunAlgo3Setup {
  double alpha, beta, gamma_abs;
  double m;                  // matching point (after Newton refinement)
  double r;                  // = (beta*m + |gamma|) / (2*beta*m + |gamma|)
  double shape;              // = alpha * r
  double rate;               // = m * (beta*m + |gamma|)
  double m_betam_gam;        // = m * (beta*m + |gamma|), pre-cached
  double m_init;             // diagnostic: m before Newton refinement
  bool used_inflex_heuristic;// true when alpha > 1.1 elaborate path is taken
};

// Build the Algorithm 1 setup.  Throws Rcpp::stop on alpha <= 1 or gamma <= 0.
SunAlgo1Setup build_sun_algo1(double alpha, double beta, double gamma);

// Draw one sample using the Algorithm 1 setup.  retries_out, if non-null,
// is incremented by the number of accept/reject retries.
double sample_sun_algo1(const SunAlgo1Setup& s, int* retries_out = nullptr);

// Build the Algorithm 3 setup.  Throws Rcpp::stop on gamma > 0.
SunAlgo3Setup build_sun_algo3(double alpha, double beta, double gamma);

// Draw one sample using the Algorithm 3 setup.
double sample_sun_algo3(const SunAlgo3Setup& s, int* retries_out = nullptr);

}  // namespace mhn

#endif  // MHN_SUN_H
