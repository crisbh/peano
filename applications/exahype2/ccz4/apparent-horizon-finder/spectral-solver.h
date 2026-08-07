#ifndef APPARENT_HORIZON_FINDER_SPECTRAL_SOLVER_H_
#define APPARENT_HORIZON_FINDER_SPECTRAL_SOLVER_H_

#include <functional>
#include <string>
#include <vector>

#include "harmonics.h"

namespace exahype2::ccz4::ahf {

/**
 * Residual R(a): the spherical-harmonic coefficients of the expansion Theta for a
 * given coefficient vector a. Injected into the Newton solve so the same solver
 * serves grid-backed residuals (synthesize -> interpolate -> AssembleGlobalMatrix)
 * and analytic residuals used in tests.
 */
using SpectralResidual = std::function<std::vector<double>(const std::vector<double>&)>;

/**
 * Parameters for the damped Newton solve. Defaults mirror the prototype
 * ahf/solver.py.
 */
struct SpectralNewtonOptions {
  double tolerance = 1.0e-10;            ///< convergence threshold on ||R||_inf
  double massScale = 1.0;               ///< effective tolerance = tolerance / massScale
  int maxIterations = 40;
  double finiteDifferenceStep = 1.0e-6;  ///< relative step for the Jacobian
  int maxLineSearchSteps = 20;
  double maxRadiusFactor = 6.0;         ///< runaway cap on max(h) as a multiple of the initial mean radius
};

/**
 * Outcome of the Newton solver.
 */
struct SpectralNewtonResult {
  std::vector<double> coefficients;
  bool converged = false;
  int iterations = 0;
  std::vector<double> residualHistory;  ///< ||R||_inf per iteration
  std::string message;
};

/**
 * Damped Newton-Raphson solve of the spectral apparent-horizon equations R(a) = 0.
 *
 * Each step builds a dense Jacobian J_pq = dR_p/da_q by central finite differences
 * over the coefficients, solves J delta a = -R by LU, and applies a backtracking
 * line search on ||R||_2. A physical guard rejects any step whose synthesised h
 * leaves the valid star-shaped range (h > 0 everywhere, max h below a cap set from
 * the initial guess), which catches the spurious large-radius root where Theta -> 0
 * as the surface grows without bound.
 *
 * @param basis    the spectral basis (used by the guard to synthesise h).
 * @param residual R(a); must return a vector of length initialCoefficients.size().
 * @param initialCoefficients the starting coefficient vector (e.g. a coordinate sphere).
 */
SpectralNewtonResult SolveSpectralNewton(const SpectralBasis& basis,
                                         const SpectralResidual& residual,
                                         const std::vector<double>& initialCoefficients,
                                         const SpectralNewtonOptions& options = {});

}  // namespace exahype2::ccz4::ahf

#endif  // APPARENT_HORIZON_FINDER_SPECTRAL_SOLVER_H_
