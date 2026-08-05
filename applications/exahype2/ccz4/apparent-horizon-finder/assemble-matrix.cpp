#include "assemble-matrix.h"

#include "ahf-data.h"
#include "expansion.h"

namespace exahype2::ccz4::ahf {

void AssembleGlobalMatrix(const InterpolatedSurface& interpolatedSurface,
                          const AHGuess& ahGuess,
                          std::vector<std::vector<double>>& matrixCoefficients,
                          std::vector<double>& residuals) {
  (void)ahGuess;

  // Compute the Residual = Theta_I: the expansion scalar at each surface point. 
  // The surface point I corresponds to AHGuess::points[I] by construction (InterpolatedSurface::matches).
  // The A/B/C/D calculation for Theta is done by the discretisation-agnostic expansion kernel.
  // s_i and partial_i s_j are already present on each point (filled by the
  // shell-to-surface step).
  const std::size_t numberOfPoints = interpolatedSurface.numberOfPoints();
  residuals.resize(numberOfPoints);
  for (std::size_t i = 0; i < numberOfPoints; ++i) {
    residuals[i] = computeExpansion(interpolatedSurface.points[i]);
  }

  // Global matrix for the Newton update.
  //   FD solver:       (sparse) matrix J[h] with J delta h = -Theta[h] 
  //   Spectral solver: (dense) matrix  J[a] with J delta a = -Theta[a] 
  //
  // TODO(next kernel): assemble the Jacobian. In the finite-difference prototype
  // (ahf/fd_solver.py) this is a coloured sparse finite difference over the angular
  // stencil pattern; it needs the FD angular operators, which are not yet ported.
  // Left empty here so this pass delivers a verified residual in isolation.
  (void)matrixCoefficients;
}
}  // namespace exahype2::ccz4::ahf
