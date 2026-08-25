#ifndef APPARENT_HORIZON_FINDER_INTERPOLATED_HORIZON_SOLVER_H_
#define APPARENT_HORIZON_FINDER_INTERPOLATED_HORIZON_SOLVER_H_

#include "ahf-data.h"
#include "harmonics.h"
#include "interpolation.h"
#include "spectral-solver.h"

namespace exahype2::ccz4::ahf {

/**
 * Solve for an apparent horizon from a constant-radius spherical initial guess
 * and geometry supplied by a persistent Cartesian interpolator.
 *
 * The interpolator is referenced throughout the complete Newton solve, so its
 * spatial index is reused by every residual evaluation.
 */
SpectralNewtonResult SolveInterpolatedApparentHorizon(
  const SpectralBasis& basis,
  const CartesianPointInterpolator& interpolator,
  const ReferencePoint& referencePoint,
  double initialRadius,
  const SpectralNewtonOptions& options = {}
);

}  // namespace exahype2::ccz4::ahf

#endif  // APPARENT_HORIZON_FINDER_INTERPOLATED_HORIZON_SOLVER_H_
