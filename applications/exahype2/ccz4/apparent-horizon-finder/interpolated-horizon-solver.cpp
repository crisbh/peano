#include "interpolated-horizon-solver.h"

#include "initial-guess.h"
#include "interpolated-surface-residual.h"

#include <cmath>
#include <stdexcept>

namespace exahype2::ccz4::ahf {

SpectralNewtonResult SolveInterpolatedApparentHorizon(
  const SpectralBasis& basis,
  const CartesianPointInterpolator& interpolator,
  const ReferencePoint& referencePoint,
  double initialRadius,
  const SpectralNewtonOptions& options
) {
  if (!std::isfinite(initialRadius) || initialRadius <= 0.0) {
    throw std::invalid_argument(
      "SolveInterpolatedApparentHorizon: initialRadius must be finite and strictly positive"
    );
  }

  const std::vector<double> initialCoefficients =
    MakeSphericalInitialCoefficients(basis, initialRadius);
  const SpectralResidual residual = MakeInterpolatedSurfaceResidual(
    basis, interpolator, referencePoint
  );
  return SolveSpectralNewton(
    basis, residual, initialCoefficients, options
  );
}

}  // namespace exahype2::ccz4::ahf
