#include "ah-finder-result.h"

#include "expansion.h"
#include "interpolated-surface-residual.h"

#include <algorithm>
#include <cmath>
#include <exception>
#include <numeric>
#include <utility>
#include <vector>

namespace exahype2::ccz4::ahf {

AHFinderResult BuildAHFinderResult(
  const SpectralBasis& basis,
  const CartesianPointInterpolator& interpolator,
  const ReferencePoint& referencePoint,
  const SpectralNewtonResult& solverResult
) {
  AHFinderResult result{
    .solver = solverResult,
    .referencePoint = referencePoint,
    .lMax = basis.lmax
  };

  if (solverResult.coefficients.size() != basis.nModes
      || basis.rho.size() != basis.nPoints
      || basis.sigma.size() != basis.nPoints) {
    return result;
  }

  try {
    const SynthesizedField field = Synthesize(
      basis, solverResult.coefficients
    );
    if (field.h.size() != basis.nPoints || field.h.empty()) return result;

    AHGuess horizon;
    horizon.points.resize(basis.nPoints);
    for (std::size_t point = 0; point < basis.nPoints; ++point) {
      if (!std::isfinite(field.h[point]) || field.h[point] <= 0.0) {
        return result;
      }
      horizon.points[point] = {
        .r = field.h[point],
        .rho = basis.rho[point],
        .sigma = basis.sigma[point]
      };
    }

    const InterpolatedSurface surface = BuildInterpolatedSurface(
      basis, solverResult.coefficients, interpolator, referencePoint
    );
    if (surface.numberOfPoints() != basis.nPoints) return result;

    std::vector<double> expansion;
    expansion.reserve(basis.nPoints);
    double maxAbsExpansion = 0.0;
    for (const InterpolatedSurfacePoint& point : surface.points) {
      const double theta = computeExpansion(point);
      if (!std::isfinite(theta)) return result;
      expansion.push_back(theta);
      maxAbsExpansion = std::max(maxAbsExpansion, std::abs(theta));
    }

    const std::vector<double> spectralResidual = Project(basis, expansion);
    if (spectralResidual.size() != basis.nModes) return result;
    double spectralResidualInfinityNorm = 0.0;
    for (double component : spectralResidual) {
      if (!std::isfinite(component)) return result;
      spectralResidualInfinityNorm = std::max(
        spectralResidualInfinityNorm, std::abs(component)
      );
    }

    const auto [minimum, maximum] = std::minmax_element(
      field.h.begin(), field.h.end()
    );
    result.horizon = std::move(horizon);
    result.minRadius = *minimum;
    result.meanRadius = std::accumulate(
      field.h.begin(), field.h.end(), 0.0
    ) / static_cast<double>(field.h.size());
    result.maxRadius = *maximum;
    result.maxAbsExpansion = maxAbsExpansion;
    result.spectralResidualInfinityNorm = spectralResidualInfinityNorm;
    result.surfaceAvailable = true;
  } catch (const std::exception&) {
    // A failed or stopped solve can leave coefficients whose trial surface is
    // outside the interpolation domain or otherwise unreconstructable. Solver
    // diagnostics remain available; surface diagnostics retain their NaNs.
  }

  return result;
}

}  // namespace exahype2::ccz4::ahf
