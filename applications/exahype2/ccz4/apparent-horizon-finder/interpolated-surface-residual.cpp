#include "interpolated-surface-residual.h"

#include "assemble-matrix.h"
#include "surface-geometry.h"

#include <cmath>
#include <stdexcept>
#include <string>

namespace exahype2::ccz4::ahf {
namespace {

void ValidateSynthesizedField(
  const SpectralBasis& basis,
  const SynthesizedField& field
) {
  const auto requirePointCount = [&basis](
    const std::vector<double>& values,
    const char* name
  ) {
    if (values.size() != basis.nPoints) {
      throw std::invalid_argument(
        "BuildInterpolatedSurface: synthesized " + std::string(name)
        + " point count does not match basis.nPoints"
      );
    }
  };

  requirePointCount(field.h, "h");
  requirePointCount(field.hRho, "hRho");
  requirePointCount(field.hSigma, "hSigma");
  requirePointCount(field.hRhoRho, "hRhoRho");
  requirePointCount(field.hRhoSigma, "hRhoSigma");
  requirePointCount(field.hSigmaSigma, "hSigmaSigma");

  for (std::size_t point = 0; point < basis.nPoints; ++point) {
    if (!std::isfinite(field.h[point]) || field.h[point] <= 0.0) {
      throw std::invalid_argument(
        "BuildInterpolatedSurface: synthesized h at point "
        + std::to_string(point) + " must be finite and strictly positive"
      );
    }
  }
}

}  // namespace

InterpolatedSurface BuildInterpolatedSurface(
  const SpectralBasis& basis,
  const std::vector<double>& coefficients,
  const CartesianPointInterpolator& interpolator,
  const ReferencePoint& referencePoint
) {
  if (coefficients.size() != basis.nModes) {
    throw std::invalid_argument(
      "BuildInterpolatedSurface: coefficient count does not match basis.nModes"
    );
  }

  const SynthesizedField field = Synthesize(basis, coefficients);
  ValidateSynthesizedField(basis, field);

  AHGuess guess;
  guess.points.resize(basis.nPoints);
  for (std::size_t point = 0; point < basis.nPoints; ++point) {
    guess.points[point] = {
      .r = field.h[point],
      .rho = basis.rho[point],
      .sigma = basis.sigma[point]
    };
  }

  InterpolatedSurface surface;
  InterpolateCartesianToSurface(
    interpolator, guess, referencePoint, surface
  );
  if (surface.numberOfPoints() != basis.nPoints) {
    throw std::runtime_error(
      "BuildInterpolatedSurface: interpolated surface point count does not match basis.nPoints"
    );
  }
  for (std::size_t point = 0; point < surface.numberOfPoints(); ++point) {
    if (!surface.points[point].geometry.spatialMetricDerivativeAvailable) {
      throw std::runtime_error(
        "BuildInterpolatedSurface: physical metric derivatives are unavailable at point "
        + std::to_string(point)
      );
    }
  }

  PopulateSurfaceShapeGeometry(surface, field);
  return surface;
}

SpectralResidual MakeInterpolatedSurfaceResidual(
  const SpectralBasis& basis,
  const CartesianPointInterpolator& interpolator,
  const ReferencePoint& referencePoint
) {
  return [&basis, &interpolator, referencePoint](
    const std::vector<double>& coefficients
  ) -> std::vector<double> {
    const InterpolatedSurface surface = BuildInterpolatedSurface(
      basis, coefficients, interpolator, referencePoint
    );

    std::vector<std::vector<double>> unusedMatrixCoefficients;
    std::vector<double> residuals;
    AssembleGlobalMatrix(
      surface, basis, unusedMatrixCoefficients, residuals
    );
    if (residuals.size() != basis.nModes) {
      throw std::runtime_error(
        "MakeInterpolatedSurfaceResidual: residual count does not match basis.nModes"
      );
    }
    return residuals;
  };
}

}  // namespace exahype2::ccz4::ahf
