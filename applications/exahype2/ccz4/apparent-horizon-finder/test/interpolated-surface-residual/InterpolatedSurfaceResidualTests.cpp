#include "../../interpolated-surface-residual.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace exahype2::ccz4::ahf::tests::interpolated_surface_residual {
namespace {

GridData MakeFlatGrid() {
  GridData grid;
  for (int ix = -2; ix <= 2; ++ix) {
    for (int iy = -2; iy <= 2; ++iy) {
      for (int iz = -2; iz <= 2; ++iz) {
        grid.coordinates.push_back({
          static_cast<double>(ix),
          static_cast<double>(iy),
          static_cast<double>(iz)
        });
        const std::array<double, NumberOfCompactFields> row{
          1.0, 0.0, 0.0, 1.0, 0.0, 1.0,
          0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.0, 1.0
        };
        grid.values.insert(grid.values.end(), row.begin(), row.end());
      }
    }
  }
  assert(grid.isConsistent());
  return grid;
}

CartesianToSurfaceOptions DerivativeOptions() {
  return {
    .method = CartesianToSurfaceMethod::UnweightedPolynomial,
    .numberOfNeighbours = 64,
    .polynomialDegree = 1,
    .allowOrderReduction = false
  };
}

bool Different(
  const std::vector<double>& first,
  const std::vector<double>& second
) {
  for (std::size_t i = 0; i < first.size(); ++i) {
    if (std::abs(first[i] - second[i]) > 1.0e-10) return true;
  }
  return false;
}

void TestStructureOrderingAndShapeGeometry() {
  const GridData grid = MakeFlatGrid();
  const CartesianPointInterpolator interpolator(grid, DerivativeOptions());
  const SpectralBasis basis = BuildSpectralBasis(2);
  std::vector<double> coefficients = SphereCoefficients(basis, 0.8);
  coefficients[1] = 0.03;
  const SynthesizedField field = Synthesize(basis, coefficients);

  const InterpolatedSurface surface = BuildInterpolatedSurface(
    basis, coefficients, interpolator, {0.0, 0.0, 0.0}
  );
  assert(surface.numberOfPoints() == basis.nPoints);
  bool hasNonzeroNormalDerivative = false;
  for (std::size_t point = 0; point < basis.nPoints; ++point) {
    assert(surface.points[point].r == field.h[point]);
    assert(surface.points[point].rho == basis.rho[point]);
    assert(surface.points[point].sigma == basis.sigma[point]);
    assert(surface.points[point].geometry.spatialMetricDerivativeAvailable);
    for (double normal : surface.points[point].surfaceNormal) {
      assert(std::isfinite(normal));
    }
    for (const auto& row : surface.points[point].surfaceNormalDerivative) {
      for (double derivative : row) {
        assert(std::isfinite(derivative));
        hasNonzeroNormalDerivative = hasNonzeroNormalDerivative
          || derivative != 0.0;
      }
    }
  }
  assert(hasNonzeroNormalDerivative);

  const SpectralResidual residual = MakeInterpolatedSurfaceResidual(
    basis, interpolator, {0.0, 0.0, 0.0}
  );
  assert(residual(coefficients).size() == basis.nModes);
}

void TestInputValidation() {
  const GridData grid = MakeFlatGrid();
  const SpectralBasis basis = BuildSpectralBasis(1);
  const CartesianPointInterpolator derivativeInterpolator(
    grid, DerivativeOptions()
  );
  const SpectralResidual derivativeResidual = MakeInterpolatedSurfaceResidual(
    basis, derivativeInterpolator, {0.0, 0.0, 0.0}
  );

  bool wrongCountRejected = false;
  try {
    (void)derivativeResidual(std::vector<double>(basis.nModes - 1, 0.0));
  } catch (const std::invalid_argument& error) {
    wrongCountRejected =
      std::string(error.what()).find("coefficient count") != std::string::npos;
  }
  assert(wrongCountRejected);

  std::vector<double> zeroRadius(basis.nModes, 0.0);
  bool radiusRejected = false;
  try {
    (void)derivativeResidual(zeroRadius);
  } catch (const std::invalid_argument& error) {
    radiusRejected = std::string(error.what()).find("point 0") != std::string::npos
      && std::string(error.what()).find("strictly positive") != std::string::npos;
  }
  assert(radiusRejected);

  std::vector<double> nonFinite = SphereCoefficients(basis, 0.8);
  nonFinite[0] = std::numeric_limits<double>::quiet_NaN();
  bool nonFiniteRejected = false;
  try {
    (void)derivativeResidual(nonFinite);
  } catch (const std::invalid_argument& error) {
    nonFiniteRejected =
      std::string(error.what()).find("must be finite") != std::string::npos;
  }
  assert(nonFiniteRejected);

  const CartesianPointInterpolator nearestInterpolator(
    grid,
    {.method = CartesianToSurfaceMethod::NearestNeighbour,
     .numberOfNeighbours = 1}
  );
  const SpectralResidual nearestResidual = MakeInterpolatedSurfaceResidual(
    basis, nearestInterpolator, {0.0, 0.0, 0.0}
  );
  bool derivativesRejected = false;
  try {
    (void)nearestResidual(SphereCoefficients(basis, 0.8));
  } catch (const std::runtime_error& error) {
    derivativesRejected =
      std::string(error.what()).find("derivatives are unavailable at point 0")
      != std::string::npos;
  }
  assert(derivativesRejected);
}

void TestRepeatabilityAndPerturbation() {
  const GridData grid = MakeFlatGrid();
  const CartesianPointInterpolator interpolator(grid, DerivativeOptions());
  const SpectralBasis basis = BuildSpectralBasis(2);
  const SpectralResidual residual = MakeInterpolatedSurfaceResidual(
    basis, interpolator, {0.0, 0.0, 0.0}
  );
  const std::vector<double> coefficients = SphereCoefficients(basis, 0.8);

  const std::vector<double> first = residual(coefficients);
  const std::vector<double> second = residual(coefficients);
  assert(first == second);

  std::vector<double> perturbed = coefficients;
  perturbed[1] += 0.02;
  const SynthesizedField originalField = Synthesize(basis, coefficients);
  const SynthesizedField perturbedField = Synthesize(basis, perturbed);
  assert(Different(originalField.h, perturbedField.h));
  assert(Different(first, residual(perturbed)));
}

}  // namespace
}  // namespace exahype2::ccz4::ahf::tests::interpolated_surface_residual

int main() {
  using namespace exahype2::ccz4::ahf::tests::interpolated_surface_residual;
  TestStructureOrderingAndShapeGeometry();
  TestInputValidation();
  TestRepeatabilityAndPerturbation();
  std::cout << "All interpolated-surface residual tests completed.\n";
}
