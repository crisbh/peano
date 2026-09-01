#include "../../ah-finder-result.h"
#include "../../initial-guess.h"
#include "../../iteration.h"

#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>

namespace exahype2::ccz4::ahf::tests::result {
namespace {

GridData MakeSphericalHorizonGrid() {
  GridData grid;
  for (int x = -2; x <= 2; ++x) {
    for (int y = -2; y <= 2; ++y) {
      for (int z = -2; z <= 2; ++z) {
        grid.coordinates.push_back({
          static_cast<double>(x),
          static_cast<double>(y),
          static_cast<double>(z)
        });
        const std::array<double, NumberOfCompactFields> row{
          1.0, 0.0, 0.0, 1.0, 0.0, 1.0,
          1.0, 0.0, 0.0, 1.0, 0.0, 1.0,
          3.0, 1.0
        };
        grid.values.insert(grid.values.end(), row.begin(), row.end());
      }
    }
  }
  return grid;
}

CartesianToSurfaceOptions InterpolationOptions() {
  return {
    .method = CartesianToSurfaceMethod::UnweightedPolynomial,
    .numberOfNeighbours = 64,
    .polynomialDegree = 1,
    .allowOrderReduction = false
  };
}

void TestConvergedResultAndCanonicalSurface() {
  const GridData grid = MakeSphericalHorizonGrid();
  const CartesianPointInterpolator interpolator(grid, InterpolationOptions());
  const SpectralBasis basis = BuildSpectralBasis(0);
  const ReferencePoint referencePoint{0.0, 0.0, 0.0};
  const SpectralNewtonResult solver = SolvePDEforAH(
    basis, interpolator, referencePoint, 1.3
  );
  assert(solver.converged);

  const AHFinderResult result = BuildAHFinderResult(
    basis, interpolator, referencePoint, solver
  );
  assert(result.surfaceAvailable);
  assert(result.solver.coefficients == solver.coefficients);
  assert(result.solver.residualHistory == solver.residualHistory);
  assert(result.solver.converged == solver.converged);
  assert(result.solver.iterations == solver.iterations);
  assert(result.solver.message == solver.message);
  assert(result.referencePoint.x == referencePoint.x);
  assert(result.referencePoint.y == referencePoint.y);
  assert(result.referencePoint.z == referencePoint.z);
  assert(result.lMax == basis.lmax);

  const SynthesizedField field = Synthesize(basis, solver.coefficients);
  assert(result.horizon.points.size() == basis.nPoints);
  for (std::size_t point = 0; point < basis.nPoints; ++point) {
    assert(result.horizon.points[point].r == field.h[point]);
    assert(result.horizon.points[point].rho == basis.rho[point]);
    assert(result.horizon.points[point].sigma == basis.sigma[point]);
  }

  assert(result.minRadius == field.h.front());
  assert(result.meanRadius == field.h.front());
  assert(result.maxRadius == field.h.front());
  assert(std::isfinite(result.maxAbsExpansion));
  assert(result.maxAbsExpansion < 1.0e-8);
  assert(std::isfinite(result.spectralResidualInfinityNorm));
  assert(result.spectralResidualInfinityNorm < 1.0e-8);
}

void TestNonConvergedButReconstructableResult() {
  const GridData grid = MakeSphericalHorizonGrid();
  const CartesianPointInterpolator interpolator(grid, InterpolationOptions());
  const SpectralBasis basis = BuildSpectralBasis(2);
  const ReferencePoint referencePoint{0.25, -0.5, 0.75};
  SpectralNewtonResult solver{
    .coefficients = MakeSphericalInitialCoefficients(basis, 1.2),
    .converged = false,
    .iterations = 3,
    .residualHistory = {2.0, 1.0, 0.5},
    .message = "controlled test failure"
  };

  const AHFinderResult result = BuildAHFinderResult(
    basis, interpolator, referencePoint, solver
  );
  assert(!result.solver.converged);
  assert(result.solver.residualHistory == solver.residualHistory);
  assert(result.solver.message == solver.message);
  assert(result.surfaceAvailable);
  assert(result.referencePoint.x == referencePoint.x);
  assert(result.referencePoint.y == referencePoint.y);
  assert(result.referencePoint.z == referencePoint.z);
  assert(std::abs(result.meanRadius - 1.2) < 1.0e-12);
  assert(std::isfinite(result.maxAbsExpansion));
  assert(std::isfinite(result.spectralResidualInfinityNorm));
  const SynthesizedField field = Synthesize(basis, solver.coefficients);
  assert(result.horizon.points.size() == basis.nPoints);
  for (std::size_t point = 0; point < basis.nPoints; ++point) {
    assert(result.horizon.points[point].r == field.h[point]);
    assert(result.horizon.points[point].rho == basis.rho[point]);
    assert(result.horizon.points[point].sigma == basis.sigma[point]);
  }
}

void TestUnavailableSurfacePreservesSolverDiagnostics() {
  const GridData grid = MakeSphericalHorizonGrid();
  const CartesianPointInterpolator interpolator(grid, InterpolationOptions());
  const SpectralBasis basis = BuildSpectralBasis(2);
  const SpectralNewtonResult solver{
    .coefficients = {},
    .converged = false,
    .iterations = 4,
    .residualHistory = {3.0, 2.0, 1.0},
    .message = "no usable coefficients"
  };

  const AHFinderResult result = BuildAHFinderResult(
    basis, interpolator, {1.0, 2.0, 3.0}, solver
  );
  assert(!result.surfaceAvailable);
  assert(result.horizon.points.empty());
  assert(result.solver.coefficients.empty());
  assert(result.solver.residualHistory == solver.residualHistory);
  assert(result.solver.message == solver.message);
  assert(std::isnan(result.minRadius));
  assert(std::isnan(result.meanRadius));
  assert(std::isnan(result.maxRadius));
  assert(std::isnan(result.maxAbsExpansion));
  assert(std::isnan(result.spectralResidualInfinityNorm));
}

}  // namespace
}  // namespace exahype2::ccz4::ahf::tests::result

int main() {
  using namespace exahype2::ccz4::ahf::tests::result;
  TestConvergedResultAndCanonicalSurface();
  TestNonConvergedButReconstructableResult();
  TestUnavailableSurfacePreservesSolverDiagnostics();
  std::cout << "All AHFinderResult tests completed.\n";
}
