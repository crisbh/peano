#include "../../initial-guess.h"
#include "../../interpolated-horizon-solver.h"
#include "../../interpolated-surface-residual.h"

#include <array>
#include <cassert>
#include <iostream>
#include <type_traits>
#include <utility>
#include <vector>

namespace exahype2::ccz4::ahf::tests::interpolated_horizon_solver {
namespace {

static_assert(std::is_same_v<
  decltype(SolveInterpolatedApparentHorizon(
    std::declval<const SpectralBasis&>(),
    std::declval<const CartesianPointInterpolator&>(),
    std::declval<const ReferencePoint&>(),
    std::declval<double>(),
    std::declval<const SpectralNewtonOptions&>()
  )),
  SpectralNewtonResult
>);

GridData MakeConstantHorizonGrid() {
  GridData grid;
  for (int ix = -2; ix <= 2; ++ix) {
    for (int iy = -2; iy <= 2; ++iy) {
      for (int iz = -2; iz <= 2; ++iz) {
        grid.coordinates.push_back({
          static_cast<double>(ix),
          static_cast<double>(iy),
          static_cast<double>(iz)
        });
        // gamma_ij = delta_ij, K_ij = delta_ij, traceK = 3. Therefore a
        // coordinate sphere has Theta = 2/r - 2 and a root at r = 1.
        const std::array<double, NumberOfCompactFields> row{
          1.0, 0.0, 0.0, 1.0, 0.0, 1.0,
          0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
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

SpectralNewtonResult SolveDirectly(
  const SpectralBasis& basis,
  const CartesianPointInterpolator& interpolator,
  double initialRadius,
  const SpectralNewtonOptions& options
) {
  const std::vector<double> initial =
    MakeSphericalInitialCoefficients(basis, initialRadius);
  assert(initial.size() == basis.nModes);
  const SpectralResidual residual = MakeInterpolatedSurfaceResidual(
    basis, interpolator, {0.0, 0.0, 0.0}
  );
  return SolveSpectralNewton(basis, residual, initial, options);
}

void AssertSameResult(
  const SpectralNewtonResult& actual,
  const SpectralNewtonResult& expected
) {
  assert(actual.coefficients == expected.coefficients);
  assert(actual.converged == expected.converged);
  assert(actual.iterations == expected.iterations);
  assert(actual.residualHistory == expected.residualHistory);
  assert(actual.message == expected.message);
}

void TestSuccessfulOrchestrationAndResultPropagation() {
  const GridData grid = MakeConstantHorizonGrid();
  const CartesianPointInterpolator interpolator(grid, InterpolationOptions());
  const SpectralBasis basis = BuildSpectralBasis(0);
  constexpr double initialRadius = 1.3;

  const SpectralNewtonResult expected = SolveDirectly(
    basis, interpolator, initialRadius, {}
  );
  const SpectralNewtonResult actual = SolveInterpolatedApparentHorizon(
    basis, interpolator, {0.0, 0.0, 0.0}, initialRadius
  );

  AssertSameResult(actual, expected);
  assert(actual.converged);
  assert(actual.coefficients.size() == basis.nModes);
  assert(actual.iterations > 0);
  assert(actual.residualHistory.size()
         == static_cast<std::size_t>(actual.iterations + 1));
}

void TestFailureResultPropagation() {
  const GridData grid = MakeConstantHorizonGrid();
  const CartesianPointInterpolator interpolator(grid, InterpolationOptions());
  const SpectralBasis basis = BuildSpectralBasis(0);
  constexpr double initialRadius = 1.3;
  const SpectralNewtonOptions options{.maxIterations = 0};

  const SpectralNewtonResult expected = SolveDirectly(
    basis, interpolator, initialRadius, options
  );
  const SpectralNewtonResult actual = SolveInterpolatedApparentHorizon(
    basis, interpolator, {0.0, 0.0, 0.0}, initialRadius, options
  );

  AssertSameResult(actual, expected);
  assert(!actual.converged);
  assert(actual.message == "max_iter reached");
}

}  // namespace
}  // namespace exahype2::ccz4::ahf::tests::interpolated_horizon_solver

int main() {
  using namespace exahype2::ccz4::ahf::tests::interpolated_horizon_solver;
  TestSuccessfulOrchestrationAndResultPropagation();
  TestFailureResultPropagation();
  std::cout << "All interpolated-horizon solver tests completed.\n";
}
