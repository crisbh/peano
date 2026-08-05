#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <numbers>
#include <string>

#include "../..//apparent-horizon-finder.h"

namespace exahype2::ccz4::ahf::tests::interpolation {

static_assert(requires(const GridData& gridData,
                       const AHGuess& ahGuess,
                       const ReferencePoint& referencePoint,
                       InterpolatedShell& interpolatedShell,
                       const std::vector<AngularSample>& angularSamples) {
  InterpolateFromCartesianToShell(
    gridData,
    ahGuess,
    referencePoint,
    interpolatedShell
  );
  InterpolateFromCartesianToShell(
    gridData,
    ahGuess,
    referencePoint,
    interpolatedShell,
    angularSamples
  );
});

bool approximatelyEqual(double actual, double expected) {
  constexpr double AbsoluteTolerance = 1.0e-14;
  constexpr double RelativeTolerance = 1.0e-13;
  return std::abs(actual - expected)
      <= AbsoluteTolerance
       + RelativeTolerance * std::max(std::abs(actual), std::abs(expected));
}

void assertCoordinate(const GridCoordinate& actual,
                      double expectedX,
                      double expectedY,
                      double expectedZ) {
  assert(approximatelyEqual(actual.x, expectedX));
  assert(approximatelyEqual(actual.y, expectedY));
  assert(approximatelyEqual(actual.z, expectedZ));
}

void expectSurfaceCoordinateRejected(const AHGuessPoint& surfacePoint,
                                     const ReferencePoint& referencePoint,
                                     const std::string& expectedMessage) {
  bool rejected = false;
  try {
    (void)SurfacePointToCartesian(surfacePoint, referencePoint);
  } catch (const std::invalid_argument& error) {
    rejected = true;
    assert(std::string(error.what()).find(expectedMessage) != std::string::npos);
  }
  assert(rejected);
}

void testSurfacePointToCartesian() {
  constexpr double Pi = std::numbers::pi_v<double>;
  constexpr ReferencePoint Origin{0.0, 0.0, 0.0};

  assertCoordinate(SurfacePointToCartesian({2.0, 0.0, 0.0}, Origin), 2.0, 0.0, 0.0);
  assertCoordinate(SurfacePointToCartesian({2.0, Pi / 2.0, 0.0}, Origin), 0.0, 2.0, 0.0);
  assertCoordinate(SurfacePointToCartesian({2.0, Pi, 0.0}, Origin), -2.0, 0.0, 0.0);

  const double sqrtThree = std::sqrt(3.0);
  assertCoordinate(
    SurfacePointToCartesian({2.0, 0.0, Pi / 6.0}, Origin),
    sqrtThree, 0.0, 1.0
  );
  assertCoordinate(
    SurfacePointToCartesian({2.0, 0.0, -Pi / 6.0}, Origin),
    sqrtThree, 0.0, -1.0
  );

  constexpr ReferencePoint Shift{1.5, -2.0, 0.25};
  assertCoordinate(
    SurfacePointToCartesian({2.0, Pi / 2.0, 0.0}, Shift),
    1.5, 0.0, 0.25
  );

  const GridCoordinate periodicBase =
    SurfacePointToCartesian({1.7, 0.37, -0.21}, Shift);
  const GridCoordinate periodicShift =
    SurfacePointToCartesian({1.7, 0.37 + 2.0 * Pi, -0.21}, Shift);
  assertCoordinate(periodicShift, periodicBase.x, periodicBase.y, periodicBase.z);

  std::cout << "Surface-point spherical coordinates convert to Cartesian targets"
            << std::endl;
}

void testInvalidSurfacePointToCartesianInputs() {
  constexpr ReferencePoint Origin{0.0, 0.0, 0.0};
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double infinity = std::numeric_limits<double>::infinity();

  expectSurfaceCoordinateRejected({0.0, 0.0, 0.0}, Origin, "strictly positive");
  expectSurfaceCoordinateRejected({-1.0, 0.0, 0.0}, Origin, "strictly positive");
  expectSurfaceCoordinateRejected({nan, 0.0, 0.0}, Origin, "must be finite");
  expectSurfaceCoordinateRejected({infinity, 0.0, 0.0}, Origin, "must be finite");
  expectSurfaceCoordinateRejected({1.0, nan, 0.0}, Origin, "must be finite");
  expectSurfaceCoordinateRejected({1.0, 0.0, infinity}, Origin, "must be finite");
  expectSurfaceCoordinateRejected({1.0, 0.0, 0.0}, {nan, 0.0, 0.0}, "must be finite");
  expectSurfaceCoordinateRejected({1.0, 0.0, 0.0}, {0.0, infinity, 0.0}, "must be finite");
}

void expectNeighbourSearchRejected(
  const GridData& gridData,
  const GridCoordinate& target,
  std::size_t neighbourCount,
  const std::string& expectedMessage
) {
  bool rejected = false;
  try {
    (void)FindNearestGridCells(gridData, target, neighbourCount);
  } catch (const std::invalid_argument& error) {
    rejected = true;
    assert(std::string(error.what()).find(expectedMessage) != std::string::npos);
  }
  assert(rejected);
}

void testFindNearestGridCells() {
  GridData gridData;
  gridData.coordinates = {
    {0.0, 0.0, 0.0},  // index 0
    {2.0, 0.0, 0.0},  // index 1
    {1.0, 0.0, 0.0},  // index 2
    {0.0, 2.0, 0.0},  // index 3
    {0.0, 0.0, 3.0}   // index 4
  };

  // Exact coincidence is the unique nearest point.
  assert(FindNearestGridCells(gridData, {1.0, 0.0, 0.0}, 1)
         == std::vector<std::size_t>{2});

  // Distances from this target are 0.04, 0.64, 3.24, 4.04 and 9.04.
  assert(FindNearestGridCells(gridData, {0.2, 0.0, 0.0}, 4)
         == (std::vector<std::size_t>{0, 2, 1, 3}));

  // Index 2 is nearest; indices 0, 1 and 3 tie, so cell index orders the tie.
  assert(FindNearestGridCells(gridData, {1.0, 1.0, 0.0}, 4)
         == (std::vector<std::size_t>{2, 0, 1, 3}));

  std::cout << "Brute-force neighbour search returns sorted deterministic cell indices"
            << std::endl;
}

void testInvalidNearestGridCellInputs() {
  GridData empty;
  GridData twoCells;
  twoCells.coordinates = {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}};
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double infinity = std::numeric_limits<double>::infinity();

  expectNeighbourSearchRejected(empty, {0.0, 0.0, 0.0}, 1, "empty");
  expectNeighbourSearchRejected(twoCells, {0.0, 0.0, 0.0}, 0, "greater than zero");
  expectNeighbourSearchRejected(twoCells, {0.0, 0.0, 0.0}, 3, "only 2 cells");
  expectNeighbourSearchRejected(twoCells, {nan, 0.0, 0.0}, 1, "must be finite");
  expectNeighbourSearchRejected(twoCells, {0.0, infinity, 0.0}, 1, "must be finite");
  expectNeighbourSearchRejected(twoCells, {0.0, 0.0, -infinity}, 1, "must be finite");
}

std::vector<double> evaluateLinearField(
  const std::vector<GridCoordinate>& coordinates,
  double constant,
  const Vector3& derivative
) {
  std::vector<double> values;
  values.reserve(coordinates.size());
  for (const GridCoordinate& coordinate : coordinates) {
    values.push_back(
      constant + derivative[0] * coordinate.x
      + derivative[1] * coordinate.y + derivative[2] * coordinate.z
    );
  }
  return values;
}

void assertLinearFit(
  const std::vector<GridCoordinate>& coordinates,
  double constant,
  const Vector3& derivative,
  const GridCoordinate& target
) {
  const std::vector<double> values =
    evaluateLinearField(coordinates, constant, derivative);
  const LinearPolynomialFit fit = FitLinearPolynomial(coordinates, values, target);
  const double expectedValue =
    constant + derivative[0] * target.x
    + derivative[1] * target.y + derivative[2] * target.z;
  assert(approximatelyEqual(fit.value, expectedValue));
  for (std::size_t direction = 0; direction < 3; ++direction) {
    assert(approximatelyEqual(fit.derivative[direction], derivative[direction]));
  }
  assert(std::isfinite(fit.lengthScale));
  assert(fit.lengthScale > 0.0);
}

void testFitLinearPolynomialExactReproduction() {
  const std::vector<GridCoordinate> regular{
    {-1.0, -1.0, -1.0}, {1.0, -1.0, -1.0},
    {-1.0, 1.0, -1.0}, {-1.0, -1.0, 1.0},
    {1.0, 1.0, 1.0}, {2.0, -0.5, 0.75}
  };
  const GridCoordinate origin{0.0, 0.0, 0.0};

  assertLinearFit(regular, 7.0, {0.0, 0.0, 0.0}, origin);
  assertLinearFit(regular, 2.0, {3.0, -4.0, 5.0}, origin);
  assertLinearFit(regular, -1.0, {6.0, 0.0, 0.0}, origin);
  assertLinearFit(regular, -1.0, {0.0, 6.0, 0.0}, origin);
  assertLinearFit(regular, -1.0, {0.0, 0.0, 6.0}, origin);
  assertLinearFit(regular, 2.0, {3.0, -4.0, 5.0}, {4.0, -3.0, 2.0});

  const std::vector<GridCoordinate> irregular{
    {-2.1, 0.3, 1.7}, {0.2, -1.4, 0.1}, {3.0, 2.2, -0.8},
    {-0.4, 1.1, -2.5}, {1.7, -0.6, 2.9}, {-3.2, -2.0, 0.4},
    {0.9, 3.1, 1.2}, {2.5, 0.4, -3.3}
  };
  assertLinearFit(irregular, 2.0, {3.0, -4.0, 5.0}, {0.7, -0.2, 1.3});

  std::cout << "Degree-1 QR reconstruction exactly reproduces analytic linear fields"
            << std::endl;
}

void expectLinearFitRejected(
  const std::vector<GridCoordinate>& coordinates,
  const std::vector<double>& values,
  const std::string& expectedMessage
) {
  bool rejected = false;
  try {
    (void)FitLinearPolynomial(coordinates, values, {0.0, 0.0, 0.0});
  } catch (const std::invalid_argument& error) {
    rejected = true;
    assert(std::string(error.what()).find(expectedMessage) != std::string::npos);
  }
  assert(rejected);
}

void testFitLinearPolynomialInvalidStencils() {
  const std::vector<GridCoordinate> tooFew{
    {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}
  };
  expectLinearFitRejected(tooFew, {1.0, 2.0, 3.0}, "at least four");

  const std::vector<GridCoordinate> coplanar{
    {-1.0, -1.0, 0.0}, {1.0, -1.0, 0.0},
    {-1.0, 1.0, 0.0}, {1.0, 1.0, 0.0}, {0.25, -0.5, 0.0}
  };
  expectLinearFitRejected(
    coplanar,
    evaluateLinearField(coplanar, 2.0, {3.0, -4.0, 5.0}),
    "rank deficient"
  );
}

void setCompactCell(
  GridData& gridData,
  std::size_t cellIndex,
  const std::array<double, NumberOfCompactCCZ4Fields>& values
) {
  std::copy(
    values.begin(),
    values.end(),
    gridData.values.begin() + cellIndex * GridData::numberOfVariables
  );
}

void assertTensor(const Tensor3& actual, const Tensor3& expected) {
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      assert(approximatelyEqual(actual[i][j], expected[i][j]));
    }
  }
}

void testInterpolateToSurfaceNearestNeighbour() {
  constexpr double Pi = std::numbers::pi_v<double>;
  constexpr ReferencePoint Reference{10.0, -4.0, 2.0};
  GridData gridData;
  gridData.resize(3);
  gridData.coordinates = {
    {11.0, -4.0, 2.0},  // target for first AH point
    {10.0, -3.0, 2.0},  // target for second AH point
    {30.0, 30.0, 30.0}  // deliberately irrelevant cell
  };
  setCompactCell(gridData, 0, {
    4.0, 2.0, 6.0, 8.0, 10.0, 12.0,
    1.0, 2.0, 3.0, 4.0, 5.0, 6.0,
    3.0, 2.0
  });
  setCompactCell(gridData, 1, {
    9.0, 3.0, 6.0, 12.0, 15.0, 18.0,
    2.0, 4.0, 6.0, 8.0, 10.0, 12.0,
    6.0, 3.0
  });
  setCompactCell(gridData, 2, {
    101.0, 102.0, 103.0, 104.0, 105.0, 106.0,
    107.0, 108.0, 109.0, 110.0, 111.0, 112.0,
    113.0, 4.0
  });

  const AHGuess ahGuess{.points = {
    {.r = 1.0, .rho = 0.0, .sigma = 0.0},
    {.r = 1.0, .rho = Pi / 2.0, .sigma = 0.0}
  }};
  InterpolatedSurface surface;
  InterpolateToSurfaceNearestNeighbour(gridData, ahGuess, Reference, surface);

  assert(surface.matches(ahGuess));
  assert(surface.points.size() == 2);

  const Tensor3 expectedMetric0{{
    {{1.0, 0.5, 1.5}}, {{0.5, 2.0, 2.5}}, {{1.5, 2.5, 3.0}}
  }};
  const Tensor3 expectedCurvature0{{
    {{1.25, 1.0, 2.25}}, {{1.0, 3.0, 3.75}}, {{2.25, 3.75, 4.5}}
  }};
  const Tensor3 expectedMetric1{{
    {{1.0, 1.0 / 3.0, 2.0 / 3.0}},
    {{1.0 / 3.0, 4.0 / 3.0, 5.0 / 3.0}},
    {{2.0 / 3.0, 5.0 / 3.0, 2.0}}
  }};
  const Tensor3 expectedCurvature1{{
    {{20.0 / 9.0, 10.0 / 9.0, 2.0}},
    {{10.0 / 9.0, 32.0 / 9.0, 40.0 / 9.0}},
    {{2.0, 40.0 / 9.0, 16.0 / 3.0}}
  }};

  assertTensor(surface.points[0].geometry.spatialMetric, expectedMetric0);
  assertTensor(surface.points[0].geometry.extrinsicCurvature, expectedCurvature0);
  assert(surface.points[0].geometry.traceK == 3.0);
  assertTensor(surface.points[1].geometry.spatialMetric, expectedMetric1);
  assertTensor(surface.points[1].geometry.extrinsicCurvature, expectedCurvature1);
  assert(surface.points[1].geometry.traceK == 6.0);

  const Tensor3Derivative zeroDerivative{};
  for (const InterpolatedSurfacePoint& point : surface.points) {
    assert(point.geometry.spatialMetricDerivative == zeroDerivative);
    assert(!point.geometry.spatialMetricDerivativeAvailable);
    assert(point.surfaceNormal == Vector3{});
    assert(point.surfaceNormalDerivative == Tensor3{});
  }

  std::cout << "Nearest-neighbour sampling populates surface physical geometry in AH order"
            << std::endl;
}

void testInterpolateFromCartesianToShell() {
  std::cout << "Placeholder test for InterpolateFromCartesianToShell" << std::endl;
}

void testInterpolateFromShellToSurface() {
  std::cout << "Placeholder test for InterpolateFromShellToSurface" << std::endl;
}

}  // namespace exahype2::ccz4::ahf::tests::interpolation

int main() {
  using namespace exahype2::ccz4::ahf::tests::interpolation;

  std::cout << "Apparent Horizon Finder interpolation tests begin..." << std::endl;

  testSurfacePointToCartesian();
  testInvalidSurfacePointToCartesianInputs();
  testFindNearestGridCells();
  testInvalidNearestGridCellInputs();
  testFitLinearPolynomialExactReproduction();
  testFitLinearPolynomialInvalidStencils();
  testInterpolateToSurfaceNearestNeighbour();
  testInterpolateFromCartesianToShell();
  testInterpolateFromShellToSurface();

  std::cout << "All interpolation tests completed." << std::endl;
  return 0;
}
