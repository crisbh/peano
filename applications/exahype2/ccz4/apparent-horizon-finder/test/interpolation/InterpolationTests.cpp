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
  testInterpolateFromCartesianToShell();
  testInterpolateFromShellToSurface();

  std::cout << "All interpolation tests completed." << std::endl;
  return 0;
}
