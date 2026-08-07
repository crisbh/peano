#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <numbers>
#include <random>
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

void testAnalyticCoordinateGeometry() {
  const Tensor3 expectedX{{
    {{0.0, 0.0, 0.0}}, {{0.0, 4.0, 0.0}}, {{0.0, 0.0, 4.0}}
  }};
  const Tensor3 expectedY{{
    {{4.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}, {{0.0, 0.0, 4.0}}
  }};
  const Tensor3 expectedZ{{
    {{4.0, 0.0, 0.0}}, {{0.0, 4.0, 0.0}}, {{0.0, 0.0, 0.0}}
  }};
  assert(ComputeRadialProjectionTensor({2.0, 0.0, 0.0}) == expectedX);
  assert(ComputeRadialProjectionTensor({0.0, 2.0, 0.0}) == expectedY);
  assert(ComputeRadialProjectionTensor({0.0, 0.0, 2.0}) == expectedZ);

  const Tensor3 expectedGeneral{{
    {{13.0, -2.0, -3.0}}, {{-2.0, 10.0, -6.0}}, {{-3.0, -6.0, 5.0}}
  }};
  const Tensor3 general = ComputeRadialProjectionTensor({1.0, 2.0, 3.0});
  assert(general == expectedGeneral);
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      assert(general[i][j] == general[j][i]);
    }
  }

  const AnalyticCoordinateGeometry originGeometry =
    ComputeAnalyticCoordinateGeometry({1.0, 2.0, 3.0}, {0.0, 0.0, 0.0});
  const AnalyticCoordinateGeometry shiftedGeometry =
    ComputeAnalyticCoordinateGeometry({5.0, 1.0, 5.0}, {4.0, -1.0, 2.0});
  assert(originGeometry.T == shiftedGeometry.T);
  assert(originGeometry.gradient == shiftedGeometry.gradient);
  assert(originGeometry.hessian == shiftedGeometry.hessian);
}

void testAnalyticCoordinateDerivativesAgainstFiniteDifferences() {
  constexpr ReferencePoint Reference{1.5, -0.75, 0.4};
  const GridCoordinate global =
    SurfacePointToCartesian({2.0, 0.6, 0.35}, Reference);
  const AnalyticCoordinateGeometry analytic =
    ComputeAnalyticCoordinateGeometry(global, Reference);
  constexpr double Step = 1.0e-5;

  const auto angles = [&Reference](const GridCoordinate& point) {
    const double x = point.x - Reference.x;
    const double y = point.y - Reference.y;
    const double z = point.z - Reference.z;
    return std::array<double, 2>{
      std::atan2(y, x), std::atan2(z, std::hypot(x, y))
    };
  };
  for (std::size_t direction = 0; direction < 3; ++direction) {
    GridCoordinate plus = global;
    GridCoordinate minus = global;
    if (direction == 0) {
      plus.x += Step;
      minus.x -= Step;
    } else if (direction == 1) {
      plus.y += Step;
      minus.y -= Step;
    } else {
      plus.z += Step;
      minus.z -= Step;
    }
    const auto plusAngles = angles(plus);
    const auto minusAngles = angles(minus);
    for (std::size_t u = 0; u < 2; ++u) {
      const double finiteDifference =
        (plusAngles[u] - minusAngles[u]) / (2.0 * Step);
      assert(std::abs(analytic.gradient[u][direction] - finiteDifference) < 1.0e-9);
    }

    const AnalyticCoordinateGeometry plusGeometry =
      ComputeAnalyticCoordinateGeometry(plus, Reference);
    const AnalyticCoordinateGeometry minusGeometry =
      ComputeAnalyticCoordinateGeometry(minus, Reference);
    for (std::size_t u = 0; u < 2; ++u) {
      for (std::size_t component = 0; component < 3; ++component) {
        const double finiteDifference =
          (plusGeometry.gradient[u][component]
           - minusGeometry.gradient[u][component]) / (2.0 * Step);
        assert(std::abs(
          analytic.hessian[u][component][direction] - finiteDifference) < 1.0e-8);
      }
    }
  }

  bool poleRejected = false;
  try {
    (void)ComputeAnalyticCoordinateGeometry(
      {Reference.x + 1.0e-14, Reference.y, Reference.z + 1.0}, Reference);
  } catch (const std::domain_error& error) {
    poleRejected = std::string(error.what()).find("cos(sigma)=0") != std::string::npos;
  }
  assert(poleRejected);

  std::cout << "Analytic angular coordinate geometry matches Cartesian finite differences"
            << std::endl;
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

void testAcceleratedNeighbourSearchMatchesBruteForce() {
  std::mt19937_64 generator(4815162342ULL);
  std::uniform_real_distribution<double> coordinate(-5.0, 5.0);
  GridData randomGrid;
  randomGrid.coordinates.resize(2000);
  for (GridCoordinate& point : randomGrid.coordinates) {
    point = {coordinate(generator), coordinate(generator), coordinate(generator)};
  }
  const CartesianSpatialIndex randomIndex(randomGrid);
  for (std::size_t query = 0; query < 30; ++query) {
    const GridCoordinate target{
      coordinate(generator), coordinate(generator), coordinate(generator)};
    for (std::size_t k : {
           std::size_t{1}, std::size_t{2}, std::size_t{5},
           std::size_t{16}, std::size_t{64}}) {
      assert(FindNearestGridCells(randomIndex, target, k)
             == FindNearestGridCellsBruteForce(randomGrid, target, k));
    }
  }

  GridData structured;
  for (int z = -2; z <= 2; ++z) {
    for (int y = -2; y <= 2; ++y) {
      for (int x = -2; x <= 2; ++x) {
        structured.coordinates.push_back(
          {static_cast<double>(x), static_cast<double>(y), static_cast<double>(z)});
      }
    }
  }
  const CartesianSpatialIndex structuredIndex(structured);
  for (const GridCoordinate target : {
         GridCoordinate{0.0, 0.0, 0.0},
         GridCoordinate{0.5, 0.5, 0.5},
         GridCoordinate{1.0, -1.0, 2.0}}) {
    for (std::size_t k : {
           std::size_t{1}, std::size_t{4}, std::size_t{12},
           structured.coordinates.size()}) {
      assert(FindNearestGridCells(structuredIndex, target, k)
             == FindNearestGridCellsBruteForce(structured, target, k));
    }
  }

  // The origin is an exact hit; the six axis neighbours then tie by cell index.
  const auto exact = FindNearestGridCells(structuredIndex, {0.0, 0.0, 0.0}, 7);
  const auto exactReference =
    FindNearestGridCellsBruteForce(structured, {0.0, 0.0, 0.0}, 7);
  assert(exact == exactReference);
  assert(structured.coordinates[exact.front()].x == 0.0);
  assert(structured.coordinates[exact.front()].y == 0.0);
  assert(structured.coordinates[exact.front()].z == 0.0);

  std::cout << "K-d tree matches brute force for random, structured, exact, and tied queries"
            << std::endl;
}

void reportNeighbourSearchPerformance() {
  std::mt19937_64 generator(8675309ULL);
  std::uniform_real_distribution<double> coordinate(-10.0, 10.0);
  GridData gridData;
  gridData.coordinates.resize(20000);
  for (GridCoordinate& point : gridData.coordinates) {
    point = {coordinate(generator), coordinate(generator), coordinate(generator)};
  }
  std::vector<GridCoordinate> targets(50);
  for (GridCoordinate& target : targets) {
    target = {coordinate(generator), coordinate(generator), coordinate(generator)};
  }

  const auto buildStart = std::chrono::steady_clock::now();
  const CartesianSpatialIndex spatialIndex(gridData);
  const auto buildEnd = std::chrono::steady_clock::now();
  const auto acceleratedStart = std::chrono::steady_clock::now();
  for (const GridCoordinate& target : targets) {
    (void)FindNearestGridCells(spatialIndex, target, 8);
  }
  const auto acceleratedEnd = std::chrono::steady_clock::now();
  const auto bruteStart = std::chrono::steady_clock::now();
  for (const GridCoordinate& target : targets) {
    (void)FindNearestGridCellsBruteForce(gridData, target, 8);
  }
  const auto bruteEnd = std::chrono::steady_clock::now();

  const auto microseconds = [](auto start, auto end) {
    return std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
  };
  std::cout << "Neighbour-search timing (informational only): build="
            << microseconds(buildStart, buildEnd) << " us, k-d queries="
            << microseconds(acceleratedStart, acceleratedEnd) << " us, brute queries="
            << microseconds(bruteStart, bruteEnd) << " us" << std::endl;
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

void testWeightedLinearPolynomialReconstruction() {
  const std::vector<GridCoordinate> irregular{
    {-2.1, 0.3, 1.7}, {0.2, -1.4, 0.1}, {3.0, 2.2, -0.8},
    {-0.4, 1.1, -2.5}, {1.7, -0.6, 2.9}, {-3.2, -2.0, 0.4},
    {0.9, 3.1, 1.2}, {2.5, 0.4, -3.3}
  };
  const GridCoordinate target{0.7, -0.2, 1.3};
  const std::vector<double> values =
    evaluateLinearField(irregular, 2.0, {3.0, -4.0, 5.0});
  const LinearReconstructionOptions weighted{
    .method = LinearReconstructionMethod::WeightedMovingLeastSquares,
    .weightSupport = 0.75
  };
  const LinearPolynomialFit weightedFit =
    FitLinearPolynomial(irregular, values, target, weighted);
  const LinearPolynomialFit unweightedFit =
    FitLinearPolynomial(irregular, values, target);

  assert(approximatelyEqual(weightedFit.value, 11.4));
  assert(approximatelyEqual(weightedFit.derivative[0], 3.0));
  assert(approximatelyEqual(weightedFit.derivative[1], -4.0));
  assert(approximatelyEqual(weightedFit.derivative[2], 5.0));
  assert(approximatelyEqual(weightedFit.value, unweightedFit.value));
  for (std::size_t direction = 0; direction < 3; ++direction) {
    assert(approximatelyEqual(
      weightedFit.derivative[direction], unweightedFit.derivative[direction]));
  }

  // A coincident sample has distance zero and therefore exactly unit weight.
  const std::vector<GridCoordinate> coincident{
    target, {1.2, -0.2, 1.3}, {0.7, 0.3, 1.3}, {0.7, -0.2, 1.8},
    {0.1, -0.8, 0.7}
  };
  assert(LinearReconstructionWeight(0.0, 1.0, weighted) == 1.0);
  assertLinearFit(coincident, 2.0, {3.0, -4.0, 5.0}, target);
  const auto coincidentValues =
    evaluateLinearField(coincident, 2.0, {3.0, -4.0, 5.0});
  const LinearPolynomialFit coincidentFit =
    FitLinearPolynomial(coincident, coincidentValues, target, weighted);
  assert(approximatelyEqual(coincidentFit.value, 11.4));
  assert(approximatelyEqual(coincidentFit.derivative[0], 3.0));
  assert(approximatelyEqual(coincidentFit.derivative[1], -4.0));
  assert(approximatelyEqual(coincidentFit.derivative[2], 5.0));
}

void testWeightedFitFavoursNearbySamples() {
  const std::vector<GridCoordinate> coordinates{
    {0.0, 0.0, 0.0}, {0.1, 0.0, 0.0},
    {0.0, 0.1, 0.0}, {0.0, 0.0, 0.1},
    {-1.0, 0.0, 0.0}, {0.0, -1.0, 0.0},
    {0.0, 0.0, -1.0}, {1.0, 1.0, 1.0}
  };
  const std::vector<double> noisyValues{0.0, 0.0, 0.0, 0.0, 10.0, 10.0, 10.0, 10.0};
  const LinearPolynomialFit unweighted =
    FitLinearPolynomial(coordinates, noisyValues, {0.0, 0.0, 0.0});
  const LinearPolynomialFit weighted = FitLinearPolynomial(
    coordinates,
    noisyValues,
    {0.0, 0.0, 0.0},
    {.method = LinearReconstructionMethod::WeightedMovingLeastSquares,
     .weightSupport = 0.5}
  );
  assert(std::abs(weighted.value) < std::abs(unweighted.value));
}

void testWeightedLinearPolynomialFailures() {
  const std::vector<GridCoordinate> fullRank{
    {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0},
    {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}
  };
  const std::vector<double> values{1.0, 2.0, 3.0, 4.0};
  for (double support : {0.0, -1.0, std::numeric_limits<double>::quiet_NaN()}) {
    bool rejected = false;
    try {
      (void)FitLinearPolynomial(
        fullRank,
        values,
        {0.0, 0.0, 0.0},
        {.method = LinearReconstructionMethod::WeightedMovingLeastSquares,
         .weightSupport = support}
      );
    } catch (const std::invalid_argument&) {
      rejected = true;
    }
    assert(rejected);
  }

  const std::vector<GridCoordinate> coplanar{
    {-1.0, -1.0, 0.0}, {1.0, -1.0, 0.0},
    {-1.0, 1.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 0.0, 0.0}
  };
  bool rankRejected = false;
  try {
    (void)FitLinearPolynomial(
      coplanar,
      evaluateLinearField(coplanar, 2.0, {3.0, -4.0, 5.0}),
      {0.0, 0.0, 0.0},
      {.method = LinearReconstructionMethod::WeightedMovingLeastSquares}
    );
  } catch (const std::invalid_argument& error) {
    rankRejected = std::string(error.what()).find("rank deficient") != std::string::npos;
  }
  assert(rankRejected);

  std::cout << "Weighted MLS handles locality, coincidence, and invalid stencils"
            << std::endl;
}

void testFitCompactLinearPolynomial() {
  const GridCoordinate target{0.7, -0.2, 1.3};
  GridData gridData;
  gridData.coordinates = {
    {-2.1, 0.3, 1.7}, {0.2, -1.4, 0.1}, {3.0, 2.2, -0.8},
    {-0.4, 1.1, -2.5}, {1.7, -0.6, 2.9}, {-3.2, -2.0, 0.4},
    {0.9, 3.1, 1.2}, {2.5, 0.4, -3.3}
  };
  gridData.resize(gridData.coordinates.size());

  std::array<double, NumberOfCompactFields> constants{};
  std::array<Vector3, NumberOfCompactFields> derivatives{};
  for (std::size_t field = 0; field < NumberOfCompactFields; ++field) {
    constants[field] = 100.0 + 7.0 * static_cast<double>(field);
    derivatives[field] = {
      1.0 + static_cast<double>(field),
      -2.0 - 2.0 * static_cast<double>(field),
      3.0 + 0.5 * static_cast<double>(field)
    };
  }
  for (std::size_t cell = 0; cell < gridData.coordinates.size(); ++cell) {
    const GridCoordinate& coordinate = gridData.coordinates[cell];
    for (std::size_t field = 0; field < NumberOfCompactFields; ++field) {
      gridData.values[cell * GridData::numberOfVariables + field] =
        constants[field] + derivatives[field][0] * coordinate.x
        + derivatives[field][1] * coordinate.y
        + derivatives[field][2] * coordinate.z;
    }
  }

  const std::vector<std::size_t> neighbours{7, 0, 5, 2, 4, 1, 6, 3};
  const CompactLinearFit fit = FitCompactLinearPolynomial(
    gridData,
    neighbours,
    target,
    {.method = LinearReconstructionMethod::WeightedMovingLeastSquares,
     .weightSupport = 0.75}
  );
  const CompactLinearFit unweightedFit =
    FitCompactLinearPolynomial(gridData, neighbours, target);

  assert(fit.numberOfNeighbours == neighbours.size());
  assert(fit.polynomialDegree == 1);
  assert(fit.rank == 4);
  assert(fit.lengthScale > 0.0);
  for (std::size_t field = 0; field < NumberOfCompactFields; ++field) {
    const double expectedValue =
      constants[field] + derivatives[field][0] * target.x
      + derivatives[field][1] * target.y
      + derivatives[field][2] * target.z;
    assert(approximatelyEqual(fit.values[field], expectedValue));
    assert(approximatelyEqual(fit.values[field], unweightedFit.values[field]));
    for (std::size_t direction = 0; direction < 3; ++direction) {
      assert(approximatelyEqual(
        fit.gradients[field][direction], derivatives[field][direction]));
      assert(approximatelyEqual(
        fit.gradients[field][direction], unweightedFit.gradients[field][direction]));
    }
  }

  // Named compact indices must still identify K, direct phi, and tensor order.
  assert(approximatelyEqual(
    fit.values[static_cast<std::size_t>(CompactField::TraceK)],
    constants[12] + derivatives[12][0] * target.x
      + derivatives[12][1] * target.y + derivatives[12][2] * target.z));
  assert(approximatelyEqual(
    fit.values[static_cast<std::size_t>(CompactField::Phi)],
    constants[13] + derivatives[13][0] * target.x
      + derivatives[13][1] * target.y + derivatives[13][2] * target.z));
  for (std::size_t component = 0;
       component < NumberOfSymmetricTensorComponents;
       ++component) {
    assert(ConformalMetricCompactFields[component]
           == static_cast<CompactField>(component));
    assert(ConformalExtrinsicCurvatureCompactFields[component]
           == static_cast<CompactField>(component + 6));
  }

  std::cout << "One QR factorisation reconstructs all fourteen compact fields"
            << std::endl;
}

void setCompactCell(
  GridData& gridData,
  std::size_t cellIndex,
  const std::array<double, NumberOfCompactFields>& values
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

void testInterpolateCartesianToSurfaceDegree1() {
  constexpr double Pi = std::numbers::pi_v<double>;
  constexpr ReferencePoint Reference{1.0, -2.0, 0.5};
  GridData gridData;
  gridData.coordinates = {
    {-1.0, -4.0, -1.0}, {3.0, -4.0, -1.0},
    {-1.0, 0.0, -1.0}, {3.0, 0.0, -1.0},
    {-1.0, -4.0, 2.0}, {3.0, -4.0, 2.0},
    {-1.0, 0.0, 2.0}, {3.0, 0.0, 2.0}
  };
  gridData.resize(gridData.coordinates.size());

  std::array<double, NumberOfCompactFields> constants{};
  std::array<Vector3, NumberOfCompactFields> gradients{};
  constants = {
    4.0, 0.2, 0.3, 5.0, 0.4, 6.0,
    0.1, 0.2, 0.3, 0.4, 0.5, 0.6,
    3.0, 2.0
  };
  for (std::size_t field = 0; field < NumberOfCompactFields; ++field) {
    const double scale = static_cast<double>(field + 1);
    gradients[field] = {0.01 * scale, -0.005 * scale, 0.002 * scale};
  }
  gradients[12] = {0.2, 0.1, -0.1};
  gradients[13] = {0.05, -0.04, 0.03};

  for (std::size_t cell = 0; cell < gridData.coordinates.size(); ++cell) {
    const GridCoordinate& coordinate = gridData.coordinates[cell];
    for (std::size_t field = 0; field < NumberOfCompactFields; ++field) {
      gridData.values[cell * GridData::numberOfVariables + field] =
        constants[field] + gradients[field][0] * coordinate.x
        + gradients[field][1] * coordinate.y
        + gradients[field][2] * coordinate.z;
    }
  }

  const AHGuess ahGuess{.points = {
    {.r = 0.5, .rho = 0.0, .sigma = 0.0},
    {.r = 0.5, .rho = Pi / 2.0, .sigma = 0.0},
    {.r = 0.5, .rho = 0.0, .sigma = Pi / 6.0}
  }};
  InterpolatedSurface mlsSurface;
  InterpolateCartesianToSurface(
    gridData,
    ahGuess,
    Reference,
    mlsSurface,
    {.method = CartesianToSurfaceMethod::WeightedDegree1,
     .numberOfNeighbours = 8,
     .weightSupport = 0.8}
  );
  assert(mlsSurface.matches(ahGuess));

  constexpr std::array<std::array<std::size_t, 3>, 3> tensorComponent{{
    {{0, 1, 2}}, {{1, 3, 4}}, {{2, 4, 5}}
  }};
  for (std::size_t pointIndex = 0; pointIndex < ahGuess.points.size(); ++pointIndex) {
    const GridCoordinate target =
      SurfacePointToCartesian(ahGuess.points[pointIndex], Reference);
    std::array<double, NumberOfCompactFields> targetValues{};
    for (std::size_t field = 0; field < NumberOfCompactFields; ++field) {
      targetValues[field] = constants[field] + gradients[field][0] * target.x
        + gradients[field][1] * target.y + gradients[field][2] * target.z;
    }
    const double phi = targetValues[13];
    const double inversePhiSquared = 1.0 / (phi * phi);
    const double inversePhiCubed = inversePhiSquared / phi;
    const auto& geometry = mlsSurface.points[pointIndex].geometry;
    assert(geometry.spatialMetricDerivativeAvailable);
    assert(approximatelyEqual(geometry.traceK, targetValues[12]));
    for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 3; ++j) {
        const std::size_t component = tensorComponent[i][j];
        const double gammaTilde = targetValues[component];
        const double ATilde = targetValues[component + 6];
        assert(approximatelyEqual(
          geometry.spatialMetric[i][j], gammaTilde * inversePhiSquared));
        assert(approximatelyEqual(
          geometry.extrinsicCurvature[i][j],
          (ATilde + targetValues[12] * gammaTilde / 3.0) * inversePhiSquared));
        for (std::size_t k = 0; k < 3; ++k) {
          const double expectedDerivative =
            gradients[component][k] * inversePhiSquared
            - 2.0 * gammaTilde * gradients[13][k] * inversePhiCubed;
          assert(approximatelyEqual(
            geometry.spatialMetricDerivative[k][i][j], expectedDerivative));
        }
      }
    }
  }

  InterpolatedSurface nearestSurface;
  InterpolateCartesianToSurface(
    gridData,
    ahGuess,
    Reference,
    nearestSurface,
    {.method = CartesianToSurfaceMethod::NearestNeighbour, .numberOfNeighbours = 1}
  );
  assert(!nearestSurface.points[0].geometry.spatialMetricDerivativeAvailable);
  assert(!approximatelyEqual(
    nearestSurface.points[0].geometry.spatialMetric[0][0],
    mlsSurface.points[0].geometry.spatialMetric[0][0]));

  std::cout << "Degree-1 MLS reconstructs physical surface geometry and Cartesian derivatives"
            << std::endl;
}

void testInterpolateCartesianToSurfaceContextualFailures() {
  GridData gridData;
  gridData.resize(4);
  gridData.coordinates = {
    {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0},
    {0.0, 1.0, 0.0}, {1.0, 1.0, 0.0}
  };
  for (std::size_t cell = 0; cell < gridData.coordinates.size(); ++cell) {
    for (std::size_t field = 0; field < NumberOfCompactFields; ++field) {
      gridData.values[cell * GridData::numberOfVariables + field] = 1.0;
    }
    gridData.values[cell * GridData::numberOfVariables + 13] = 2.0;
  }

  bool rankFailureHasContext = false;
  try {
    InterpolatedSurface surface;
    InterpolateCartesianToSurface(
      gridData,
      {.points = {{.r = 1.0, .rho = 0.0, .sigma = 0.0}}},
      {0.0, 0.0, 0.0},
      surface,
      {.method = CartesianToSurfaceMethod::UnweightedDegree1,
       .numberOfNeighbours = 4}
    );
  } catch (const std::runtime_error& error) {
    const std::string message = error.what();
    rankFailureHasContext = message.find("surface point 0") != std::string::npos
                         && message.find("target (1") != std::string::npos;
  }
  assert(rankFailureHasContext);

  bool invalidTargetHasContext = false;
  try {
    InterpolatedSurface surface;
    InterpolateCartesianToSurface(
      gridData,
      {.points = {
        {.r = 1.0, .rho = 0.0, .sigma = 0.0},
        {.r = std::numeric_limits<double>::quiet_NaN(), .rho = 0.0, .sigma = 0.0}
      }},
      {0.0, 0.0, 0.0},
      surface,
      {.method = CartesianToSurfaceMethod::NearestNeighbour,
       .numberOfNeighbours = 1}
    );
  } catch (const std::runtime_error& error) {
    const std::string message = error.what();
    invalidTargetHasContext = message.find("surface point 1") != std::string::npos
                           && message.find("target (") != std::string::npos;
  }
  assert(invalidTargetHasContext);
}

GridData makeValidationGrid() {
  GridData gridData;
  gridData.resize(4);
  gridData.coordinates = {
    {0.0, 0.0, 0.0}, {2.0, 0.0, 0.0},
    {1.0, 1.0, 0.0}, {1.0, 0.0, 1.0}
  };
  for (std::size_t cell = 0; cell < gridData.coordinates.size(); ++cell) {
    const std::array<double, NumberOfCompactFields> values{
      1.0, 0.0, 0.0, 1.0, 0.0, 1.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 2.0
    };
    setCompactCell(gridData, cell, values);
  }
  return gridData;
}

std::string interpolationFailure(
  const GridData& gridData,
  const CartesianToSurfaceOptions& options
) {
  try {
    InterpolatedSurface surface;
    InterpolateCartesianToSurface(
      gridData,
      {.points = {{.r = 1.0, .rho = 0.0, .sigma = 0.0}}},
      {0.0, 0.0, 0.0},
      surface,
      options
    );
  } catch (const std::exception& error) {
    return error.what();
  }
  return {};
}

void testInterpolationValidationAndDiagnostics() {
  const CartesianToSurfaceOptions degreeOne{
    .method = CartesianToSurfaceMethod::UnweightedDegree1,
    .numberOfNeighbours = 4
  };

  GridData malformed = makeValidationGrid();
  malformed.values.pop_back();
  assert(interpolationFailure(malformed, degreeOne).find("times fourteen")
         != std::string::npos);

  GridData nonFinite = makeValidationGrid();
  nonFinite.values[2 * GridData::numberOfVariables + 4] =
    std::numeric_limits<double>::quiet_NaN();
  const std::string nonFiniteMessage = interpolationFailure(nonFinite, degreeOne);
  assert(nonFiniteMessage.find("surface point 0") != std::string::npos);
  assert(nonFiniteMessage.find("non-finite compact value") != std::string::npos);

  GridData invalidPhi = makeValidationGrid();
  for (std::size_t cell = 0; cell < invalidPhi.coordinates.size(); ++cell) {
    invalidPhi.values[cell * GridData::numberOfVariables + 13] = -1.0;
  }
  const std::string phiMessage = interpolationFailure(
    invalidPhi,
    {.method = CartesianToSurfaceMethod::NearestNeighbour,
     .numberOfNeighbours = 1}
  );
  assert(phiMessage.find("surface point 0") != std::string::npos);
  assert(phiMessage.find("strictly positive") != std::string::npos);

  GridData insufficient = makeValidationGrid();
  insufficient.coordinates.resize(3);
  insufficient.values.resize(3 * GridData::numberOfVariables);
  const CartesianToSurfaceOptions fallbackDisabled{
    .method = CartesianToSurfaceMethod::WeightedDegree1,
    .numberOfNeighbours = 3,
    .allowOrderReduction = false
  };
  const std::string insufficientMessage =
    interpolationFailure(insufficient, fallbackDisabled);
  assert(insufficientMessage.find("insufficient neighbours") != std::string::npos);
  assert(insufficientMessage.find("surface point 0") != std::string::npos);

  InterpolatedSurface reducedSurface;
  InterpolateCartesianToSurface(
    insufficient,
    {.points = {{.r = 1.0, .rho = 0.0, .sigma = 0.0}}},
    {0.0, 0.0, 0.0},
    reducedSurface,
    {.method = CartesianToSurfaceMethod::WeightedDegree1,
     .numberOfNeighbours = 3,
     .allowOrderReduction = true}
  );
  const InterpolationDiagnostics& reduced =
    reducedSurface.points[0].interpolationDiagnostics;
  assert(reduced.requestedMethod == CartesianToSurfaceMethod::WeightedDegree1);
  assert(reduced.requestedDegree == 1);
  assert(reduced.degreeUsed == 0);
  assert(reduced.neighbourCount == 1);
  assert(reduced.rank == 1);
  assert(reduced.conditionEstimate == 1.0);
  assert(reduced.orderReduced);
  assertCoordinate(reduced.target, 1.0, 0.0, 0.0);
  assert(!reducedSurface.points[0].geometry.spatialMetricDerivativeAvailable);

  GridData rankDeficient = makeValidationGrid();
  rankDeficient.coordinates = {
    {0.0, 0.0, 0.0}, {2.0, 0.0, 0.0},
    {0.0, 1.0, 0.0}, {2.0, 1.0, 0.0}
  };
  const std::string rankMessage = interpolationFailure(rankDeficient, degreeOne);
  assert(rankMessage.find("rank deficient") != std::string::npos);
  assert(rankMessage.find("surface point 0") != std::string::npos);

  InterpolatedSurface successfulSurface;
  InterpolateCartesianToSurface(
    makeValidationGrid(),
    {.points = {{.r = 1.0, .rho = 0.0, .sigma = 0.0}}},
    {0.0, 0.0, 0.0},
    successfulSurface,
    degreeOne
  );
  const InterpolationDiagnostics& successful =
    successfulSurface.points[0].interpolationDiagnostics;
  assert(successful.requestedMethod == CartesianToSurfaceMethod::UnweightedDegree1);
  assert(successful.requestedDegree == 1);
  assert(successful.degreeUsed == 1);
  assert(successful.neighbourCount == 4);
  assert(successful.rank == 4);
  assert(std::isfinite(successful.conditionEstimate));
  assert(successful.conditionEstimate >= 1.0);
  assert(!successful.orderReduced);
  assertCoordinate(successful.target, 1.0, 0.0, 0.0);

  std::cout << "Interpolation validation records successful and reduced-order diagnostics"
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
  testAnalyticCoordinateGeometry();
  testAnalyticCoordinateDerivativesAgainstFiniteDifferences();
  testFindNearestGridCells();
  testAcceleratedNeighbourSearchMatchesBruteForce();
  reportNeighbourSearchPerformance();
  testInvalidNearestGridCellInputs();
  testFitLinearPolynomialExactReproduction();
  testFitLinearPolynomialInvalidStencils();
  testWeightedLinearPolynomialReconstruction();
  testWeightedFitFavoursNearbySamples();
  testWeightedLinearPolynomialFailures();
  testFitCompactLinearPolynomial();
  testInterpolateToSurfaceNearestNeighbour();
  testInterpolateCartesianToSurfaceDegree1();
  testInterpolateCartesianToSurfaceContextualFailures();
  testInterpolationValidationAndDiagnostics();
  testInterpolateFromCartesianToShell();
  testInterpolateFromShellToSurface();

  std::cout << "All interpolation tests completed." << std::endl;
  return 0;
}
