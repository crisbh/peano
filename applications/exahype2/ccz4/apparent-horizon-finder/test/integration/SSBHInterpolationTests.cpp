#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <numbers>
#include <string>

#include "../../file-input.h"
#include "../../interpolation.h"

namespace exahype2::ccz4::ahf::tests::integration {

std::string snapshotPath() {
  const std::string sourceFile = __FILE__;
  const std::string sourceDirectory =
    sourceFile.substr(0, sourceFile.find_last_of("/\\"));
  return sourceDirectory
    + "/../../ccz4Results/SSBH/solution-CCZ4.CCZ4Q-1.vtu";
}

bool fileExists(const std::string& path) {
  std::ifstream input(path, std::ios::binary);
  return input.good();
}

SourceVariableMapping sourceMapping() {
  return SourceVariableMapping{
    .conformalMetric = {0, 1, 2, 3, 4, 5},
    .conformalExtrinsicCurvature = {6, 7, 8, 9, 10, 11},
    .traceK = 53,
    .conformalFactor = 54
  };
}

double determinant(const Tensor3& matrix) {
  return matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
       - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
       + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);
}

/** Eigenvalues of a symmetric 3x3 matrix from cyclic Jacobi rotations. */
Vector3 symmetricEigenvalues(Tensor3 matrix) {
  for (int sweep = 0; sweep < 20; ++sweep) {
    const std::array<std::array<std::size_t, 2>, 3> pairs{{
      {{0, 1}}, {{0, 2}}, {{1, 2}}
    }};
    for (const auto& pair : pairs) {
      const std::size_t p = pair[0];
      const std::size_t q = pair[1];
      if (std::abs(matrix[p][q]) < 1.0e-15) {
        continue;
      }
      const double angle =
        0.5 * std::atan2(2.0 * matrix[p][q], matrix[q][q] - matrix[p][p]);
      const double cosine = std::cos(angle);
      const double sine = std::sin(angle);
      const double app = matrix[p][p];
      const double aqq = matrix[q][q];
      const double apq = matrix[p][q];
      matrix[p][p] = cosine * cosine * app - 2.0 * sine * cosine * apq
                   + sine * sine * aqq;
      matrix[q][q] = sine * sine * app + 2.0 * sine * cosine * apq
                   + cosine * cosine * aqq;
      matrix[p][q] = matrix[q][p] = 0.0;
      for (std::size_t k = 0; k < 3; ++k) {
        if (k == p || k == q) {
          continue;
        }
        const double akp = matrix[k][p];
        const double akq = matrix[k][q];
        matrix[k][p] = matrix[p][k] = cosine * akp - sine * akq;
        matrix[k][q] = matrix[q][k] = sine * akp + cosine * akq;
      }
    }
  }
  Vector3 eigenvalues{matrix[0][0], matrix[1][1], matrix[2][2]};
  std::sort(eigenvalues.begin(), eigenvalues.end());
  return eigenvalues;
}

void assertFiniteAndPhysical(const InterpolatedSurface& surface) {
  for (const InterpolatedSurfacePoint& point : surface.points) {
    assert(std::isfinite(point.r));
    assert(std::isfinite(point.rho));
    assert(std::isfinite(point.sigma));
    assert(std::isfinite(point.geometry.traceK));
    for (std::size_t i = 0; i < 3; ++i) {
      assert(std::isfinite(point.surfaceNormal[i]));
      for (std::size_t u = 0; u < 2; ++u) {
        assert(std::isfinite(point.geometry.coordinateGradient[u][i]));
      }
      for (std::size_t j = 0; j < 3; ++j) {
        assert(std::isfinite(point.geometry.T[i][j]));
        assert(std::isfinite(point.surfaceNormalDerivative[i][j]));
        assert(std::isfinite(point.geometry.spatialMetric[i][j]));
        assert(std::isfinite(point.geometry.extrinsicCurvature[i][j]));
        assert(std::abs(point.geometry.spatialMetric[i][j]
                        - point.geometry.spatialMetric[j][i]) < 1.0e-11);
        assert(std::abs(point.geometry.extrinsicCurvature[i][j]
                        - point.geometry.extrinsicCurvature[j][i]) < 1.0e-11);
        for (std::size_t u = 0; u < 2; ++u) {
          assert(std::isfinite(point.geometry.coordinateHessian[u][i][j]));
        }
        for (std::size_t k = 0; k < 3; ++k) {
          assert(std::isfinite(point.geometry.spatialMetricDerivative[k][i][j]));
        }
      }
    }
    assert(determinant(point.geometry.spatialMetric) > 0.0);
    for (double eigenvalue : symmetricEigenvalues(point.geometry.spatialMetric)) {
      assert(eigenvalue > 1.0e-10);
    }
  }
}

InterpolatedSurface interpolate(
  const GridData& gridData,
  const AHGuess& guess,
  CartesianToSurfaceMethod method,
  std::size_t neighbours
) {
  InterpolatedSurface surface;
  InterpolateCartesianToSurface(
    gridData,
    guess,
    {0.0, 0.0, 0.0},
    surface,
    {.method = method,
     .numberOfNeighbours = neighbours,
     .weightSupport = 1.0,
     .allowOrderReduction = false}
  );
  return surface;
}

void testSSBHInterpolation(const GridData& gridData) {
  constexpr double Pi = std::numbers::pi_v<double>;
  // The first five points form an outward radial line; the remainder exercise
  // distinct, non-polar angular directions around the origin.
  const AHGuess guess{.points = {
    {.r = 0.50, .rho = 0.0, .sigma = 0.0},
    {.r = 0.75, .rho = 0.0, .sigma = 0.0},
    {.r = 1.00, .rho = 0.0, .sigma = 0.0},
    {.r = 1.25, .rho = 0.0, .sigma = 0.0},
    {.r = 1.50, .rho = 0.0, .sigma = 0.0},
    {.r = 1.00, .rho = Pi / 4.0, .sigma = 0.0},
    {.r = 1.00, .rho = Pi / 2.0, .sigma = 0.0},
    {.r = 1.00, .rho = Pi, .sigma = 0.0},
    {.r = 1.00, .rho = 0.7, .sigma = 0.25}
  }};

  const InterpolatedSurface mls64 = interpolate(
    gridData, guess, CartesianToSurfaceMethod::WeightedPolynomial, 64);
  const InterpolatedSurface mls128 = interpolate(
    gridData, guess, CartesianToSurfaceMethod::WeightedPolynomial, 128);
  assertFiniteAndPhysical(mls64);
  assertFiniteAndPhysical(mls128);

  // The initial SSBH data have direct phi increasing outward, so tr(gamma)
  // changes smoothly and decreases along this radial line.
  double previousTrace = std::numeric_limits<double>::infinity();
  for (std::size_t point = 0; point < 5; ++point) {
    const Tensor3& gamma = mls128.points[point].geometry.spatialMetric;
    const double trace = gamma[0][0] + gamma[1][1] + gamma[2][2];
    assert(trace <= previousTrace * (1.0 + 1.0e-3));
    if (std::isfinite(previousTrace)) {
      assert(std::abs(trace - previousTrace) / previousTrace < 0.75);
    }
    previousTrace = trace;
  }

  // Doubling the stencil should not substantially change resolved values.
  for (std::size_t point = 0; point < guess.points.size(); ++point) {
    for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 3; ++j) {
        const double first = mls64.points[point].geometry.spatialMetric[i][j];
        const double second = mls128.points[point].geometry.spatialMetric[i][j];
        const double scale = std::max({1.0, std::abs(first), std::abs(second)});
        assert(std::abs(first - second) / scale < 0.20);
      }
    }
  }

  const AHGuess resolvedPoint{.points = {
    {.r = 1.0, .rho = 0.0, .sigma = 0.0}
  }};
  const InterpolatedSurface nearest = interpolate(
    gridData, resolvedPoint, CartesianToSurfaceMethod::NearestNeighbour, 1);
  const InterpolatedSurface mls = interpolate(
    gridData, resolvedPoint, CartesianToSurfaceMethod::WeightedPolynomial, 128);
  assertFiniteAndPhysical(nearest);
  assertFiniteAndPhysical(mls);
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      const double first = nearest.points[0].geometry.spatialMetric[i][j];
      const double second = mls.points[0].geometry.spatialMetric[i][j];
      const double scale = std::max({1.0, std::abs(first), std::abs(second)});
      assert(std::abs(first - second) / scale < 0.25);
    }
  }

  // Independently reconstruct compact phi and check gamma=gammaTilde/phi^2.
  const GridCoordinate target =
    SurfacePointToCartesian(resolvedPoint.points[0], {0.0, 0.0, 0.0});
  const std::vector<std::size_t> indices = FindNearestGridCells(gridData, target, 128);
  const CompactLinearFit compact = FitCompactLinearPolynomial(
    gridData,
    indices,
    target,
    {.method = LinearReconstructionMethod::WeightedMovingLeastSquares,
     .numberOfNeighbours = 128,
     .weightSupport = 1.0}
  );
  const double gammaTildeXX =
    compact.values[static_cast<std::size_t>(CompactField::GammaTildeXX)];
  const double phi = compact.values[static_cast<std::size_t>(CompactField::Phi)];
  const double actualGammaXX = mls.points[0].geometry.spatialMetric[0][0];
  const double directPhiExpected = gammaTildeXX / (phi * phi);
  const double logarithmicPhiResult = gammaTildeXX / std::exp(2.0 * phi);
  assert(std::abs(actualGammaXX - directPhiExpected)
         < 1.0e-10 * std::max(1.0, std::abs(directPhiExpected)));
  assert(std::abs(actualGammaXX - logarithmicPhiResult)
         > 0.1 * std::abs(actualGammaXX));
}

}  // namespace exahype2::ccz4::ahf::tests::integration

int main() {
  using namespace exahype2::ccz4::ahf;
  using namespace exahype2::ccz4::ahf::tests::integration;
  const std::string input = snapshotPath();
  if (!fileExists(input)) {
    std::cout << "Skipping SSBH integration test: snapshot not found at "
              << input << std::endl;
    return 0;
  }

  std::cout << "SSBH interpolation integration test reading " << input << std::endl;
  GridData gridData;
  ReadGridFromFile(input, gridData, sourceMapping());
  testSSBHInterpolation(gridData);
  std::cout << "SSBH interpolation integration test completed." << std::endl;
  return 0;
}
