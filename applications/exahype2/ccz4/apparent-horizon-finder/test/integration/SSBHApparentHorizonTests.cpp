#include "../../apparent-horizon-finder.h"
#include "../../dense-linear-algebra.h"
#include "../../file-input.h"
#include "../../interpolated-surface-residual.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <optional>
#include <string>
#include <vector>

namespace exahype2::ccz4::ahf::tests::integration {
namespace {

using Clock = std::chrono::steady_clock;

std::string DatasetPath() {
  const std::string source = __FILE__;
  const std::string directory = source.substr(0, source.find_last_of("/\\"));
  return directory + "/../../ccz4Results/sample_data_S.vtu";
}

SourceVariableMapping SourceMapping() {
  return {
    .conformalMetric = {0, 1, 2, 3, 4, 5},
    .conformalExtrinsicCurvature = {6, 7, 8, 9, 10, 11},
    .traceK = 53,
    .conformalFactor = 54
  };
}

double Seconds(Clock::time_point start, Clock::time_point finish) {
  return std::chrono::duration<double>(finish - start).count();
}

double Determinant(const Tensor3& matrix) {
  return matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
       - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
       + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);
}

double InfinityNorm(const std::vector<double>& values) {
  double result = 0.0;
  for (double value : values) result = std::max(result, std::abs(value));
  return result;
}

double TwoNorm(const std::vector<double>& values) {
  double result = 0.0;
  for (double value : values) result += value * value;
  return std::sqrt(result);
}

bool AllFinite(const std::vector<double>& values) {
  return std::all_of(values.begin(), values.end(), [](double value) {
    return std::isfinite(value);
  });
}

double MatrixInfinityNorm(const std::vector<std::vector<double>>& matrix) {
  double norm = 0.0;
  for (const auto& row : matrix) {
    double rowSum = 0.0;
    for (double value : row) rowSum += std::abs(value);
    norm = std::max(norm, rowSum);
  }
  return norm;
}

std::size_t EstimateRank(std::vector<std::vector<double>> matrix) {
  const std::size_t rows = matrix.size();
  const std::size_t columns = rows == 0 ? 0 : matrix.front().size();
  const double scale = MatrixInfinityNorm(matrix);
  const double threshold = std::numeric_limits<double>::epsilon()
    * static_cast<double>(std::max(rows, columns)) * scale;
  std::size_t rank = 0;
  for (std::size_t column = 0; column < columns && rank < rows; ++column) {
    std::size_t pivot = rank;
    for (std::size_t row = rank + 1; row < rows; ++row) {
      if (std::abs(matrix[row][column]) > std::abs(matrix[pivot][column])) {
        pivot = row;
      }
    }
    if (std::abs(matrix[pivot][column]) <= threshold) continue;
    std::swap(matrix[pivot], matrix[rank]);
    for (std::size_t row = rank + 1; row < rows; ++row) {
      const double factor = matrix[row][column] / matrix[rank][column];
      for (std::size_t entry = column; entry < columns; ++entry) {
        matrix[row][entry] -= factor * matrix[rank][entry];
      }
    }
    ++rank;
  }
  return rank;
}

double EstimateConditionInfinity(
  const std::vector<std::vector<double>>& matrix
) {
  const std::size_t size = matrix.size();
  std::vector<std::vector<double>> inverse(
    size, std::vector<double>(size, 0.0)
  );
  for (std::size_t column = 0; column < size; ++column) {
    std::vector<double> unit(size, 0.0);
    unit[column] = 1.0;
    const std::vector<double> inverseColumn = SolveLinearSystem(matrix, unit);
    for (std::size_t row = 0; row < size; ++row) {
      inverse[row][column] = inverseColumn[row];
    }
  }
  return MatrixInfinityNorm(matrix) * MatrixInfinityNorm(inverse);
}

void PrintVector(const char* name, const std::vector<double>& values) {
  std::cout << name << " = [";
  for (std::size_t i = 0; i < values.size(); ++i) {
    if (i != 0) std::cout << ", ";
    std::cout << values[i];
  }
  std::cout << "]\n";
}

std::vector<GridCoordinate> SurfaceTargets(
  const SpectralBasis& basis,
  const SynthesizedField& field,
  const ReferencePoint& referencePoint
) {
  std::vector<GridCoordinate> targets;
  targets.reserve(basis.nPoints);
  for (std::size_t point = 0; point < basis.nPoints; ++point) {
    targets.push_back(SurfacePointToCartesian(
      {field.h[point], basis.rho[point], basis.sigma[point]}, referencePoint
    ));
  }
  return targets;
}

using Stencils = std::vector<std::vector<std::size_t>>;

Stencils SelectStencils(
  const CartesianSpatialIndex& index,
  const std::vector<GridCoordinate>& targets,
  std::size_t neighbourCount
) {
  Stencils stencils;
  stencils.reserve(targets.size());
  for (const GridCoordinate& target : targets) {
    stencils.push_back(FindNearestGridCells(index, target, neighbourCount));
  }
  return stencils;
}

std::vector<double> FixedStencilResidual(
  const SpectralBasis& basis,
  const std::vector<double>& coefficients,
  const GridData& grid,
  const ReferencePoint& referencePoint,
  const Stencils& stencils,
  const CartesianToSurfaceOptions& options
) {
  const SynthesizedField field = Synthesize(basis, coefficients);
  const std::vector<GridCoordinate> targets = SurfaceTargets(
    basis, field, referencePoint
  );
  InterpolatedSurface surface;
  surface.points.resize(basis.nPoints);
  const LinearReconstructionOptions fitOptions{
    .method = LinearReconstructionMethod::WeightedMovingLeastSquares,
    .numberOfNeighbours = options.numberOfNeighbours,
    .polynomialDegree = options.polynomialDegree,
    .weightSupport = options.weightSupport,
    .rankToleranceMultiplier = options.rankToleranceMultiplier
  };
  for (std::size_t pointIndex = 0; pointIndex < basis.nPoints; ++pointIndex) {
    const CompactLinearFit fit = FitCompactLinearPolynomial(
      grid, stencils[pointIndex], targets[pointIndex], fitOptions
    );
    std::array<double, NumberOfSymmetricTensorComponents> gammaComponents{};
    std::array<double, NumberOfSymmetricTensorComponents> curvatureComponents{};
    for (std::size_t component = 0;
         component < NumberOfSymmetricTensorComponents; ++component) {
      gammaComponents[component] = fit.values[component];
      curvatureComponents[component] = fit.values[component + 6];
    }
    const Tensor3 gammaTilde = UnpackSymmetricTensor(gammaComponents);
    const Tensor3 aTilde = UnpackSymmetricTensor(curvatureComponents);
    const double traceK = fit.values[static_cast<std::size_t>(CompactField::TraceK)];
    const double phi = fit.values[static_cast<std::size_t>(CompactField::Phi)];
    Tensor3Derivative gammaTildeDerivative{};
    for (std::size_t direction = 0; direction < 3; ++direction) {
      std::array<double, NumberOfSymmetricTensorComponents> components{};
      for (std::size_t component = 0;
           component < NumberOfSymmetricTensorComponents; ++component) {
        components[component] = fit.gradients[component][direction];
      }
      gammaTildeDerivative[direction] = UnpackSymmetricTensor(components);
    }
    InterpolatedSurfacePoint& point = surface.points[pointIndex];
    point.r = field.h[pointIndex];
    point.rho = basis.rho[pointIndex];
    point.sigma = basis.sigma[pointIndex];
    point.geometry.spatialMetric = ComputePhysicalSpatialMetric(gammaTilde, phi);
    point.geometry.extrinsicCurvature = ComputePhysicalExtrinsicCurvature(
      aTilde, traceK, gammaTilde, phi
    );
    point.geometry.traceK = traceK;
    point.geometry.spatialMetricDerivative = ComputePhysicalSpatialMetricDerivative(
      gammaTilde,
      gammaTildeDerivative,
      phi,
      fit.gradients[static_cast<std::size_t>(CompactField::Phi)]
    );
    point.geometry.spatialMetricDerivativeAvailable = true;
    const AnalyticCoordinateGeometry coordinate =
      ComputeAnalyticCoordinateGeometry(targets[pointIndex], referencePoint);
    point.geometry.T = coordinate.T;
    point.geometry.coordinateGradient = coordinate.gradient;
    point.geometry.coordinateHessian = coordinate.hessian;
  }
  PopulateSurfaceShapeGeometry(surface, field);
  std::vector<std::vector<double>> unused;
  std::vector<double> residual;
  AssembleGlobalMatrix(surface, basis, unused, residual);
  return residual;
}

std::vector<double> DifferenceScaled(
  const std::vector<double>& first,
  const std::vector<double>& second,
  double scale
) {
  std::vector<double> result(first.size());
  for (std::size_t i = 0; i < first.size(); ++i) {
    result[i] = scale * (first[i] - second[i]);
  }
  return result;
}

std::vector<double> Difference(
  const std::vector<double>& first,
  const std::vector<double>& second
) {
  return DifferenceScaled(first, second, 1.0);
}

struct RadiusResidualEvaluation {
  InterpolatedSurface surface;
  std::vector<double> residual;
  std::vector<std::size_t> supportCounts;
  std::vector<std::vector<std::size_t>> supports;
  std::vector<double> fitConditions;
};

RadiusResidualEvaluation EvaluateRadiusResidual(
  const SpectralBasis& basis,
  const std::vector<double>& coefficients,
  const CartesianPointInterpolator& interpolator,
  const CartesianSpatialIndex& index,
  const ReferencePoint& referencePoint,
  double supportRadius
) {
  const SynthesizedField field = Synthesize(basis, coefficients);
  const std::vector<GridCoordinate> targets = SurfaceTargets(
    basis, field, referencePoint
  );
  RadiusResidualEvaluation evaluation;
  evaluation.surface = BuildInterpolatedSurface(
    basis, coefficients, interpolator, referencePoint
  );
  for (std::size_t pointIndex = 0; pointIndex < basis.nPoints; ++pointIndex) {
    const auto support = FindGridCellsWithinRadius(
      index, targets[pointIndex], supportRadius
    );
    evaluation.supportCounts.push_back(support.size());
    evaluation.supports.push_back(support);
    evaluation.fitConditions.push_back(
      evaluation.surface.points[pointIndex].interpolationDiagnostics.conditionEstimate
    );
  }
  std::vector<std::vector<double>> unused;
  AssembleGlobalMatrix(
    evaluation.surface, basis, unused, evaluation.residual
  );
  return evaluation;
}

double Distance(const GridCoordinate& first, const GridCoordinate& second) {
  return std::hypot(
    first.x - second.x,
    first.y - second.y,
    first.z - second.z
  );
}

double SquaredDistance(
  const GridCoordinate& first,
  const GridCoordinate& second
) {
  const double dx = first.x - second.x;
  const double dy = first.y - second.y;
  const double dz = first.z - second.z;
  return dx * dx + dy * dy + dz * dz;
}

struct Range {
  double minimum;
  double mean;
  double maximum;
};

Range Summarize(const std::vector<double>& values) {
  assert(!values.empty());
  return {
    *std::min_element(values.begin(), values.end()),
    std::accumulate(values.begin(), values.end(), 0.0) / values.size(),
    *std::max_element(values.begin(), values.end())
  };
}

void RequireFinite(double value) {
  assert(std::isfinite(value));
}

void ValidateSurfaceGeometry(
  const InterpolatedSurface& surface,
  std::vector<double>& determinants,
  std::vector<double>& expansion
) {
  determinants.clear();
  expansion.clear();
  determinants.reserve(surface.numberOfPoints());
  expansion.reserve(surface.numberOfPoints());
  constexpr double symmetryTolerance = 1.0e-11;

  for (std::size_t pointIndex = 0; pointIndex < surface.numberOfPoints(); ++pointIndex) {
    const InterpolatedSurfacePoint& point = surface.points[pointIndex];
    assert(point.geometry.spatialMetricDerivativeAvailable);
    RequireFinite(point.r);
    RequireFinite(point.rho);
    RequireFinite(point.sigma);
    RequireFinite(point.geometry.traceK);
    for (std::size_t i = 0; i < 3; ++i) {
      RequireFinite(point.surfaceNormal[i]);
      for (std::size_t u = 0; u < 2; ++u) {
        RequireFinite(point.geometry.coordinateGradient[u][i]);
      }
      for (std::size_t j = 0; j < 3; ++j) {
        RequireFinite(point.geometry.spatialMetric[i][j]);
        RequireFinite(point.geometry.extrinsicCurvature[i][j]);
        RequireFinite(point.surfaceNormalDerivative[i][j]);
        RequireFinite(point.geometry.T[i][j]);
        assert(std::abs(point.geometry.spatialMetric[i][j]
                        - point.geometry.spatialMetric[j][i])
               <= symmetryTolerance * std::max({
                    1.0,
                    std::abs(point.geometry.spatialMetric[i][j]),
                    std::abs(point.geometry.spatialMetric[j][i])
                  }));
        for (std::size_t k = 0; k < 3; ++k) {
          RequireFinite(point.geometry.spatialMetricDerivative[k][i][j]);
        }
        for (std::size_t u = 0; u < 2; ++u) {
          RequireFinite(point.geometry.coordinateHessian[u][i][j]);
        }
      }
    }
    const double determinant = Determinant(point.geometry.spatialMetric);
    RequireFinite(determinant);
    assert(determinant > 0.0);
    determinants.push_back(determinant);

    const double theta = computeExpansion(point);
    RequireFinite(theta);
    expansion.push_back(theta);
  }
}

void RequireInsideDomain(
  const SpectralBasis& basis,
  const SynthesizedField& field,
  const ReferencePoint& referencePoint,
  const GridCoordinate& minimum,
  const GridCoordinate& maximum
) {
  assert(field.h.size() == basis.nPoints);
  for (std::size_t point = 0; point < basis.nPoints; ++point) {
    const GridCoordinate target = SurfacePointToCartesian(
      {field.h[point], basis.rho[point], basis.sigma[point]}, referencePoint
    );
    assert(target.x >= minimum.x && target.x <= maximum.x);
    assert(target.y >= minimum.y && target.y <= maximum.y);
    assert(target.z >= minimum.z && target.z <= maximum.z);
  }
}

void PrintRange(const char* name, const Range& range) {
  std::cout << name << " min/mean/max = " << range.minimum << " / "
            << range.mean << " / " << range.maximum << '\n';
}

void TestSampleDataSingleBlackHolePipeline() {
  const std::string path = DatasetPath();
  std::ifstream probe(path, std::ios::binary);
  assert(probe.good());

  const auto readStart = Clock::now();
  GridData grid;
  ReadGridFromFile(path, grid, SourceMapping()); // The only VTU read in this test.
  const auto readFinish = Clock::now();

  constexpr std::size_t expectedCells = 3088448;
  static_assert(GridData::numberOfVariables == 14);
  assert(grid.numberOfPoints() == expectedCells);
  assert(grid.isConsistent());
  assert(grid.values.size() == expectedCells * GridData::numberOfVariables);

  GridCoordinate domainMinimum = grid.coordinates.front();
  GridCoordinate domainMaximum = grid.coordinates.front();
  for (const GridCoordinate& coordinate : grid.coordinates) {
    RequireFinite(coordinate.x);
    RequireFinite(coordinate.y);
    RequireFinite(coordinate.z);
    domainMinimum.x = std::min(domainMinimum.x, coordinate.x);
    domainMinimum.y = std::min(domainMinimum.y, coordinate.y);
    domainMinimum.z = std::min(domainMinimum.z, coordinate.z);
    domainMaximum.x = std::max(domainMaximum.x, coordinate.x);
    domainMaximum.y = std::max(domainMaximum.y, coordinate.y);
    domainMaximum.z = std::max(domainMaximum.z, coordinate.z);
  }

  // The data numerically resemble an M~1 isotropic puncture, whose candidate
  // horizon is near r=0.5. Estimate local centre spacing there without building
  // a second spatial index.
  const std::array<GridCoordinate, 6> resolutionProbes{{
    {0.5, 0.0, 0.0}, {-0.5, 0.0, 0.0},
    {0.0, 0.5, 0.0}, {0.0, -0.5, 0.0},
    {0.0, 0.0, 0.5}, {0.0, 0.0, -0.5}
  }};
  double localSpacingMinimum = std::numeric_limits<double>::infinity();
  double localSpacingMaximum = 0.0;
  for (const GridCoordinate& target : resolutionProbes) {
    const std::vector<std::size_t> nearest =
      FindNearestGridCellsBruteForce(grid, target, 2);
    const double spacing = Distance(
      grid.coordinates[nearest[0]], grid.coordinates[nearest[1]]
    );
    localSpacingMinimum = std::min(localSpacingMinimum, spacing);
    localSpacingMaximum = std::max(localSpacingMaximum, spacing);
  }

  const CartesianToSurfaceOptions interpolationOptions{
    .method = CartesianToSurfaceMethod::WeightedPolynomial,
    .numberOfNeighbours = 64,
    .polynomialDegree = 1,
    .weightSupport = 1.0,
    .rankToleranceMultiplier = 64.0,
    .allowOrderReduction = false,
    .domainRelativeTolerance = 0.05
  };
  const auto indexStart = Clock::now();
  const CartesianPointInterpolator interpolator(grid, interpolationOptions);
  const auto indexFinish = Clock::now();

  constexpr int lmax = 2;
  constexpr double initialRadius = 0.8; // Exploratory: numerical data are consistent with M~1 isotropic data.
  constexpr ReferencePoint referencePoint{0.0, 0.0, 0.0};
  const SpectralBasis basis = BuildSpectralBasis(lmax);
  const SpectralNewtonOptions newtonOptions{
    .tolerance = 1.0e-8,
    .massScale = 1.0,
    .maxIterations = 12,
    .finiteDifferenceStep = 1.0e-6,
    .maxLineSearchSteps = 20,
    .maxRadiusFactor = 6.0
  };

  const std::vector<double> initialCoefficients =
    MakeSphericalInitialCoefficients(basis, initialRadius);
  const SynthesizedField initialField = Synthesize(basis, initialCoefficients);
  const Range initialRadiusRange = Summarize(initialField.h);
  for (double radius : initialField.h) {
    RequireFinite(radius);
    assert(radius > 0.0);
  }
  RequireInsideDomain(
    basis, initialField, referencePoint, domainMinimum, domainMaximum
  );

  const auto initialResidualStart = Clock::now();
  const InterpolatedSurface initialSurface = BuildInterpolatedSurface(
    basis, initialCoefficients, interpolator, referencePoint
  );
  assert(initialSurface.numberOfPoints() == basis.nPoints);
  std::vector<double> initialDeterminants;
  std::vector<double> initialExpansion;
  ValidateSurfaceGeometry(initialSurface, initialDeterminants, initialExpansion);
  const SpectralResidual residual = detail::MakeInterpolatedSurfaceResidual(
    basis, interpolator, referencePoint
  );
  const std::vector<double> initialResidual = residual(initialCoefficients);
  assert(initialResidual.size() == basis.nModes);
  for (double value : initialResidual) RequireFinite(value);
  const double initialResidualNorm = InfinityNorm(initialResidual);
  const double initialResidualTwoNorm = TwoNorm(initialResidual);
  const auto initialResidualFinish = Clock::now();

  // Reproduce the first Newton iteration exactly, for diagnostics only.
  const std::size_t modeCount = initialCoefficients.size();
  std::vector<std::vector<double>> jacobian(
    modeCount, std::vector<double>(modeCount, 0.0)
  );
  bool finiteDifferenceResidualsFinite = true;
  for (std::size_t mode = 0; mode < modeCount; ++mode) {
    const double epsilon = newtonOptions.finiteDifferenceStep
      * std::max(1.0, std::abs(initialCoefficients[mode]));
    std::vector<double> plus = initialCoefficients;
    std::vector<double> minus = initialCoefficients;
    plus[mode] += epsilon;
    minus[mode] -= epsilon;
    const std::vector<double> residualPlus = residual(plus);
    const std::vector<double> residualMinus = residual(minus);
    finiteDifferenceResidualsFinite = finiteDifferenceResidualsFinite
      && AllFinite(residualPlus) && AllFinite(residualMinus);
    for (std::size_t equation = 0; equation < modeCount; ++equation) {
      jacobian[equation][mode] =
        (residualPlus[equation] - residualMinus[equation]) / (2.0 * epsilon);
    }
  }
  bool jacobianFinite = true;
  for (const auto& row : jacobian) {
    jacobianFinite = jacobianFinite && AllFinite(row);
  }
  const std::size_t jacobianRank = EstimateRank(jacobian);
  const double jacobianCondition = EstimateConditionInfinity(jacobian);
  std::vector<double> negativeResidual(initialResidual.size());
  std::transform(
    initialResidual.begin(), initialResidual.end(), negativeResidual.begin(),
    [](double value) { return -value; }
  );
  const std::vector<double> newtonDelta = SolveLinearSystem(
    jacobian, negativeResidual
  );
  const double deltaTwoNorm = TwoNorm(newtonDelta);
  const double deltaInfinityNorm = InfinityNorm(newtonDelta);

  std::vector<double> jacobianTimesDelta(modeCount, 0.0);
  for (std::size_t row = 0; row < modeCount; ++row) {
    for (std::size_t column = 0; column < modeCount; ++column) {
      jacobianTimesDelta[row] += jacobian[row][column] * newtonDelta[column];
    }
  }

  const std::vector<double> positiveScales{
    1.0e-2, 5.0e-3, 2.0e-3, 1.0e-3,
    5.0e-4, 2.0e-4, 1.0e-4, 5.0e-5,
    2.0e-5, 1.0e-5, 5.0e-6, 2.0e-6,
    1.0e-6, 5.0e-7, 2.0e-7, 1.0e-7
  };
  const CartesianSpatialIndex diagnosticIndex(grid);
  const double interpolationDomainDiagonal = std::hypot(
    domainMaximum.x - domainMinimum.x,
    domainMaximum.y - domainMinimum.y,
    domainMaximum.z - domainMinimum.z
  );
  const auto targetsInsideInterpolationDomain = [&](const auto& targets) {
    for (const GridCoordinate& target : targets) {
      const double outsideX = std::max(
        {domainMinimum.x - target.x, 0.0, target.x - domainMaximum.x});
      const double outsideY = std::max(
        {domainMinimum.y - target.y, 0.0, target.y - domainMaximum.y});
      const double outsideZ = std::max(
        {domainMinimum.z - target.z, 0.0, target.z - domainMaximum.z});
      if (std::hypot(outsideX, outsideY, outsideZ)
          > interpolationOptions.domainRelativeTolerance
            * interpolationDomainDiagonal) {
        return false;
      }
    }
    return true;
  };

  // Basic degree-one reproduction check through the production Wendland-C2 fit.
  {
    GridData affineGrid;
    for (int x = -2; x <= 2; ++x) {
      for (int y = -2; y <= 2; ++y) {
        for (int z = -2; z <= 2; ++z) {
          affineGrid.coordinates.push_back({
            static_cast<double>(x), static_cast<double>(y), static_cast<double>(z)
          });
          for (std::size_t field = 0; field < NumberOfCompactFields; ++field) {
            affineGrid.values.push_back(
              2.0 + 0.1 * field + 3.0 * x - 2.0 * y + 0.5 * z
            );
          }
        }
      }
    }
    const GridCoordinate target{0.1, 0.2, -0.1};
    const CartesianSpatialIndex affineIndex(affineGrid);
    const auto support = FindGridCellsWithinRadius(affineIndex, target, 2.5);
    const CompactLinearFit fit = FitCompactLinearPolynomial(
      affineGrid,
      support,
      target,
      {
        .method = LinearReconstructionMethod::RadiusWeightedWendlandC2,
        .numberOfNeighbours = support.size(),
        .polynomialDegree = 1,
        .supportRadius = 2.5,
        .rankToleranceMultiplier = 64.0
      }
    );
    for (std::size_t field = 0; field < NumberOfCompactFields; ++field) {
      const double expected = 2.0 + 0.1 * field
        + 3.0 * target.x - 2.0 * target.y + 0.5 * target.z;
      assert(std::abs(fit.values[field] - expected) < 1.0e-12);
      assert(std::abs(fit.gradients[field][0] - 3.0) < 1.0e-12);
      assert(std::abs(fit.gradients[field][1] + 2.0) < 1.0e-12);
      assert(std::abs(fit.gradients[field][2] - 0.5) < 1.0e-12);
    }
    assert(fit.rank == 4);
    assert(std::isfinite(fit.conditionEstimate));
  }

  constexpr double radiusSupport = 0.1;
  const CartesianToSurfaceOptions radiusInterpolationOptions{
    .method = CartesianToSurfaceMethod::RadiusWeightedPolynomial,
    .numberOfNeighbours = 0,
    .polynomialDegree = 1,
    .weightSupport = 1.0,
    .supportRadius = radiusSupport,
    .rankToleranceMultiplier = interpolationOptions.rankToleranceMultiplier,
    .allowOrderReduction = false,
    .domainRelativeTolerance = interpolationOptions.domainRelativeTolerance
  };
  const CartesianPointInterpolator radiusInterpolator(
    grid, radiusInterpolationOptions
  );
  const RadiusResidualEvaluation radiusInitial = EvaluateRadiusResidual(
    basis,
    initialCoefficients,
    radiusInterpolator,
    diagnosticIndex,
    referencePoint,
    radiusSupport
  );
  const auto radiusResidualAt = [&](const std::vector<double>& coefficients) {
    return EvaluateRadiusResidual(
      basis,
      coefficients,
      radiusInterpolator,
      diagnosticIndex,
      referencePoint,
      radiusSupport
    );
  };

  std::vector<std::vector<double>> radiusJacobian(
    modeCount, std::vector<double>(modeCount, 0.0)
  );
  for (std::size_t mode = 0; mode < modeCount; ++mode) {
    const double epsilon = newtonOptions.finiteDifferenceStep
      * std::max(1.0, std::abs(initialCoefficients[mode]));
    std::vector<double> plus = initialCoefficients;
    std::vector<double> minus = initialCoefficients;
    plus[mode] += epsilon;
    minus[mode] -= epsilon;
    const std::vector<double> plusResidual = radiusResidualAt(plus).residual;
    const std::vector<double> minusResidual = radiusResidualAt(minus).residual;
    for (std::size_t equation = 0; equation < modeCount; ++equation) {
      radiusJacobian[equation][mode] =
        (plusResidual[equation] - minusResidual[equation]) / (2.0 * epsilon);
    }
  }
  std::vector<double> negativeRadiusResidual(radiusInitial.residual.size());
  std::transform(
    radiusInitial.residual.begin(),
    radiusInitial.residual.end(),
    negativeRadiusResidual.begin(),
    [](double value) { return -value; }
  );
  const std::vector<double> radiusDelta = SolveLinearSystem(
    radiusJacobian, negativeRadiusResidual
  );
  std::vector<double> radiusJacobianTimesDelta(modeCount, 0.0);
  for (std::size_t row = 0; row < modeCount; ++row) {
    for (std::size_t column = 0; column < modeCount; ++column) {
      radiusJacobianTimesDelta[row] += radiusJacobian[row][column]
        * radiusDelta[column];
    }
  }

  struct RadiusMotionSample {
    double t;
    std::vector<std::size_t> supportCounts;
    std::size_t changedSupportPoints;
    double twoNorm;
    double infinityNorm;
  };
  const std::vector<double> radiusMotionScales{
    5.0e-6, 2.5e-6, 2.0e-6, 1.75e-6, 1.5e-6,
    1.25e-6, 1.0e-6, 5.0e-7, 2.0e-7, 1.0e-7
  };
  std::vector<RadiusMotionSample> radiusMotionSamples;
  for (double magnitude : radiusMotionScales) {
    for (double sign : {1.0, -1.0}) {
      std::vector<double> coefficients(modeCount);
      for (std::size_t mode = 0; mode < modeCount; ++mode) {
        coefficients[mode] = initialCoefficients[mode]
          + sign * magnitude * newtonDelta[mode];
      }
      const RadiusResidualEvaluation evaluation = radiusResidualAt(coefficients);
      std::size_t changedPoints = 0;
      for (std::size_t point = 0; point < basis.nPoints; ++point) {
        if (evaluation.supports[point] != radiusInitial.supports[point]) {
          ++changedPoints;
        }
      }
      radiusMotionSamples.push_back({
        sign * magnitude,
        evaluation.supportCounts,
        changedPoints,
        TwoNorm(evaluation.residual),
        InfinityNorm(evaluation.residual)
      });
    }
  }

  struct RadiusDerivativeSample {
    double t;
    double plusError;
    double minusError;
    double centralError;
    double centralRelative;
    std::size_t plusChangedPoints;
    std::size_t minusChangedPoints;
  };
  std::vector<RadiusDerivativeSample> radiusDerivativeSamples;
  for (double t : radiusMotionScales) {
    std::vector<double> plus(modeCount);
    std::vector<double> minus(modeCount);
    for (std::size_t mode = 0; mode < modeCount; ++mode) {
      plus[mode] = initialCoefficients[mode] + t * radiusDelta[mode];
      minus[mode] = initialCoefficients[mode] - t * radiusDelta[mode];
    }
    const RadiusResidualEvaluation plusEvaluation = radiusResidualAt(plus);
    const RadiusResidualEvaluation minusEvaluation = radiusResidualAt(minus);
    const std::vector<double> dPlus = DifferenceScaled(
      plusEvaluation.residual, radiusInitial.residual, 1.0 / t
    );
    const std::vector<double> dMinus = DifferenceScaled(
      radiusInitial.residual, minusEvaluation.residual, 1.0 / t
    );
    const std::vector<double> dCentral = DifferenceScaled(
      plusEvaluation.residual, minusEvaluation.residual, 1.0 / (2.0 * t)
    );
    std::size_t plusChanged = 0;
    std::size_t minusChanged = 0;
    for (std::size_t point = 0; point < basis.nPoints; ++point) {
      plusChanged += plusEvaluation.supports[point] != radiusInitial.supports[point];
      minusChanged += minusEvaluation.supports[point] != radiusInitial.supports[point];
    }
    const double centralError = TwoNorm(
      Difference(dCentral, radiusJacobianTimesDelta)
    );
    radiusDerivativeSamples.push_back({
      t,
      TwoNorm(Difference(dPlus, radiusJacobianTimesDelta)),
      TwoNorm(Difference(dMinus, radiusJacobianTimesDelta)),
      centralError,
      centralError / TwoNorm(radiusJacobianTimesDelta),
      plusChanged,
      minusChanged
    });
  }
  assert(std::any_of(
    radiusMotionSamples.begin(), radiusMotionSamples.end(),
    [](const RadiusMotionSample& sample) {
      return sample.changedSupportPoints > 0;
    }
  ));
  {
    std::vector<double> atTwo(modeCount);
    std::vector<double> atOnePointFive(modeCount);
    for (std::size_t mode = 0; mode < modeCount; ++mode) {
      atTwo[mode] = initialCoefficients[mode] + 2.0e-6 * newtonDelta[mode];
      atOnePointFive[mode] = initialCoefficients[mode]
        + 1.5e-6 * newtonDelta[mode];
    }
    const double transitionChange = TwoNorm(Difference(
      radiusResidualAt(atTwo).residual,
      radiusResidualAt(atOnePointFive).residual
    ));
    // The hard-64 path jumps by approximately 1.9e-3 here.  The compact
    // production path varies continuously across the same membership event.
    assert(transitionChange < 1.0e-4);
  }
  assert(!radiusDerivativeSamples.empty());
  assert(radiusDerivativeSamples.back().centralRelative < 1.0e-3);

  struct RadiusLineTrial {
    double lambda;
    Range radius;
    bool valid;
    double twoNorm;
    double infinityNorm;
    bool accepted;
    std::vector<double> coefficients;
  };
  std::vector<RadiusLineTrial> radiusLineTrials;
  std::optional<std::size_t> radiusAccepted;
  double radiusLambda = 1.0;
  for (int step = 0; step < newtonOptions.maxLineSearchSteps; ++step) {
    std::vector<double> coefficients(modeCount);
    for (std::size_t mode = 0; mode < modeCount; ++mode) {
      coefficients[mode] = initialCoefficients[mode]
        + radiusLambda * radiusDelta[mode];
    }
    const SynthesizedField field = Synthesize(basis, coefficients);
    RadiusLineTrial trial{
      radiusLambda,
      Summarize(field.h),
      false,
      std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::quiet_NaN(),
      false,
      coefficients
    };
    const std::vector<GridCoordinate> targets = SurfaceTargets(
      basis, field, referencePoint
    );
    if (std::all_of(field.h.begin(), field.h.end(), [](double h) {
          return std::isfinite(h) && h > 0.0;
        }) && targetsInsideInterpolationDomain(targets)) {
      try {
        const RadiusResidualEvaluation evaluation = radiusResidualAt(coefficients);
        trial.valid = AllFinite(evaluation.residual);
        if (trial.valid) {
          trial.twoNorm = TwoNorm(evaluation.residual);
          trial.infinityNorm = InfinityNorm(evaluation.residual);
          if (!radiusAccepted
              && trial.twoNorm < TwoNorm(radiusInitial.residual)) {
            trial.accepted = true;
            radiusAccepted = radiusLineTrials.size();
          }
        }
      } catch (const std::exception&) {
      }
    }
    radiusLineTrials.push_back(std::move(trial));
    radiusLambda *= 0.5;
  }

  struct RadiusShortIteration {
    int iteration;
    double beforeNorm;
    double acceptedLambda;
    double afterNorm;
    Range radius;
    std::size_t rank;
    double condition;
  };
  std::vector<RadiusShortIteration> radiusShortIterations;
  if (radiusAccepted) {
    std::vector<double> center =
      radiusLineTrials[*radiusAccepted].coefficients;
    for (int iteration = 1; iteration <= 2; ++iteration) {
      const RadiusResidualEvaluation centerEvaluation = radiusResidualAt(center);
      std::vector<std::vector<double>> iterationJacobian(
        modeCount, std::vector<double>(modeCount, 0.0)
      );
      for (std::size_t mode = 0; mode < modeCount; ++mode) {
        const double epsilon = newtonOptions.finiteDifferenceStep
          * std::max(1.0, std::abs(center[mode]));
        std::vector<double> plus = center;
        std::vector<double> minus = center;
        plus[mode] += epsilon;
        minus[mode] -= epsilon;
        const auto plusResidual = radiusResidualAt(plus).residual;
        const auto minusResidual = radiusResidualAt(minus).residual;
        for (std::size_t equation = 0; equation < modeCount; ++equation) {
          iterationJacobian[equation][mode] =
            (plusResidual[equation] - minusResidual[equation])
            / (2.0 * epsilon);
        }
      }
      std::vector<double> negative(centerEvaluation.residual.size());
      std::transform(
        centerEvaluation.residual.begin(), centerEvaluation.residual.end(),
        negative.begin(), [](double value) { return -value; }
      );
      const std::vector<double> delta = SolveLinearSystem(
        iterationJacobian, negative
      );
      double lambda = 1.0;
      bool accepted = false;
      for (int step = 0; step < newtonOptions.maxLineSearchSteps; ++step) {
        std::vector<double> trial(center.size());
        for (std::size_t mode = 0; mode < center.size(); ++mode) {
          trial[mode] = center[mode] + lambda * delta[mode];
        }
        const SynthesizedField trialField = Synthesize(basis, trial);
        const auto targets = SurfaceTargets(basis, trialField, referencePoint);
        if (std::all_of(trialField.h.begin(), trialField.h.end(), [](double h) {
              return std::isfinite(h) && h > 0.0;
            }) && targetsInsideInterpolationDomain(targets)) {
          const RadiusResidualEvaluation trialEvaluation = radiusResidualAt(trial);
          if (AllFinite(trialEvaluation.residual)
              && TwoNorm(trialEvaluation.residual)
                < TwoNorm(centerEvaluation.residual)) {
            radiusShortIterations.push_back({
              iteration,
              TwoNorm(centerEvaluation.residual),
              lambda,
              TwoNorm(trialEvaluation.residual),
              Summarize(trialField.h),
              EstimateRank(iterationJacobian),
              EstimateConditionInfinity(iterationJacobian)
            });
            center = std::move(trial);
            accepted = true;
            break;
          }
        }
        lambda *= 0.5;
      }
      if (!accepted) break;
    }
  }
  const std::vector<GridCoordinate> initialTargets = SurfaceTargets(
    basis, initialField, referencePoint
  );
  const Stencils initialStencils = SelectStencils(
    diagnosticIndex, initialTargets, interpolationOptions.numberOfNeighbours
  );

  struct DirectionalSample {
    double t;
    Range radius;
    bool valid;
    std::vector<double> residual;
    double residualTwoNorm;
    double residualInfinityNorm;
    std::string failure;
    std::size_t affectedPoints;
    double meanChangedNeighbours;
    std::size_t maximumChangedNeighbours;
    std::size_t nearestChangedPoints;
    double maximumDisplacement;
    std::vector<double> fixedResidual;
  };
  std::vector<DirectionalSample> directionalSamples;
  std::vector<double> smallestStencilChange(
    basis.nPoints, std::numeric_limits<double>::infinity()
  );
  for (double magnitude : positiveScales) {
    for (double sign : {1.0, -1.0}) {
      const double t = sign * magnitude;
      std::vector<double> coefficients(modeCount);
      for (std::size_t mode = 0; mode < modeCount; ++mode) {
        coefficients[mode] = initialCoefficients[mode] + t * newtonDelta[mode];
      }
      const SynthesizedField field = Synthesize(basis, coefficients);
      const std::vector<GridCoordinate> targets = SurfaceTargets(
        basis, field, referencePoint
      );
      const Stencils stencils = SelectStencils(
        diagnosticIndex, targets, interpolationOptions.numberOfNeighbours
      );
      std::size_t affectedPoints = 0;
      std::size_t totalChanged = 0;
      std::size_t maximumChanged = 0;
      std::size_t nearestChanged = 0;
      double maximumDisplacement = 0.0;
      for (std::size_t point = 0; point < basis.nPoints; ++point) {
        std::size_t changed = 0;
        for (std::size_t cell : initialStencils[point]) {
          if (std::find(stencils[point].begin(), stencils[point].end(), cell)
              == stencils[point].end()) {
            ++changed;
          }
        }
        if (changed != 0) {
          ++affectedPoints;
          smallestStencilChange[point] = std::min(
            smallestStencilChange[point], magnitude
          );
        }
        totalChanged += changed;
        maximumChanged = std::max(maximumChanged, changed);
        nearestChanged += stencils[point].front() != initialStencils[point].front();
        maximumDisplacement = std::max(
          maximumDisplacement, Distance(targets[point], initialTargets[point])
        );
      }
      DirectionalSample sample{
        t,
        Summarize(field.h),
        false,
        {},
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN(),
        "",
        affectedPoints,
        static_cast<double>(totalChanged) / basis.nPoints,
        maximumChanged,
        nearestChanged,
        maximumDisplacement,
        {}
      };
      try {
        (void)BuildInterpolatedSurface(
          basis, coefficients, interpolator, referencePoint
        );
        sample.residual = residual(coefficients);
        if (!AllFinite(sample.residual)) {
          sample.failure = "non-finite residual";
        } else {
          sample.valid = true;
          sample.residualTwoNorm = TwoNorm(sample.residual);
          sample.residualInfinityNorm = InfinityNorm(sample.residual);
          sample.fixedResidual = FixedStencilResidual(
            basis,
            coefficients,
            grid,
            referencePoint,
            initialStencils,
            interpolationOptions
          );
        }
      } catch (const std::exception& error) {
        sample.failure = error.what();
      }
      directionalSamples.push_back(std::move(sample));
    }
  }

  // Diagnose the 64/65 boundary. Ranks below are one-based: rank 64 is
  // neighbours[63]. Point zero is an unaffected control at these scales.
  const std::array<std::size_t, 5> boundaryPoints{{2, 3, 12, 13, 0}};
  const std::vector<double> boundaryScales{
    0.0, 5.0e-7, 1.0e-6, 1.25e-6, 1.5e-6,
    1.75e-6, 2.0e-6, 2.5e-6, 5.0e-6
  };
  struct BoundaryRecord {
    double t;
    GridCoordinate target;
    std::vector<std::size_t> neighbours;
  };
  struct BoundaryPointDiagnostic {
    std::size_t point;
    std::vector<BoundaryRecord> records;
    std::size_t leavingCell = std::numeric_limits<std::size_t>::max();
    std::size_t enteringCell = std::numeric_limits<std::size_t>::max();
    double crossingT = std::numeric_limits<double>::quiet_NaN();
    double crossingDisplacement = std::numeric_limits<double>::quiet_NaN();
    double finiteDifferenceDisplacementMinimum =
      std::numeric_limits<double>::infinity();
    double finiteDifferenceDisplacementMaximum = 0.0;
  };
  std::vector<BoundaryPointDiagnostic> boundaryDiagnostics;
  for (std::size_t point : boundaryPoints) {
    BoundaryPointDiagnostic diagnostic{.point = point};
    for (double t : boundaryScales) {
      std::vector<double> coefficients(modeCount);
      for (std::size_t mode = 0; mode < modeCount; ++mode) {
        coefficients[mode] = initialCoefficients[mode] + t * newtonDelta[mode];
      }
      const SynthesizedField field = Synthesize(basis, coefficients);
      const GridCoordinate target = SurfaceTargets(
        basis, field, referencePoint
      )[point];
      diagnostic.records.push_back({
        t, target, FindNearestGridCells(diagnosticIndex, target, 70)
      });
    }
    const auto& initial70 = diagnostic.records.front().neighbours;
    const auto& final70 = diagnostic.records.back().neighbours;
    for (std::size_t rank = 0; rank < 64; ++rank) {
      if (std::find(final70.begin(), final70.begin() + 64, initial70[rank])
          == final70.begin() + 64) {
        diagnostic.leavingCell = initial70[rank];
      }
      if (std::find(initial70.begin(), initial70.begin() + 64, final70[rank])
          == initial70.begin() + 64) {
        diagnostic.enteringCell = final70[rank];
      }
    }
    if (diagnostic.leavingCell != std::numeric_limits<std::size_t>::max()
        && diagnostic.enteringCell != std::numeric_limits<std::size_t>::max()) {
      const GridCoordinate& target0 = diagnostic.records.front().target;
      const GridCoordinate& target5 = diagnostic.records.back().target;
      const GridCoordinate& leaving = grid.coordinates[diagnostic.leavingCell];
      const GridCoordinate& entering = grid.coordinates[diagnostic.enteringCell];
      const double difference0 = SquaredDistance(target0, leaving)
        - SquaredDistance(target0, entering);
      const double difference5 = SquaredDistance(target5, leaving)
        - SquaredDistance(target5, entering);
      diagnostic.crossingT = boundaryScales.back()
        * (-difference0) / (difference5 - difference0);
      const double fraction = diagnostic.crossingT / boundaryScales.back();
      const GridCoordinate crossingTarget{
        target0.x + fraction * (target5.x - target0.x),
        target0.y + fraction * (target5.y - target0.y),
        target0.z + fraction * (target5.z - target0.z)
      };
      diagnostic.crossingDisplacement = Distance(target0, crossingTarget);
    }
    for (std::size_t mode = 0; mode < modeCount; ++mode) {
      const double epsilon = newtonOptions.finiteDifferenceStep
        * std::max(1.0, std::abs(initialCoefficients[mode]));
      for (double sign : {1.0, -1.0}) {
        std::vector<double> coefficients = initialCoefficients;
        coefficients[mode] += sign * epsilon;
        const GridCoordinate target = SurfaceTargets(
          basis, Synthesize(basis, coefficients), referencePoint
        )[point];
        const double displacement = Distance(
          target, diagnostic.records.front().target
        );
        if (displacement > 0.0) {
          diagnostic.finiteDifferenceDisplacementMinimum = std::min(
            diagnostic.finiteDifferenceDisplacementMinimum, displacement
          );
        }
        diagnostic.finiteDifferenceDisplacementMaximum = std::max(
          diagnostic.finiteDifferenceDisplacementMaximum, displacement
        );
      }
    }
    boundaryDiagnostics.push_back(std::move(diagnostic));
  }

  struct FixedLineSearchTrial {
    double lambda;
    Range radius;
    bool valid;
    double fixedTwoNorm;
    double fixedInfinityNorm;
    double dynamicTwoNorm;
    double fixedDynamicDifference;
    std::size_t changedStencilPoints;
    std::string failure;
    std::vector<double> coefficients;
    std::vector<double> fixedResidual;
  };
  struct FixedNewtonPrototype {
    std::vector<double> centerCoefficients;
    std::vector<double> centerFixedResidual;
    std::vector<double> centerDynamicResidual;
    Stencils frozenStencils;
    std::vector<std::vector<double>> jacobian;
    std::size_t rank;
    double condition;
    std::vector<double> delta;
    std::vector<FixedLineSearchTrial> trials;
    std::optional<std::size_t> acceptedTrial;
    std::vector<double> rebuiltDynamicResidual;
  };
  const auto runFixedNewtonPrototype = [&](
    const std::vector<double>& centerCoefficients
  ) {
    FixedNewtonPrototype prototype;
    prototype.centerCoefficients = centerCoefficients;
    const SynthesizedField centerField = Synthesize(basis, centerCoefficients);
    const std::vector<GridCoordinate> centerTargets = SurfaceTargets(
      basis, centerField, referencePoint
    );
    prototype.frozenStencils = SelectStencils(
      diagnosticIndex, centerTargets, interpolationOptions.numberOfNeighbours
    );
    prototype.centerFixedResidual = FixedStencilResidual(
      basis,
      centerCoefficients,
      grid,
      referencePoint,
      prototype.frozenStencils,
      interpolationOptions
    );
    prototype.centerDynamicResidual = residual(centerCoefficients);
    prototype.jacobian.assign(
      modeCount, std::vector<double>(modeCount, 0.0)
    );
    for (std::size_t mode = 0; mode < modeCount; ++mode) {
      const double epsilon = newtonOptions.finiteDifferenceStep
        * std::max(1.0, std::abs(centerCoefficients[mode]));
      std::vector<double> plus = centerCoefficients;
      std::vector<double> minus = centerCoefficients;
      plus[mode] += epsilon;
      minus[mode] -= epsilon;
      const std::vector<double> residualPlus = FixedStencilResidual(
        basis, plus, grid, referencePoint, prototype.frozenStencils,
        interpolationOptions
      );
      const std::vector<double> residualMinus = FixedStencilResidual(
        basis, minus, grid, referencePoint, prototype.frozenStencils,
        interpolationOptions
      );
      for (std::size_t equation = 0; equation < modeCount; ++equation) {
        prototype.jacobian[equation][mode] =
          (residualPlus[equation] - residualMinus[equation]) / (2.0 * epsilon);
      }
    }
    prototype.rank = EstimateRank(prototype.jacobian);
    prototype.condition = EstimateConditionInfinity(prototype.jacobian);
    std::vector<double> negative(prototype.centerFixedResidual.size());
    std::transform(
      prototype.centerFixedResidual.begin(),
      prototype.centerFixedResidual.end(),
      negative.begin(),
      [](double value) { return -value; }
    );
    prototype.delta = SolveLinearSystem(prototype.jacobian, negative);
    double trialLambda = 1.0;
    for (int step = 0; step < newtonOptions.maxLineSearchSteps; ++step) {
      std::vector<double> coefficients(modeCount);
      for (std::size_t mode = 0; mode < modeCount; ++mode) {
        coefficients[mode] = centerCoefficients[mode]
          + trialLambda * prototype.delta[mode];
      }
      const SynthesizedField field = Synthesize(basis, coefficients);
      const std::vector<GridCoordinate> targets = SurfaceTargets(
        basis, field, referencePoint
      );
      FixedLineSearchTrial trial{
        trialLambda,
        Summarize(field.h),
        false,
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN(),
        0,
        "",
        coefficients,
        {}
      };
      if (std::any_of(field.h.begin(), field.h.end(), [](double h) {
            return !std::isfinite(h) || h <= 0.0;
          })) {
        trial.failure = "h <= 0 or non-finite";
      } else if (!targetsInsideInterpolationDomain(targets)) {
        trial.failure = "target outside interpolation domain tolerance";
      } else {
        try {
          trial.fixedResidual = FixedStencilResidual(
            basis,
            coefficients,
            grid,
            referencePoint,
            prototype.frozenStencils,
            interpolationOptions
          );
          trial.valid = AllFinite(trial.fixedResidual);
          if (trial.valid) {
            trial.fixedTwoNorm = TwoNorm(trial.fixedResidual);
            trial.fixedInfinityNorm = InfinityNorm(trial.fixedResidual);
            const std::vector<double> dynamic = residual(coefficients);
            if (AllFinite(dynamic)) {
              trial.dynamicTwoNorm = TwoNorm(dynamic);
              trial.fixedDynamicDifference = TwoNorm(
                Difference(trial.fixedResidual, dynamic)
              );
            }
            const Stencils dynamicStencils = SelectStencils(
              diagnosticIndex,
              targets,
              interpolationOptions.numberOfNeighbours
            );
            for (std::size_t point = 0; point < basis.nPoints; ++point) {
              bool sameIdentities = true;
              for (std::size_t cell : prototype.frozenStencils[point]) {
                sameIdentities = sameIdentities
                  && std::find(
                    dynamicStencils[point].begin(),
                    dynamicStencils[point].end(),
                    cell
                  ) != dynamicStencils[point].end();
              }
              if (!sameIdentities) {
                ++trial.changedStencilPoints;
              }
            }
            if (!prototype.acceptedTrial
                && trial.fixedTwoNorm
                  < TwoNorm(prototype.centerFixedResidual)) {
              prototype.acceptedTrial = prototype.trials.size();
            }
          } else {
            trial.failure = "non-finite fixed residual";
          }
        } catch (const std::exception& error) {
          trial.failure = error.what();
        }
      }
      prototype.trials.push_back(std::move(trial));
      trialLambda *= 0.5;
    }
    if (prototype.acceptedTrial) {
      prototype.rebuiltDynamicResidual = residual(
        prototype.trials[*prototype.acceptedTrial].coefficients
      );
    }
    return prototype;
  };

  const FixedNewtonPrototype fixedIteration0 = runFixedNewtonPrototype(
    initialCoefficients
  );
  std::optional<FixedNewtonPrototype> fixedIteration1;
  if (fixedIteration0.acceptedTrial
      && AllFinite(fixedIteration0.rebuiltDynamicResidual)
      && TwoNorm(fixedIteration0.rebuiltDynamicResidual)
        < TwoNorm(fixedIteration0.centerDynamicResidual)) {
    fixedIteration1 = runFixedNewtonPrototype(
      fixedIteration0.trials[*fixedIteration0.acceptedTrial].coefficients
    );
  }

  struct TrialDiagnostic {
    double lambda;
    Range radius;
    bool valid;
    double residualTwoNorm;
    double residualInfinityNorm;
    std::string failure;
  };
  std::vector<TrialDiagnostic> trialDiagnostics;
  double lambda = 1.0;
  for (int step = 0; step < newtonOptions.maxLineSearchSteps; ++step) {
    std::vector<double> trialCoefficients(modeCount);
    for (std::size_t mode = 0; mode < modeCount; ++mode) {
      trialCoefficients[mode] = initialCoefficients[mode]
        + lambda * newtonDelta[mode];
    }
    const SynthesizedField trialField = Synthesize(basis, trialCoefficients);
    TrialDiagnostic diagnostic{
      lambda,
      Summarize(trialField.h),
      false,
      std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::quiet_NaN(),
      ""
    };
    const auto invalidRadius = std::find_if(
      trialField.h.begin(), trialField.h.end(), [](double radius) {
        return !std::isfinite(radius) || radius <= 0.0;
      }
    );
    if (invalidRadius != trialField.h.end()) {
      diagnostic.failure = "h <= 0 or non-finite";
    } else {
      try {
        (void)BuildInterpolatedSurface(
          basis, trialCoefficients, interpolator, referencePoint
        );
        const std::vector<double> trialResidual = residual(trialCoefficients);
        if (AllFinite(trialResidual)) {
          diagnostic.valid = true;
          diagnostic.residualTwoNorm = TwoNorm(trialResidual);
          diagnostic.residualInfinityNorm = InfinityNorm(trialResidual);
        } else {
          diagnostic.failure = "non-finite residual after valid surface construction";
        }
      } catch (const std::exception& error) {
        diagnostic.failure = error.what();
      }
    }
    trialDiagnostics.push_back(diagnostic);
    lambda *= 0.5;
  }

  const auto solveStart = Clock::now();
  std::optional<SpectralNewtonResult> result;
  std::string solveException;
  try {
    result = SolvePDEforAH(
      basis,
      interpolator,
      referencePoint,
      initialRadius,
      newtonOptions
    );
  } catch (const std::exception& error) {
    solveException = error.what();
  }
  const auto solveFinish = Clock::now();
  const auto radiusSolveStart = Clock::now();
  const SpectralNewtonResult radiusResult = SolvePDEforAH(
    basis,
    radiusInterpolator,
    referencePoint,
    initialRadius,
    newtonOptions
  );
  const auto radiusSolveFinish = Clock::now();
  assert(radiusResult.coefficients.size() == basis.nModes);
  assert(radiusResult.converged);
  assert(std::is_sorted(
    radiusResult.residualHistory.rbegin(), radiusResult.residualHistory.rend()
  ));
  const SynthesizedField radiusFinalField = Synthesize(
    basis, radiusResult.coefficients
  );
  const Range radiusFinalRange = Summarize(radiusFinalField.h);
  const RadiusResidualEvaluation radiusFinal = radiusResidualAt(
    radiusResult.coefficients
  );
  std::vector<double> radiusFinalTheta;
  radiusFinalTheta.reserve(basis.nPoints);
  for (const auto& point : radiusFinal.surface.points) {
    radiusFinalTheta.push_back(computeExpansion(point));
  }
  std::vector<double> radiusFinalSupportCounts(
    radiusFinal.supportCounts.begin(), radiusFinal.supportCounts.end()
  );
  const Range radiusFinalSupportRange = Summarize(radiusFinalSupportCounts);
  const std::vector<double>& finalCoefficients = result
    ? result->coefficients : initialCoefficients;
  assert(finalCoefficients.size() == basis.nModes);

  const SynthesizedField finalField = Synthesize(basis, finalCoefficients);
  const Range finalRadiusRange = Summarize(finalField.h);
  for (double radius : finalField.h) {
    RequireFinite(radius);
    assert(radius > 0.0);
  }
  RequireInsideDomain(
    basis, finalField, referencePoint, domainMinimum, domainMaximum
  );
  const InterpolatedSurface finalSurface = BuildInterpolatedSurface(
    basis, finalCoefficients, interpolator, referencePoint
  );
  assert(finalSurface.numberOfPoints() == basis.nPoints);
  std::vector<double> finalDeterminants;
  std::vector<double> finalExpansion;
  ValidateSurfaceGeometry(finalSurface, finalDeterminants, finalExpansion);

  const std::vector<double> finalResidual = residual(finalCoefficients);
  assert(finalResidual.size() == basis.nModes);
  for (double value : finalResidual) RequireFinite(value);
  const double finalResidualNorm = InfinityNorm(finalResidual);
  const double effectiveTolerance = newtonOptions.tolerance
    / (newtonOptions.massScale > 0.0 ? newtonOptions.massScale : 1.0);
  if (result && result->converged) {
    assert(finalResidualNorm < effectiveTolerance);
  }
  if (result) {
    assert(!result->residualHistory.empty());
    const double historyScale = std::max(1.0, finalResidualNorm);
    assert(std::abs(result->residualHistory.back() - finalResidualNorm)
           <= 1.0e-12 * historyScale);
  }

  const Range determinantRange = Summarize(finalDeterminants);
  const Range thetaRange = Summarize(finalExpansion);
  const double nonsphericity = finalRadiusRange.mean > 0.0
    ? (finalRadiusRange.maximum - finalRadiusRange.minimum) / finalRadiusRange.mean
    : std::numeric_limits<double>::quiet_NaN();

  std::cout << std::setprecision(17)
            << "dataset = " << path << '\n'
            << "cells = " << grid.numberOfPoints() << ", compact fields = "
            << GridData::numberOfVariables << '\n'
            << "cell-centre domain = [" << domainMinimum.x << ", " << domainMaximum.x
            << "] x [" << domainMinimum.y << ", " << domainMaximum.y
            << "] x [" << domainMinimum.z << ", " << domainMaximum.z << "]\n"
            << "estimated cell-centre spacing near inferred r=0.5 = "
            << localSpacingMinimum << " .. " << localSpacingMaximum << '\n'
            << "lmax/nModes/nPoints = " << lmax << " / " << basis.nModes
            << " / " << basis.nPoints << '\n'
            << "reference point = (0, 0, 0), exploratory initial radius = "
            << initialRadius << '\n'
            << "interpolation = weighted polynomial, degree 1, neighbours 64, "
               "weightSupport 1, rankToleranceMultiplier 64, allowOrderReduction false\n"
            << "Newton tolerance/massScale/effective tolerance = "
            << newtonOptions.tolerance << " / " << newtonOptions.massScale << " / "
            << effectiveTolerance << '\n'
            << "Newton finiteDifferenceStep/maxIterations/maxLineSearchSteps/maxRadiusFactor = "
            << newtonOptions.finiteDifferenceStep << " / " << newtonOptions.maxIterations
            << " / " << newtonOptions.maxLineSearchSteps << " / "
            << newtonOptions.maxRadiusFactor << '\n';
  PrintRange("initial h", initialRadiusRange);
  PrintVector("a0", initialCoefficients);
  PrintVector("R(a0)", initialResidual);
  std::cout << "initial ||R||_2 = " << initialResidualTwoNorm << '\n'
            << "initial ||R||_inf = " << initialResidualNorm << '\n'
            << "finite-difference residuals finite = " << std::boolalpha
            << finiteDifferenceResidualsFinite << '\n'
            << "Jacobian dimensions = " << jacobian.size() << " x "
            << (jacobian.empty() ? 0 : jacobian.front().size()) << '\n'
            << "Jacobian entries finite = " << jacobianFinite << '\n'
            << "Jacobian rank estimate = " << jacobianRank << " / "
            << modeCount << '\n'
            << "Jacobian condition estimate (infinity norm) = "
            << jacobianCondition << '\n'
            << "Jacobian determinant = not reported (conditioning diagnostic is more meaningful)\n";
  PrintVector("delta_a", newtonDelta);
  PrintVector("J_delta_a", jacobianTimesDelta);
  std::cout << "||J_delta_a + R(a0)||_2 = "
            << TwoNorm(Difference(jacobianTimesDelta, negativeResidual)) << '\n';
  std::cout << "||delta_a||_2 = " << deltaTwoNorm << '\n'
            << "||delta_a||_inf = " << deltaInfinityNorm << '\n'
            << "||delta_a||_2 / ||a0||_2 = "
            << deltaTwoNorm / TwoNorm(initialCoefficients) << '\n'
            << "||delta_a||_inf / ||a0||_inf = "
            << deltaInfinityNorm / InfinityNorm(initialCoefficients) << '\n'
            << "line-search diagnostics:\n"
            << "lambda | min(h) | mean(h) | max(h) | valid | ||R||2 | ||R||inf | failure\n";
  for (const TrialDiagnostic& diagnostic : trialDiagnostics) {
    std::cout << diagnostic.lambda << " | " << diagnostic.radius.minimum
              << " | " << diagnostic.radius.mean << " | "
              << diagnostic.radius.maximum << " | " << diagnostic.valid
              << " | " << diagnostic.residualTwoNorm << " | "
              << diagnostic.residualInfinityNorm << " | "
              << diagnostic.failure << '\n';
  }
  std::cout << "directional residual diagnostics:\n"
            << "t | min(h) | mean(h) | max(h) | ||R||2 | ||R||inf | valid/failure | R components\n";
  for (const DirectionalSample& sample : directionalSamples) {
    std::cout << sample.t << " | " << sample.radius.minimum << " | "
              << sample.radius.mean << " | " << sample.radius.maximum << " | "
              << sample.residualTwoNorm << " | " << sample.residualInfinityNorm
              << " | " << (sample.valid ? "valid" : sample.failure) << " | [";
    for (std::size_t component = 0; component < sample.residual.size(); ++component) {
      if (component != 0) std::cout << ", ";
      std::cout << sample.residual[component];
    }
    std::cout << "]\n";
  }
  const double derivativeReferenceNorm = TwoNorm(jacobianTimesDelta);
  std::cout << "directional derivative diagnostics:\n"
            << "t | plus error2 | plus relative | minus error2 | minus relative | central error2 | central relative | positive linear-model error2\n";
  for (std::size_t scale = 0; scale < positiveScales.size(); ++scale) {
    const DirectionalSample& plus = directionalSamples[2 * scale];
    const DirectionalSample& minus = directionalSamples[2 * scale + 1];
    if (!plus.valid || !minus.valid) continue;
    const double t = positiveScales[scale];
    const std::vector<double> dPlus = DifferenceScaled(
      plus.residual, initialResidual, 1.0 / t
    );
    const std::vector<double> dMinus = DifferenceScaled(
      initialResidual, minus.residual, 1.0 / t
    );
    const std::vector<double> dCentral = DifferenceScaled(
      plus.residual, minus.residual, 1.0 / (2.0 * t)
    );
    const double plusError = TwoNorm(Difference(dPlus, jacobianTimesDelta));
    const double minusError = TwoNorm(Difference(dMinus, jacobianTimesDelta));
    const double centralError = TwoNorm(Difference(dCentral, jacobianTimesDelta));
    std::vector<double> linearResidual(initialResidual.size());
    for (std::size_t component = 0;
         component < initialResidual.size(); ++component) {
      linearResidual[component] = (1.0 - t) * initialResidual[component];
    }
    std::cout << t << " | " << plusError << " | "
              << plusError / derivativeReferenceNorm << " | " << minusError
              << " | " << minusError / derivativeReferenceNorm << " | "
              << centralError << " | "
              << centralError / derivativeReferenceNorm << " | "
              << TwoNorm(Difference(plus.residual, linearResidual)) << '\n';
  }
  std::cout << "stencil-change diagnostics (changed count is identities replaced from the 64-cell a0 stencil):\n"
            << "t | affected points | mean changed neighbours | max changed neighbours | nearest changed points | max target displacement\n";
  for (const DirectionalSample& sample : directionalSamples) {
    std::cout << sample.t << " | " << sample.affectedPoints << " | "
              << sample.meanChangedNeighbours << " | "
              << sample.maximumChangedNeighbours << " | "
              << sample.nearestChangedPoints << " | "
              << sample.maximumDisplacement << '\n';
  }
  std::cout << "smallest |t| with a stencil change by collocation point = [";
  for (std::size_t point = 0; point < smallestStencilChange.size(); ++point) {
    if (point != 0) std::cout << ", ";
    if (std::isfinite(smallestStencilChange[point])) {
      std::cout << smallestStencilChange[point];
    } else {
      std::cout << "none";
    }
  }
  std::cout << "]\n"
            << "fixed-versus-dynamic stencil diagnostics:\n"
            << "t | ||R_dynamic-R_fixed||2 (+) | ||R_dynamic-R_fixed||2 (-) | fixed central derivative error2 | fixed central relative\n";
  for (std::size_t scale = 0; scale < positiveScales.size(); ++scale) {
    const DirectionalSample& plus = directionalSamples[2 * scale];
    const DirectionalSample& minus = directionalSamples[2 * scale + 1];
    if (!plus.valid || !minus.valid
        || plus.fixedResidual.empty() || minus.fixedResidual.empty()) continue;
    const double t = positiveScales[scale];
    const std::vector<double> fixedCentral = DifferenceScaled(
      plus.fixedResidual, minus.fixedResidual, 1.0 / (2.0 * t)
    );
    const double fixedCentralError = TwoNorm(
      Difference(fixedCentral, jacobianTimesDelta)
    );
    std::cout << t << " | "
              << TwoNorm(Difference(plus.residual, plus.fixedResidual)) << " | "
              << TwoNorm(Difference(minus.residual, minus.fixedResidual)) << " | "
              << fixedCentralError << " | "
              << fixedCentralError / derivativeReferenceNorm << '\n';
  }
  std::cout << "64/65 boundary diagnostics (one-based ranks):\n";
  for (const BoundaryPointDiagnostic& diagnostic : boundaryDiagnostics) {
    const bool affected = diagnostic.leavingCell
      != std::numeric_limits<std::size_t>::max();
    std::cout << "boundary point p=" << diagnostic.point
              << (diagnostic.point == 0 ? " (unaffected control)" : "") << '\n'
              << "leaving cell = ";
    if (affected) std::cout << diagnostic.leavingCell;
    else std::cout << "none";
    std::cout << ", entering cell = ";
    if (affected) std::cout << diagnostic.enteringCell;
    else std::cout << "none";
    std::cout << ", estimated crossing t = " << diagnostic.crossingT
              << ", crossing displacement = "
              << diagnostic.crossingDisplacement << '\n'
              << "finite-difference target displacement min/max = "
              << diagnostic.finiteDifferenceDisplacementMinimum << " / "
              << diagnostic.finiteDifferenceDisplacementMaximum << '\n'
              << "t | target x y z | rank64 cell d64 | rank65 cell d65 | gap | relative gap | leaving distance | entering distance\n";
    for (const BoundaryRecord& record : diagnostic.records) {
      const std::size_t cell64 = record.neighbours[63];
      const std::size_t cell65 = record.neighbours[64];
      const double d64 = Distance(record.target, grid.coordinates[cell64]);
      const double d65 = Distance(record.target, grid.coordinates[cell65]);
      const double leavingDistance = affected
        ? Distance(record.target, grid.coordinates[diagnostic.leavingCell])
        : std::numeric_limits<double>::quiet_NaN();
      const double enteringDistance = affected
        ? Distance(record.target, grid.coordinates[diagnostic.enteringCell])
        : std::numeric_limits<double>::quiet_NaN();
      std::cout << record.t << " | " << record.target.x << " "
                << record.target.y << " " << record.target.z << " | "
                << cell64 << " " << d64 << " | " << cell65 << " "
                << d65 << " | " << d65 - d64 << " | "
                << (d65 - d64) / d64 << " | " << leavingDistance
                << " | " << enteringDistance << '\n';
      std::cout << "  ranks 60..68: ";
      for (std::size_t rank = 59; rank < 68; ++rank) {
        const std::size_t cell = record.neighbours[rank];
        const double squared = SquaredDistance(
          record.target, grid.coordinates[cell]
        );
        if (rank != 59) std::cout << "; ";
        std::cout << (rank + 1) << ":cell=" << cell << ",d2=" << squared
                  << ",d=" << std::sqrt(squared);
      }
      std::cout << '\n'
                << "  d63/d64/d65/d66 = ";
      for (std::size_t rank = 62; rank <= 65; ++rank) {
        if (rank != 62) std::cout << " / ";
        std::cout << Distance(
          record.target, grid.coordinates[record.neighbours[rank]]
        );
      }
      std::cout << '\n';
    }
    if (affected) {
      std::cout << "crossing displacement / local spacing = "
                << diagnostic.crossingDisplacement / localSpacingMinimum
                << '\n';
    }
  }
  const auto printFixedPrototype = [&](
    const char* name,
    const FixedNewtonPrototype& prototype
  ) {
    std::cout << name << " fixed-stencil Newton prototype:\n"
              << "rank = " << prototype.rank << " / " << modeCount
              << ", condition(infinity norm) = " << prototype.condition << '\n';
    PrintVector("delta_fixed", prototype.delta);
    std::cout << "||delta_fixed||2 = " << TwoNorm(prototype.delta)
              << ", ||delta_fixed||inf = " << InfinityNorm(prototype.delta)
              << '\n'
              << "||R_dynamic(center)||2 = "
              << TwoNorm(prototype.centerDynamicResidual)
              << ", ||R_fixed(center)||2 = "
              << TwoNorm(prototype.centerFixedResidual)
              << ", ||R_dynamic-R_fixed||2 at center = "
              << TwoNorm(Difference(
                   prototype.centerDynamicResidual,
                   prototype.centerFixedResidual
                 )) << '\n'
              << "lambda | min(h) | mean(h) | max(h) | valid | ||R_fixed||2 | ||R_fixed||inf | ||R_dynamic||2 | ||fixed-dynamic||2 | changed stencil points | accepted | failure\n";
    for (std::size_t index = 0; index < prototype.trials.size(); ++index) {
      const FixedLineSearchTrial& trial = prototype.trials[index];
      std::cout << trial.lambda << " | " << trial.radius.minimum << " | "
                << trial.radius.mean << " | " << trial.radius.maximum << " | "
                << trial.valid << " | " << trial.fixedTwoNorm << " | "
                << trial.fixedInfinityNorm << " | " << trial.dynamicTwoNorm
                << " | " << trial.fixedDynamicDifference << " | "
                << trial.changedStencilPoints << " | "
                << (prototype.acceptedTrial
                    && *prototype.acceptedTrial == index)
                << " | " << trial.failure << '\n';
    }
    if (prototype.acceptedTrial) {
      const FixedLineSearchTrial& accepted =
        prototype.trials[*prototype.acceptedTrial];
      std::cout << "first accepted lambda = " << accepted.lambda << '\n'
                << "||R_fixed(accepted)||2 = "
                << accepted.fixedTwoNorm << '\n'
                << "||R_dynamic(rebuilt accepted)||2 = "
                << TwoNorm(prototype.rebuiltDynamicResidual) << '\n'
                << "rebuilt production residual decreases = "
                << (TwoNorm(prototype.rebuiltDynamicResidual)
                    < TwoNorm(prototype.centerDynamicResidual)) << '\n';
    } else {
      std::cout << "first accepted lambda = none\n";
    }
  };
  printFixedPrototype("iteration 0", fixedIteration0);
  if (fixedIteration1) {
    printFixedPrototype("iteration 1", *fixedIteration1);
  } else {
    std::cout << "optional iteration 1 not run because rebuilt production residual did not decrease\n";
  }
  std::vector<double> radiusSupportCountsAsDouble(
    radiusInitial.supportCounts.begin(), radiusInitial.supportCounts.end()
  );
  const Range radiusSupportRange = Summarize(radiusSupportCountsAsDouble);
  const Range radiusFitConditionRange = Summarize(radiusInitial.fitConditions);
  double maximumMetricDifference = 0.0;
  double maximumMetricDerivativeDifference = 0.0;
  double maximumCurvatureDifference = 0.0;
  double maximumTraceKDifference = 0.0;
  double maximumThetaDifference = 0.0;
  for (std::size_t point = 0; point < basis.nPoints; ++point) {
    const auto& dynamicPoint = initialSurface.points[point];
    const auto& radiusPoint = radiusInitial.surface.points[point];
    maximumTraceKDifference = std::max(
      maximumTraceKDifference,
      std::abs(dynamicPoint.geometry.traceK - radiusPoint.geometry.traceK)
    );
    maximumThetaDifference = std::max(
      maximumThetaDifference,
      std::abs(computeExpansion(dynamicPoint) - computeExpansion(radiusPoint))
    );
    for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 3; ++j) {
        maximumMetricDifference = std::max(
          maximumMetricDifference,
          std::abs(dynamicPoint.geometry.spatialMetric[i][j]
                   - radiusPoint.geometry.spatialMetric[i][j])
        );
        maximumCurvatureDifference = std::max(
          maximumCurvatureDifference,
          std::abs(dynamicPoint.geometry.extrinsicCurvature[i][j]
                   - radiusPoint.geometry.extrinsicCurvature[i][j])
        );
        for (std::size_t k = 0; k < 3; ++k) {
          maximumMetricDerivativeDifference = std::max(
            maximumMetricDerivativeDifference,
            std::abs(dynamicPoint.geometry.spatialMetricDerivative[k][i][j]
                     - radiusPoint.geometry.spatialMetricDerivative[k][i][j])
          );
        }
      }
    }
  }
  std::cout << "radius-supported Wendland-C2 MLS prototype:\n"
            << "kernel w(q)=(1-q)^4(4q+1) for 0<=q<1, zero otherwise\n"
            << "support radius/query = " << radiusSupport
            << " / true k-d-tree radius query\n"
            << "initial support min/mean/max = " << radiusSupportRange.minimum
            << " / " << radiusSupportRange.mean << " / "
            << radiusSupportRange.maximum << '\n'
            << "initial support p2/p3/p12/p13 = "
            << radiusInitial.supportCounts[2] << " / "
            << radiusInitial.supportCounts[3] << " / "
            << radiusInitial.supportCounts[12] << " / "
            << radiusInitial.supportCounts[13] << '\n'
            << "fit rank = 4, fit condition min/mean/max = "
            << radiusFitConditionRange.minimum << " / "
            << radiusFitConditionRange.mean << " / "
            << radiusFitConditionRange.maximum << '\n'
            << "basic constant/linear/gradient reproduction = passed\n"
            << "dynamic-vs-radius initial max differences gamma/dgamma/Kij/K/Theta = "
            << maximumMetricDifference << " / "
            << maximumMetricDerivativeDifference << " / "
            << maximumCurvatureDifference << " / "
            << maximumTraceKDifference << " / "
            << maximumThetaDifference << '\n'
            << "||R_dynamic(a0)-R_radius(a0)||2 = "
            << TwoNorm(Difference(initialResidual, radiusInitial.residual)) << '\n'
            << "||R_radius(a0)||2/inf = " << TwoNorm(radiusInitial.residual)
            << " / " << InfinityNorm(radiusInitial.residual) << '\n'
            << "target-motion support/residual study along production delta:\n"
            << "t | support counts | changed support points | ||R_radius||2 | ||R_radius||inf\n";
  for (const RadiusMotionSample& sample : radiusMotionSamples) {
    std::cout << sample.t << " | [";
    for (std::size_t point = 0; point < sample.supportCounts.size(); ++point) {
      if (point != 0) std::cout << ",";
      std::cout << sample.supportCounts[point];
    }
    std::cout << "] | " << sample.changedSupportPoints << " | "
              << sample.twoNorm << " | " << sample.infinityNorm << '\n';
  }
  std::cout << "radius directional derivative study along delta_radius:\n"
            << "t | plus error2 | minus error2 | central error2 | central relative | changed supports +/-\n";
  for (const RadiusDerivativeSample& sample : radiusDerivativeSamples) {
    std::cout << sample.t << " | " << sample.plusError << " | "
              << sample.minusError << " | " << sample.centralError << " | "
              << sample.centralRelative << " | " << sample.plusChangedPoints
              << " / " << sample.minusChangedPoints << '\n';
  }
  std::cout << "radius Jacobian rank/condition = "
            << EstimateRank(radiusJacobian) << " / " << modeCount << " / "
            << EstimateConditionInfinity(radiusJacobian) << '\n';
  PrintVector("delta_radius", radiusDelta);
  std::cout << "||delta_radius||2/inf = " << TwoNorm(radiusDelta) << " / "
            << InfinityNorm(radiusDelta) << '\n'
            << "radius line search:\n"
            << "lambda | min(h) | mean(h) | max(h) | valid | ||R||2 | ||R||inf | accepted\n";
  for (const RadiusLineTrial& trial : radiusLineTrials) {
    std::cout << trial.lambda << " | " << trial.radius.minimum << " | "
              << trial.radius.mean << " | " << trial.radius.maximum << " | "
              << trial.valid << " | " << trial.twoNorm << " | "
              << trial.infinityNorm << " | " << trial.accepted << '\n';
  }
  std::cout << "optional radius Newton iterations:\n"
            << "iteration | before ||R||2 | accepted lambda | after ||R||2 | h min/mean/max | rank | condition\n";
  for (const RadiusShortIteration& iteration : radiusShortIterations) {
    std::cout << iteration.iteration << " | " << iteration.beforeNorm << " | "
              << iteration.acceptedLambda << " | " << iteration.afterNorm
              << " | " << iteration.radius.minimum << " / "
              << iteration.radius.mean << " / " << iteration.radius.maximum
              << " | " << iteration.rank << " | " << iteration.condition
              << '\n';
  }
  std::cout << "production radius-supported SolvePDEforAH result:\n"
            << "converged/iterations/message = " << radiusResult.converged
            << " / " << radiusResult.iterations << " / " << radiusResult.message
            << '\n' << "residualHistory = [";
  for (std::size_t i = 0; i < radiusResult.residualHistory.size(); ++i) {
    if (i != 0) std::cout << ", ";
    std::cout << radiusResult.residualHistory[i];
  }
  std::cout << "]\n"
            << "final h min/mean/max = " << radiusFinalRange.minimum << " / "
            << radiusFinalRange.mean << " / " << radiusFinalRange.maximum << '\n'
            << "final ||R||inf / max(abs(Theta)) = "
            << InfinityNorm(radiusFinal.residual) << " / "
            << InfinityNorm(radiusFinalTheta) << '\n'
            << "final support min/mean/max = "
            << radiusFinalSupportRange.minimum << " / "
            << radiusFinalSupportRange.mean << " / "
            << radiusFinalSupportRange.maximum << '\n'
            << "radius Newton seconds = "
            << Seconds(radiusSolveStart, radiusSolveFinish) << '\n';
  std::cout << "Newton returned result = " << std::boolalpha << result.has_value() << '\n';
  if (result) {
    std::cout << "converged = " << result->converged
              << ", iterations = " << result->iterations
              << ", message = " << result->message << '\n'
              << "residualHistory = [";
    for (std::size_t i = 0; i < result->residualHistory.size(); ++i) {
      if (i != 0) std::cout << ", ";
      std::cout << result->residualHistory[i];
    }
    std::cout << "]\n";
  } else {
    std::cout << "first Newton failure stage = trial Cartesian interpolation\n"
              << "exception = " << solveException << '\n'
              << "SpectralNewtonResult/residualHistory unavailable because the residual exception propagated\n";
  }
  PrintRange("final h", finalRadiusRange);
  std::cout << "relative nonsphericity = " << nonsphericity << '\n'
            << "det(gamma) min/max = " << determinantRange.minimum << " / "
            << determinantRange.maximum << '\n'
            << "Theta min/max/max(abs) = " << thetaRange.minimum << " / "
            << thetaRange.maximum << " / " << InfinityNorm(finalExpansion) << '\n'
            << "final ||R||_inf = " << finalResidualNorm << '\n'
            << "timing seconds: read/index/initial-residual/Newton = "
            << Seconds(readStart, readFinish) << " / "
            << Seconds(indexStart, indexFinish) << " / "
            << Seconds(initialResidualStart, initialResidualFinish) << " / "
            << Seconds(solveStart, solveFinish) << '\n'
            << "persistent interpolator reused for initial residual, Newton, and final surface\n";
}

}  // namespace
}  // namespace exahype2::ccz4::ahf::tests::integration

int main() {
  exahype2::ccz4::ahf::tests::integration::TestSampleDataSingleBlackHolePipeline();
  return 0;
}
