#include "interpolation.h"

#include <algorithm>
#include <array>
#include <limits>
#include <numeric>
#include <sstream>

namespace exahype2::ccz4::ahf {

double LinearReconstructionWeight(
  double distance,
  double lengthScale,
  const LinearReconstructionOptions& options
) {
  if (!std::isfinite(distance) || distance < 0.0
      || !std::isfinite(lengthScale) || lengthScale <= 0.0) {
    throw std::invalid_argument(
      "LinearReconstructionWeight: distance and length scale must be finite and valid");
  }
  if (options.method == LinearReconstructionMethod::Unweighted) {
    return 1.0;
  }
  if (!std::isfinite(options.weightSupport) || options.weightSupport <= 0.0) {
    throw std::invalid_argument(
      "LinearReconstructionWeight: weight support must be finite and strictly positive");
  }
  const double scaledDistance = distance / (options.weightSupport * lengthScale);
  const double weight = std::exp(-scaledDistance * scaledDistance);
  if (!std::isfinite(weight) || weight <= 0.0) {
    throw std::invalid_argument(
      "LinearReconstructionWeight: computed weight must be finite and strictly positive");
  }
  return weight;
}

GridCoordinate SurfacePointToCartesian(
  const AHGuessPoint& surfacePoint,
  const ReferencePoint& referencePoint
) {
  if (!std::isfinite(surfacePoint.r)
      || !std::isfinite(surfacePoint.rho)
      || !std::isfinite(surfacePoint.sigma)
      || !std::isfinite(referencePoint.x)
      || !std::isfinite(referencePoint.y)
      || !std::isfinite(referencePoint.z)) {
    throw std::invalid_argument(
      "SurfacePointToCartesian: radius, angles, and reference point must be finite");
  }
  if (surfacePoint.r <= 0.0) {
    throw std::invalid_argument(
      "SurfacePointToCartesian: apparent-horizon radius must be strictly positive");
  }

  const double cosSigma = std::cos(surfacePoint.sigma);
  return GridCoordinate{
    referencePoint.x + surfacePoint.r * cosSigma * std::cos(surfacePoint.rho),
    referencePoint.y + surfacePoint.r * cosSigma * std::sin(surfacePoint.rho),
    referencePoint.z + surfacePoint.r * std::sin(surfacePoint.sigma)
  };
}

std::vector<std::size_t> FindNearestGridCells(
  const GridData& gridData,
  const GridCoordinate& target,
  std::size_t neighbourCount
) {
  if (gridData.coordinates.empty()) {
    throw std::invalid_argument(
      "FindNearestGridCells: cannot search an empty Cartesian grid");
  }
  if (neighbourCount == 0) {
    throw std::invalid_argument(
      "FindNearestGridCells: neighbour count must be greater than zero");
  }
  if (neighbourCount > gridData.coordinates.size()) {
    throw std::invalid_argument(
      "FindNearestGridCells: requested " + std::to_string(neighbourCount)
      + " neighbours from only " + std::to_string(gridData.coordinates.size())
      + " cells");
  }
  if (!std::isfinite(target.x)
      || !std::isfinite(target.y)
      || !std::isfinite(target.z)) {
    throw std::invalid_argument(
      "FindNearestGridCells: target coordinate must be finite");
  }

  std::vector<std::size_t> indices(gridData.coordinates.size());
  std::iota(indices.begin(), indices.end(), std::size_t{0});

  const auto squaredDistance = [&gridData, &target](std::size_t cellIndex) {
    const GridCoordinate& coordinate = gridData.coordinates[cellIndex];
    const double dx = coordinate.x - target.x;
    const double dy = coordinate.y - target.y;
    const double dz = coordinate.z - target.z;
    return dx * dx + dy * dy + dz * dz;
  };
  const auto nearer = [&squaredDistance](std::size_t first, std::size_t second) {
    const double firstDistance = squaredDistance(first);
    const double secondDistance = squaredDistance(second);
    return firstDistance < secondDistance
        || (firstDistance == secondDistance && first < second);
  };

  std::partial_sort(
    indices.begin(),
    indices.begin() + neighbourCount,
    indices.end(),
    nearer
  );
  indices.resize(neighbourCount);
  return indices;
}

LinearPolynomialFit FitLinearPolynomial(
  std::span<const GridCoordinate> sampleCoordinates,
  std::span<const double> sampleValues,
  const GridCoordinate& target,
  const LinearReconstructionOptions& options
) {
  constexpr std::size_t BasisSize = 4;
  if (sampleCoordinates.size() != sampleValues.size()) {
    throw std::invalid_argument(
      "FitLinearPolynomial: coordinate and value counts must match");
  }
  const std::size_t requestedCount = options.numberOfNeighbours == 0
    ? sampleCoordinates.size() : options.numberOfNeighbours;
  if (requestedCount > sampleCoordinates.size()) {
    throw std::invalid_argument(
      "FitLinearPolynomial: requested neighbour count exceeds supplied samples");
  }
  sampleCoordinates = sampleCoordinates.first(requestedCount);
  sampleValues = sampleValues.first(requestedCount);
  if (sampleCoordinates.size() < BasisSize) {
    throw std::invalid_argument(
      "FitLinearPolynomial: degree-1 reconstruction requires at least four samples");
  }
  if (!std::isfinite(target.x) || !std::isfinite(target.y)
      || !std::isfinite(target.z)) {
    throw std::invalid_argument("FitLinearPolynomial: target must be finite");
  }
  if (!std::isfinite(options.rankToleranceMultiplier)
      || options.rankToleranceMultiplier <= 0.0) {
    throw std::invalid_argument(
      "FitLinearPolynomial: rank tolerance multiplier must be finite and positive");
  }

  double lengthScale = 0.0;
  for (std::size_t row = 0; row < sampleCoordinates.size(); ++row) {
    const GridCoordinate& sample = sampleCoordinates[row];
    if (!std::isfinite(sample.x) || !std::isfinite(sample.y)
        || !std::isfinite(sample.z) || !std::isfinite(sampleValues[row])) {
      throw std::invalid_argument(
        "FitLinearPolynomial: sample coordinates and values must be finite");
    }
    lengthScale = std::max(
      lengthScale,
      std::hypot(sample.x - target.x, sample.y - target.y, sample.z - target.z)
    );
  }
  if (!(lengthScale > 0.0) || !std::isfinite(lengthScale)) {
    throw std::invalid_argument(
      "FitLinearPolynomial: samples are geometrically rank deficient");
  }

  std::vector<double> squareRootWeights(sampleCoordinates.size());
  for (std::size_t row = 0; row < sampleCoordinates.size(); ++row) {
    const GridCoordinate& sample = sampleCoordinates[row];
    const double distance = std::hypot(
      sample.x - target.x, sample.y - target.y, sample.z - target.z);
    squareRootWeights[row] = std::sqrt(
      LinearReconstructionWeight(distance, lengthScale, options));
  }

  std::array<std::vector<double>, BasisSize> columns;
  for (auto& column : columns) {
    column.resize(sampleCoordinates.size());
  }
  for (std::size_t row = 0; row < sampleCoordinates.size(); ++row) {
    columns[0][row] = 1.0;
    columns[1][row] = (sampleCoordinates[row].x - target.x) / lengthScale;
    columns[2][row] = (sampleCoordinates[row].y - target.y) / lengthScale;
    columns[3][row] = (sampleCoordinates[row].z - target.z) / lengthScale;
    for (std::size_t column = 0; column < BasisSize; ++column) {
      columns[column][row] *= squareRootWeights[row];
    }
  }

  const auto dot = [](const std::vector<double>& first,
                      const std::vector<double>& second) {
    return std::inner_product(first.begin(), first.end(), second.begin(), 0.0);
  };
  std::array<double, BasisSize> squaredNorm{};
  for (std::size_t column = 0; column < BasisSize; ++column) {
    squaredNorm[column] = dot(columns[column], columns[column]);
  }
  const double largestInitialNorm =
    std::sqrt(*std::max_element(squaredNorm.begin(), squaredNorm.end()));
  const double rankTolerance =
    options.rankToleranceMultiplier * std::numeric_limits<double>::epsilon()
    * static_cast<double>(sampleCoordinates.size()) * largestInitialNorm;

  std::array<std::array<double, BasisSize>, BasisSize> upperTriangular{};
  std::array<std::size_t, BasisSize> permutation{0, 1, 2, 3};
  for (std::size_t pivot = 0; pivot < BasisSize; ++pivot) {
    const auto largest = std::max_element(
      squaredNorm.begin() + pivot,
      squaredNorm.end()
    );
    const std::size_t selected =
      static_cast<std::size_t>(largest - squaredNorm.begin());
    if (selected != pivot) {
      std::swap(columns[pivot], columns[selected]);
      std::swap(squaredNorm[pivot], squaredNorm[selected]);
      std::swap(permutation[pivot], permutation[selected]);
      for (std::size_t previous = 0; previous < pivot; ++previous) {
        std::swap(
          upperTriangular[previous][pivot],
          upperTriangular[previous][selected]
        );
      }
    }

    const double diagonal = std::sqrt(squaredNorm[pivot]);
    if (!(diagonal > rankTolerance) || !std::isfinite(diagonal)) {
      throw std::invalid_argument(
        "FitLinearPolynomial: samples are geometrically rank deficient");
    }
    upperTriangular[pivot][pivot] = diagonal;
    for (double& entry : columns[pivot]) {
      entry /= diagonal;
    }

    for (std::size_t column = pivot + 1; column < BasisSize; ++column) {
      upperTriangular[pivot][column] = dot(columns[pivot], columns[column]);
      for (std::size_t row = 0; row < sampleCoordinates.size(); ++row) {
        columns[column][row] -=
          upperTriangular[pivot][column] * columns[pivot][row];
      }
      squaredNorm[column] = dot(columns[column], columns[column]);
    }
  }

  std::array<double, BasisSize> projectedValues{};
  for (std::size_t column = 0; column < BasisSize; ++column) {
    for (std::size_t row = 0; row < sampleValues.size(); ++row) {
      projectedValues[column] +=
        columns[column][row] * sampleValues[row] * squareRootWeights[row];
    }
  }
  std::array<double, BasisSize> pivotedCoefficients{};
  for (std::size_t reverse = BasisSize; reverse-- > 0;) {
    double residual = projectedValues[reverse];
    for (std::size_t column = reverse + 1; column < BasisSize; ++column) {
      residual -= upperTriangular[reverse][column] * pivotedCoefficients[column];
    }
    pivotedCoefficients[reverse] = residual / upperTriangular[reverse][reverse];
  }
  std::array<double, BasisSize> coefficients{};
  for (std::size_t column = 0; column < BasisSize; ++column) {
    coefficients[permutation[column]] = pivotedCoefficients[column];
  }

  LinearPolynomialFit result{
    .value = coefficients[0],
    .derivative = {
      coefficients[1] / lengthScale,
      coefficients[2] / lengthScale,
      coefficients[3] / lengthScale
    },
    .lengthScale = lengthScale
  };
  if (!std::isfinite(result.value)
      || !std::isfinite(result.derivative[0])
      || !std::isfinite(result.derivative[1])
      || !std::isfinite(result.derivative[2])) {
    throw std::runtime_error("FitLinearPolynomial: QR solve produced a non-finite result");
  }
  return result;
}

CompactLinearFit FitCompactLinearPolynomial(
  const GridData& gridData,
  std::span<const std::size_t> neighbourIndices,
  const GridCoordinate& target,
  const LinearReconstructionOptions& options
) {
  constexpr std::size_t BasisSize = 4;
  constexpr std::size_t FieldCount = NumberOfCompactFields;

  // 1. Validate the stencil, target, and fixed 14-column grid layout.
  const std::size_t requestedCount = options.numberOfNeighbours == 0
    ? neighbourIndices.size() : options.numberOfNeighbours;
  if (requestedCount > neighbourIndices.size()) {
    throw std::invalid_argument(
      "FitCompactLinearPolynomial: requested neighbour count exceeds supplied neighbours");
  }
  neighbourIndices = neighbourIndices.first(requestedCount);
  if (neighbourIndices.size() < BasisSize) {
    throw std::invalid_argument(
      "FitCompactLinearPolynomial: degree-1 reconstruction requires at least four neighbours");
  }
  if (!std::isfinite(target.x) || !std::isfinite(target.y)
      || !std::isfinite(target.z)) {
    throw std::invalid_argument(
      "FitCompactLinearPolynomial: target must be finite");
  }
  if (!gridData.isConsistent()) {
    throw std::invalid_argument(
      "FitCompactLinearPolynomial: GridData must contain fourteen values per cell");
  }
  if (!std::isfinite(options.rankToleranceMultiplier)
      || options.rankToleranceMultiplier <= 0.0) {
    throw std::invalid_argument(
      "FitCompactLinearPolynomial: rank tolerance multiplier must be finite and positive");
  }

  // 2. Set ell to the farthest-neighbour distance and validate coordinates.
  double lengthScale = 0.0;
  for (std::size_t cellIndex : neighbourIndices) {
    if (cellIndex >= gridData.coordinates.size()) {
      throw std::out_of_range(
        "FitCompactLinearPolynomial: neighbour cell index is out of range");
    }
    const GridCoordinate& sample = gridData.coordinates[cellIndex];
    if (!std::isfinite(sample.x) || !std::isfinite(sample.y)
        || !std::isfinite(sample.z)) {
      throw std::invalid_argument(
        "FitCompactLinearPolynomial: neighbour coordinates must be finite");
    }
    lengthScale = std::max(
      lengthScale,
      std::hypot(sample.x - target.x, sample.y - target.y, sample.z - target.z)
    );
  }
  if (!(lengthScale > 0.0) || !std::isfinite(lengthScale)) {
    throw std::invalid_argument(
      "FitCompactLinearPolynomial: neighbours are geometrically rank deficient");
  }

  std::vector<double> squareRootWeights(neighbourIndices.size());
  for (std::size_t row = 0; row < neighbourIndices.size(); ++row) {
    const GridCoordinate& sample = gridData.coordinates[neighbourIndices[row]];
    const double distance = std::hypot(
      sample.x - target.x, sample.y - target.y, sample.z - target.z);
    squareRootWeights[row] = std::sqrt(
      LinearReconstructionWeight(distance, lengthScale, options));
  }

  // 3. Build the shared n-by-4 design matrix [1, X, Y, Z].
  std::array<std::vector<double>, BasisSize> columns;
  for (auto& column : columns) {
    column.resize(neighbourIndices.size());
  }
  for (std::size_t row = 0; row < neighbourIndices.size(); ++row) {
    const GridCoordinate& sample = gridData.coordinates[neighbourIndices[row]];
    columns[0][row] = 1.0;
    columns[1][row] = (sample.x - target.x) / lengthScale;
    columns[2][row] = (sample.y - target.y) / lengthScale;
    columns[3][row] = (sample.z - target.z) / lengthScale;
    for (std::size_t column = 0; column < BasisSize; ++column) {
      columns[column][row] *= squareRootWeights[row];
    }
  }

  // 4. Factor the design matrix once using column-pivoted modified Gram-Schmidt QR.
  const auto dot = [](const std::vector<double>& first,
                      const std::vector<double>& second) {
    return std::inner_product(first.begin(), first.end(), second.begin(), 0.0);
  };
  std::array<double, BasisSize> squaredNorm{};
  for (std::size_t column = 0; column < BasisSize; ++column) {
    squaredNorm[column] = dot(columns[column], columns[column]);
  }
  const double largestInitialNorm =
    std::sqrt(*std::max_element(squaredNorm.begin(), squaredNorm.end()));
  const double rankTolerance =
    options.rankToleranceMultiplier * std::numeric_limits<double>::epsilon()
    * static_cast<double>(neighbourIndices.size()) * largestInitialNorm;

  std::array<std::array<double, BasisSize>, BasisSize> upperTriangular{};
  std::array<std::size_t, BasisSize> permutation{0, 1, 2, 3};
  for (std::size_t pivot = 0; pivot < BasisSize; ++pivot) {
    const auto largest = std::max_element(
      squaredNorm.begin() + pivot, squaredNorm.end());
    const std::size_t selected =
      static_cast<std::size_t>(largest - squaredNorm.begin());
    if (selected != pivot) {
      std::swap(columns[pivot], columns[selected]);
      std::swap(squaredNorm[pivot], squaredNorm[selected]);
      std::swap(permutation[pivot], permutation[selected]);
      for (std::size_t previous = 0; previous < pivot; ++previous) {
        std::swap(
          upperTriangular[previous][pivot],
          upperTriangular[previous][selected]);
      }
    }

    const double diagonal = std::sqrt(squaredNorm[pivot]);
    if (!(diagonal > rankTolerance) || !std::isfinite(diagonal)) {
      throw std::invalid_argument(
        "FitCompactLinearPolynomial: neighbours are geometrically rank deficient");
    }
    upperTriangular[pivot][pivot] = diagonal;
    for (double& entry : columns[pivot]) {
      entry /= diagonal;
    }
    for (std::size_t column = pivot + 1; column < BasisSize; ++column) {
      upperTriangular[pivot][column] = dot(columns[pivot], columns[column]);
      for (std::size_t row = 0; row < neighbourIndices.size(); ++row) {
        columns[column][row] -=
          upperTriangular[pivot][column] * columns[pivot][row];
      }
      squaredNorm[column] = dot(columns[column], columns[column]);
    }
  }

  // 5. Apply Q^T to all 14 right-hand sides at the same time.
  std::array<std::array<double, FieldCount>, BasisSize> projected{};
  for (std::size_t basis = 0; basis < BasisSize; ++basis) {
    for (std::size_t row = 0; row < neighbourIndices.size(); ++row) {
      const auto compactRow = gridData.cell(neighbourIndices[row]).row();
      for (std::size_t field = 0; field < FieldCount; ++field) {
        if (!std::isfinite(compactRow[field])) {
          throw std::invalid_argument(
            "FitCompactLinearPolynomial: compact values must be finite");
        }
        projected[basis][field] +=
          columns[basis][row] * compactRow[field] * squareRootWeights[row];
      }
    }
  }

  // 6. Solve R C = Q^T B and undo the QR column permutation.
  std::array<std::array<double, FieldCount>, BasisSize> pivotedCoefficients{};
  for (std::size_t reverse = BasisSize; reverse-- > 0;) {
    for (std::size_t field = 0; field < FieldCount; ++field) {
      double residual = projected[reverse][field];
      for (std::size_t column = reverse + 1; column < BasisSize; ++column) {
        residual -= upperTriangular[reverse][column]
                  * pivotedCoefficients[column][field];
      }
      pivotedCoefficients[reverse][field] =
        residual / upperTriangular[reverse][reverse];
    }
  }
  std::array<std::array<double, FieldCount>, BasisSize> coefficients{};
  for (std::size_t column = 0; column < BasisSize; ++column) {
    coefficients[permutation[column]] = pivotedCoefficients[column];
  }

  // 7. Return row 0 as target values and rescale rows 1-3 into Cartesian gradients.
  CompactLinearFit result{
    .numberOfNeighbours = neighbourIndices.size(),
    .polynomialDegree = 1,
    .rank = BasisSize,
    .lengthScale = lengthScale
  };
  for (std::size_t field = 0; field < FieldCount; ++field) {
    result.values[field] = coefficients[0][field];
    for (std::size_t direction = 0; direction < 3; ++direction) {
      result.gradients[field][direction] =
        coefficients[direction + 1][field] / lengthScale;
    }
  }
  return result;
}

void InterpolateToSurfaceNearestNeighbour(
  const GridData& gridData,
  const AHGuess& ahGuess,
  const ReferencePoint& referencePoint,
  InterpolatedSurface& interpolatedSurface
) {
  InterpolateCartesianToSurface(
    gridData,
    ahGuess,
    referencePoint,
    interpolatedSurface,
    {.method = CartesianToSurfaceMethod::NearestNeighbour, .numberOfNeighbours = 1}
  );
}

void InterpolateCartesianToSurface(
  const GridData& gridData,
  const AHGuess& ahGuess,
  const ReferencePoint& referencePoint,
  InterpolatedSurface& interpolatedSurface,
  const CartesianToSurfaceOptions& options
) {
  // 1. Validate options and prepare a replacement output surface.
  if (options.numberOfNeighbours == 0) {
    throw std::invalid_argument(
      "InterpolateCartesianToSurface: neighbour count must be greater than zero");
  }

  std::vector<InterpolatedSurfacePoint> points;
  points.reserve(ahGuess.points.size());

  // 2. Reconstruct the Cartesian fields independently at each AH surface point.
  for (std::size_t pointIndex = 0; pointIndex < ahGuess.points.size(); ++pointIndex) {
    const AHGuessPoint& guessPoint = ahGuess.points[pointIndex];
    GridCoordinate target{
      std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::quiet_NaN()
    };
    try {
      // 3. Convert the surface point to Cartesian coordinates and select neighbours.
      target = SurfacePointToCartesian(guessPoint, referencePoint);
      const std::size_t neighbourCount =
        options.method == CartesianToSurfaceMethod::NearestNeighbour
        ? 1 : options.numberOfNeighbours;
      const std::vector<std::size_t> neighbours =
        FindNearestGridCells(gridData, target, neighbourCount);

      std::array<double, NumberOfCompactFields> compactValues{};
      std::array<Vector3, NumberOfCompactFields> compactGradients{};
      bool derivativesAvailable = false;

      // 4. Sample one cell directly, or reconstruct all fields with degree-1 MLS.
      if (options.method == CartesianToSurfaceMethod::NearestNeighbour) {
        const auto row = gridData.cell(neighbours.front()).row();
        std::copy(row.begin(), row.end(), compactValues.begin());
      } else {
        const LinearReconstructionOptions reconstructionOptions{
          .method = options.method == CartesianToSurfaceMethod::WeightedDegree1
            ? LinearReconstructionMethod::WeightedMovingLeastSquares
            : LinearReconstructionMethod::Unweighted,
          .numberOfNeighbours = neighbourCount,
          .weightSupport = options.weightSupport,
          .rankToleranceMultiplier = options.rankToleranceMultiplier
        };
        const CompactLinearFit fit = FitCompactLinearPolynomial(
          gridData, neighbours, target, reconstructionOptions);
        compactValues = fit.values;
        compactGradients = fit.gradients;
        derivativesAvailable = true;
      }

      // 5. Unpack the fixed compact layout into conformal tensors, K, and direct phi.
      std::array<double, NumberOfSymmetricTensorComponents> gammaTildeComponents{};
      std::array<double, NumberOfSymmetricTensorComponents> ATildeComponents{};
      for (std::size_t component = 0;
           component < NumberOfSymmetricTensorComponents;
           ++component) {
        gammaTildeComponents[component] = compactValues[component];
        ATildeComponents[component] = compactValues[component + 6];
      }
      const Tensor3 gammaTilde = UnpackSymmetricTensor(gammaTildeComponents);
      const Tensor3 ATilde = UnpackSymmetricTensor(ATildeComponents);
      const double traceK = compactValues[static_cast<std::size_t>(CompactField::TraceK)];
      const double phi = compactValues[static_cast<std::size_t>(CompactField::Phi)];

      // 6. Convert reconstructed conformal fields into physical geometry.
      InterpolatedSurfacePoint point{
        .r = guessPoint.r, .rho = guessPoint.rho, .sigma = guessPoint.sigma
      };
      point.geometry.spatialMetric = ComputePhysicalSpatialMetric(gammaTilde, phi);
      point.geometry.extrinsicCurvature =
        ComputePhysicalExtrinsicCurvature(ATilde, traceK, gammaTilde, phi);
      point.geometry.traceK = traceK;
      point.geometry.spatialMetricDerivativeAvailable = derivativesAvailable;

      // 7. MLS gradients become physical metric derivatives in [k][i][j] order.
      if (derivativesAvailable) {
        Tensor3Derivative gammaTildeDerivative{};
        for (std::size_t k = 0; k < 3; ++k) {
          std::array<double, NumberOfSymmetricTensorComponents> components{};
          for (std::size_t component = 0;
               component < NumberOfSymmetricTensorComponents;
               ++component) {
            components[component] = compactGradients[component][k];
          }
          gammaTildeDerivative[k] = UnpackSymmetricTensor(components);
        }
        point.geometry.spatialMetricDerivative =
          ComputePhysicalSpatialMetricDerivative(
            gammaTilde,
            gammaTildeDerivative,
            phi,
            compactGradients[static_cast<std::size_t>(CompactField::Phi)]
          );
      }

      // 8. Reject non-finite geometry before committing this surface point.
      const auto requireFinite = [pointIndex, &target](double value, const char* name) {
        if (!std::isfinite(value)) {
          std::ostringstream message;
          message << "InterpolateCartesianToSurface: surface point " << pointIndex
                  << " target (" << target.x << ", " << target.y << ", " << target.z
                  << ") produced non-finite " << name;
          throw std::runtime_error(message.str());
        }
      };
      requireFinite(point.geometry.traceK, "trace K");
      for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
          requireFinite(point.geometry.spatialMetric[i][j], "physical metric");
          requireFinite(point.geometry.extrinsicCurvature[i][j], "extrinsic curvature");
          for (std::size_t k = 0; k < 3; ++k) {
            requireFinite(
              point.geometry.spatialMetricDerivative[k][i][j], "metric derivative");
          }
        }
      }
      points.push_back(point);
    } catch (const std::exception& error) {
      // Attach the failing point and target to every lower-level error.
      std::ostringstream message;
      message << "InterpolateCartesianToSurface: surface point " << pointIndex
              << " target (" << target.x << ", " << target.y << ", " << target.z
              << ") failed: " << error.what();
      throw std::runtime_error(message.str());
    }
  }

  // Commit only after every point has been reconstructed successfully.
  interpolatedSurface.points = std::move(points);
}

void InterpolateFromCartesianToShell(const GridData& gridData,
                                     const AHGuess& ahGuess,
                                     const ReferencePoint& referencePoint,
                                     InterpolatedShell& interpolatedShell,
                                     std::span<const AngularSample> angularSamples) {
  (void)gridData;
  (void)ahGuess;
  (void)referencePoint;
  (void)angularSamples;
  (void)interpolatedShell;
  // Interpolate and calculate all necessary quantities from the 3D Cartesian grid to the 3D spherical
  //          coordinates using the given angular coordinates sample ($\{\rho_I,\sigma_I\}$) within a radius range
  //          (thickness) determined by AH guess input, i.e., we prepare a shell-like subdomain around
  //          the AH guess with spherical coordinates. The angular coordinate sample is given as input.
  // The radial coordinate sample can be determined by this function. Notice we need to calculate derivatives 
  //          here, so the subdomain may require halos if we do the process in parallel.
  // – Input: Full grid dataset ($Q\{x_I,y_I,z_I\}$), AH Guess function ($ \{h^{(k)}(\rho _I, \sigma _I)\}$),
  //          Reference point($s_0$), optional angular coordinate samples, Variable mapping (see below).
  //          If angularSamples is empty, use the angular coordinates stored in ahGuess.
  //– Output: the data structure InterpolatedShell that contains $T_{ij}$, $X^u_i$, $X^u_{ij}$, $ \gamma_{ij}$,
  //          $\partial_i \gamma_{jk}$, $K$ for each sample point ($\{r_I,\rho_I,\sigma_I\}$) within the
  //          shell-like subdomain, the structure also contains the shell outer and inner radii and the 
  //          angular coordinate sample.
  // The variable mapping input is required for generality; the user needs to specify how their quantities
  //          are ordered in the grid solution. The output quantities are calculated from the interpolated
  //          variables (the expressions of $T_{ij}$, $X^u_i$, and $X^u_{ij}$ are given below).
  // The data structure InterpolatedShell would be created (initialization) when this function is called
  //          for the first time during iteration. Then, during the following iteration (see below), the
  //          structure only needs to be updated in two cases: 
  //          1. the guess radius of the AH($h$) is approaching the radial boundary of the shell, 
  //             we need to thicken the shell in the corresponding direction.
  //          2. the resolution of the angular coordinates samples gets changed, we need to re-initialize
  //             the structure correspondingly.
  //          Except for these two cases, the structure remains unchanged.
  // We need derivatives of $Q$ over Cartesian coordinates and the derivative of $h$ over surface coordinates
  //          during the calculation.
}

void InterpolateFromShellToSurface(const InterpolatedShell& interpolatedShell,
                                   const AHGuess& ahGuess,
                                   const ReferencePoint& referencePoint,
                                   InterpolatedSurface& interpolatedSurface) {
  (void)interpolatedShell;
  (void)ahGuess;
  (void)referencePoint;
  (void)interpolatedSurface;
  // nterpolate and calculate all necessary quantities from the subdomain stored in InterpolatedShell
  //          to the surface points $\left\{\left( {h^{(k)}(\rho _I, \sigma _I)},\rho_I,\sigma_I \right)\right\}$.
  // Notice the angular sample needs to be identical, thus we actually only need to do the 1D radial
  //          interpolation and could increase the interpolation order for accuracy.
  // – Input: the data structure InterpolatedShell, AH Guess function($ \{h^{(k)}(\rho _I, \sigma _I)\}$),
  //          Reference point($s_0$).
  // – Output: the data structure InterpolatedSurface that contains $s_i$, $\partial_i s_j$, $T_{ij}$,
  //          $X^u_i$, $X^u_{ij}$, $ \gamma_{ij}$, $\partial_i \gamma_{jk}$, $K$ for each sample point
  //          ($\left\{\left( {h^{(k)}(\rho _I, \sigma _I)},\rho_I,\sigma_I \right)\right\}$) on the
  //          surface, and the angular coordinate sample
}

}  // namespace exahype2::ccz4::ahf
