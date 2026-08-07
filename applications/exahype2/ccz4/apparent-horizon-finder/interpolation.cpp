#include "interpolation.h"

#include <algorithm>
#include <array>
#include <limits>
#include <numeric>
#include <queue>
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

Tensor3 ComputeRadialProjectionTensor(const GridCoordinate& relativeCoordinate) {
  if (!std::isfinite(relativeCoordinate.x)
      || !std::isfinite(relativeCoordinate.y)
      || !std::isfinite(relativeCoordinate.z)) {
    throw std::invalid_argument(
      "ComputeRadialProjectionTensor: relative coordinate must be finite");
  }
  const Vector3 x{relativeCoordinate.x, relativeCoordinate.y, relativeCoordinate.z};
  const double radiusSquared = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
  Tensor3 T{};
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = i; j < 3; ++j) {
      T[i][j] = (i == j ? radiusSquared : 0.0) - x[i] * x[j];
      T[j][i] = T[i][j];
    }
  }
  return T;
}

AnalyticCoordinateGeometry ComputeAnalyticCoordinateGeometry(
  const GridCoordinate& globalCoordinate,
  const ReferencePoint& referencePoint
) {
  constexpr double MinimumAbsoluteCosSigma = 1.0e-12;
  if (!std::isfinite(globalCoordinate.x) || !std::isfinite(globalCoordinate.y)
      || !std::isfinite(globalCoordinate.z) || !std::isfinite(referencePoint.x)
      || !std::isfinite(referencePoint.y) || !std::isfinite(referencePoint.z)) {
    throw std::invalid_argument(
      "ComputeAnalyticCoordinateGeometry: coordinates must be finite");
  }
  const GridCoordinate relative{
    globalCoordinate.x - referencePoint.x,
    globalCoordinate.y - referencePoint.y,
    globalCoordinate.z - referencePoint.z
  };
  const double cylindricalRadius = std::hypot(relative.x, relative.y);
  const double radius = std::hypot(cylindricalRadius, relative.z);
  if (!(radius > 0.0) || !std::isfinite(radius)) {
    throw std::invalid_argument(
      "ComputeAnalyticCoordinateGeometry: radius must be finite and positive");
  }
  const double cosSigma = cylindricalRadius / radius;
  if (std::abs(cosSigma) < MinimumAbsoluteCosSigma) {
    throw std::domain_error(
      "ComputeAnalyticCoordinateGeometry: longitude coordinate is singular near cos(sigma)=0");
  }

  const double rho = std::atan2(relative.y, relative.x);
  const double sigma = std::atan2(relative.z, cylindricalRadius);
  const double sinRho = std::sin(rho);
  const double cosRho = std::cos(rho);
  const double sinSigma = std::sin(sigma);
  const double inverseRadius = 1.0 / radius;
  const double inverseRadiusSquared = inverseRadius * inverseRadius;
  const double inverseRadiusCosSigma = inverseRadius / cosSigma;
  const double inverseRadiusSquaredCosSigmaSquared =
    inverseRadiusSquared / (cosSigma * cosSigma);

  AnalyticCoordinateGeometry result;
  result.T = ComputeRadialProjectionTensor(relative);

  // X^u_i is stored as [u][i], with u=(rho,sigma) and i=(x,y,z).
  result.gradient[0] = {
    -sinRho * inverseRadiusCosSigma,
    cosRho * inverseRadiusCosSigma,
    0.0
  };
  result.gradient[1] = {
    -cosRho * sinSigma * inverseRadius,
    -sinRho * sinSigma * inverseRadius,
    cosSigma * inverseRadius
  };

  // X^u_ij is stored symmetrically as [u][i][j].
  Tensor3& rhoHessian = result.hessian[0];
  rhoHessian[0][0] = std::sin(2.0 * rho) * inverseRadiusSquaredCosSigmaSquared;
  rhoHessian[0][1] =
    -std::cos(2.0 * rho) * inverseRadiusSquaredCosSigmaSquared;
  rhoHessian[1][0] = rhoHessian[0][1];
  rhoHessian[1][1] = -rhoHessian[0][0];

  Tensor3& sigmaHessian = result.hessian[1];
  const double tanSigma = std::tan(sigma);
  const double cosTwoSigma = std::cos(2.0 * sigma);
  sigmaHessian[0][0] = tanSigma * inverseRadiusSquared
    * (2.0 * cosSigma * cosSigma * cosRho * cosRho - sinRho * sinRho);
  sigmaHessian[0][1] = 0.5 * std::sin(2.0 * rho) * inverseRadiusSquared
    * (std::sin(2.0 * sigma) + tanSigma);
  sigmaHessian[1][0] = sigmaHessian[0][1];
  sigmaHessian[0][2] = -cosRho * cosTwoSigma * inverseRadiusSquared;
  sigmaHessian[2][0] = sigmaHessian[0][2];
  sigmaHessian[1][1] = tanSigma * inverseRadiusSquared
    * (2.0 * cosSigma * cosSigma * sinRho * sinRho - cosRho * cosRho);
  sigmaHessian[1][2] = -sinRho * cosTwoSigma * inverseRadiusSquared;
  sigmaHessian[2][1] = sigmaHessian[1][2];
  sigmaHessian[2][2] = -std::sin(2.0 * sigma) * inverseRadiusSquared;
  return result;
}

std::vector<std::size_t> FindNearestGridCellsBruteForce(
  const GridData& gridData,
  const GridCoordinate& target,
  std::size_t neighbourCount
) {
  if (gridData.coordinates.empty()) {
    throw std::invalid_argument(
      "FindNearestGridCellsBruteForce: cannot search an empty Cartesian grid");
  }
  if (neighbourCount == 0) {
    throw std::invalid_argument(
      "FindNearestGridCellsBruteForce: neighbour count must be greater than zero");
  }
  if (neighbourCount > gridData.coordinates.size()) {
    throw std::invalid_argument(
      "FindNearestGridCellsBruteForce: requested " + std::to_string(neighbourCount)
      + " neighbours from only " + std::to_string(gridData.coordinates.size())
      + " cells");
  }
  if (!std::isfinite(target.x)
      || !std::isfinite(target.y)
      || !std::isfinite(target.z)) {
    throw std::invalid_argument(
      "FindNearestGridCellsBruteForce: target coordinate must be finite");
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

struct CartesianSpatialIndex::Implementation {
  // 1. Store one grid cell and the links/splitting axis of its k-d tree node.
  struct Node {
    std::size_t cellIndex = 0;
    std::size_t left = std::numeric_limits<std::size_t>::max();
    std::size_t right = std::numeric_limits<std::size_t>::max();
    unsigned char axis = 0;
  };

  static constexpr std::size_t NoNode = std::numeric_limits<std::size_t>::max();

  // 2. Validate and copy coordinates, then build the immutable tree once.
  explicit Implementation(const GridData& gridData): coordinates(gridData.coordinates) {
    if (coordinates.empty()) {
      throw std::invalid_argument("CartesianSpatialIndex: cannot search an empty grid");
    }
    for (const GridCoordinate& coordinate : coordinates) {
      if (!std::isfinite(coordinate.x) || !std::isfinite(coordinate.y)
          || !std::isfinite(coordinate.z)) {
        throw std::invalid_argument("CartesianSpatialIndex: coordinates must be finite");
      }
    }
    std::vector<std::size_t> indices(coordinates.size());
    std::iota(indices.begin(), indices.end(), std::size_t{0});
    nodes.reserve(coordinates.size());
    root = build(indices, 0, indices.size(), 0);
  }

  // Read one Cartesian component selected by axis=(x,y,z).
  double component(std::size_t cellIndex, std::size_t axis) const {
    const GridCoordinate& coordinate = coordinates[cellIndex];
    return axis == 0 ? coordinate.x : axis == 1 ? coordinate.y : coordinate.z;
  }

  // 3. Recursively median-split the cells, cycling through x, y, and z.
  std::size_t build(
    std::vector<std::size_t>& indices,
    std::size_t begin,
    std::size_t end,
    std::size_t depth
  ) {
    if (begin == end) {
      return NoNode;
    }
    const std::size_t axis = depth % 3;
    const std::size_t middle = begin + (end - begin) / 2;
    std::nth_element(
      indices.begin() + begin,
      indices.begin() + middle,
      indices.begin() + end,
      [this, axis](std::size_t first, std::size_t second) {
        const double firstCoordinate = component(first, axis);
        const double secondCoordinate = component(second, axis);
        return firstCoordinate < secondCoordinate
            || (firstCoordinate == secondCoordinate && first < second);
      }
    );
    const std::size_t nodeIndex = nodes.size();
    nodes.push_back(Node{.cellIndex = indices[middle], .axis = static_cast<unsigned char>(axis)});
    const std::size_t left = build(indices, begin, middle, depth + 1);
    const std::size_t right = build(indices, middle + 1, end, depth + 1);
    nodes[nodeIndex].left = left;
    nodes[nodeIndex].right = right;
    return nodeIndex;
  }

  // 4. Compare candidates without an unnecessary square root.
  double squaredDistance(std::size_t cellIndex, const GridCoordinate& target) const {
    const GridCoordinate& coordinate = coordinates[cellIndex];
    const double dx = coordinate.x - target.x;
    const double dy = coordinate.y - target.y;
    const double dz = coordinate.z - target.z;
    return dx * dx + dy * dy + dz * dz;
  }

  using Candidate = std::pair<double, std::size_t>;
  using CandidateHeap = std::priority_queue<Candidate>;

  // 5. Search the near branch first and prune the far branch when possible.
  void query(
    std::size_t nodeIndex,
    const GridCoordinate& target,
    std::size_t neighbourCount,
    CandidateHeap& candidates
  ) const {
    if (nodeIndex == NoNode) {
      return;
    }
    const Node& node = nodes[nodeIndex];
    const Candidate candidate{squaredDistance(node.cellIndex, target), node.cellIndex};
    if (candidates.size() < neighbourCount) {
      candidates.push(candidate);
    } else if (candidate < candidates.top()) {
      candidates.pop();
      candidates.push(candidate);
    }

    const double targetComponent = node.axis == 0 ? target.x
      : node.axis == 1 ? target.y : target.z;
    const double difference = targetComponent - component(node.cellIndex, node.axis);
    const std::size_t nearChild = difference <= 0.0 ? node.left : node.right;
    const std::size_t farChild = difference <= 0.0 ? node.right : node.left;
    query(nearChild, target, neighbourCount, candidates);
    const double worstDistance = candidates.size() < neighbourCount
      ? std::numeric_limits<double>::infinity() : candidates.top().first;
    if (difference * difference <= worstDistance) {
      query(farChild, target, neighbourCount, candidates);
    }
  }

  // Immutable storage shared by all read-only queries.
  std::vector<GridCoordinate> coordinates;
  std::vector<Node> nodes;
  std::size_t root = NoNode;
};

CartesianSpatialIndex::CartesianSpatialIndex(const GridData& gridData):
  _implementation(std::make_shared<const Implementation>(gridData)) {}

std::size_t CartesianSpatialIndex::size() const {
  return _implementation->coordinates.size();
}

std::vector<std::size_t> CartesianSpatialIndex::findNearest(
  const GridCoordinate& target,
  std::size_t neighbourCount
) const {
  if (!std::isfinite(target.x) || !std::isfinite(target.y)
      || !std::isfinite(target.z)) {
    throw std::invalid_argument("CartesianSpatialIndex: target must be finite");
  }
  if (neighbourCount == 0) {
    throw std::invalid_argument(
      "CartesianSpatialIndex: neighbour count must be greater than zero");
  }
  if (neighbourCount > size()) {
    throw std::invalid_argument(
      "CartesianSpatialIndex: requested " + std::to_string(neighbourCount)
      + " neighbours from only " + std::to_string(size()) + " cells");
  }
  Implementation::CandidateHeap candidates;
  _implementation->query(_implementation->root, target, neighbourCount, candidates);
  std::vector<Implementation::Candidate> ordered;
  ordered.reserve(candidates.size());
  while (!candidates.empty()) {
    ordered.push_back(candidates.top());
    candidates.pop();
  }
  std::sort(ordered.begin(), ordered.end());
  std::vector<std::size_t> indices;
  indices.reserve(ordered.size());
  for (const auto& candidate : ordered) {
    indices.push_back(candidate.second);
  }
  return indices;
}

std::vector<std::size_t> FindNearestGridCells(
  const CartesianSpatialIndex& spatialIndex,
  const GridCoordinate& target,
  std::size_t neighbourCount
) {
  return spatialIndex.findNearest(target, neighbourCount);
}

std::vector<std::size_t> FindNearestGridCells(
  const GridData& gridData,
  const GridCoordinate& target,
  std::size_t neighbourCount
) {
  return CartesianSpatialIndex(gridData).findNearest(target, neighbourCount);
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
  double minimumDiagonal = std::numeric_limits<double>::infinity();
  double maximumDiagonal = 0.0;
  for (std::size_t diagonal = 0; diagonal < BasisSize; ++diagonal) {
    const double magnitude = std::abs(upperTriangular[diagonal][diagonal]);
    minimumDiagonal = std::min(minimumDiagonal, magnitude);
    maximumDiagonal = std::max(maximumDiagonal, magnitude);
  }
  CompactLinearFit result{
    .numberOfNeighbours = neighbourIndices.size(),
    .polynomialDegree = 1,
    .rank = BasisSize,
    .lengthScale = lengthScale,
    .conditionEstimate = maximumDiagonal / minimumDiagonal
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
  if (gridData.coordinates.empty()) {
    throw std::invalid_argument(
      "InterpolateCartesianToSurface: Cartesian grid must not be empty");
  }
  if (!gridData.isConsistent()) {
    throw std::invalid_argument(
      "InterpolateCartesianToSurface: GridData value count must equal cell count times fourteen");
  }
  if (!std::isfinite(options.domainRelativeTolerance)
      || options.domainRelativeTolerance < 0.0) {
    throw std::invalid_argument(
      "InterpolateCartesianToSurface: domain relative tolerance must be finite and non-negative");
  }

  GridCoordinate domainMinimum = gridData.coordinates.front();
  GridCoordinate domainMaximum = gridData.coordinates.front();
  for (const GridCoordinate& coordinate : gridData.coordinates) {
    if (!std::isfinite(coordinate.x) || !std::isfinite(coordinate.y)
        || !std::isfinite(coordinate.z)) {
      throw std::invalid_argument(
        "InterpolateCartesianToSurface: grid coordinates must be finite");
    }
    domainMinimum.x = std::min(domainMinimum.x, coordinate.x);
    domainMinimum.y = std::min(domainMinimum.y, coordinate.y);
    domainMinimum.z = std::min(domainMinimum.z, coordinate.z);
    domainMaximum.x = std::max(domainMaximum.x, coordinate.x);
    domainMaximum.y = std::max(domainMaximum.y, coordinate.y);
    domainMaximum.z = std::max(domainMaximum.z, coordinate.z);
  }
  const double domainDiagonal = std::hypot(
    domainMaximum.x - domainMinimum.x,
    domainMaximum.y - domainMinimum.y,
    domainMaximum.z - domainMinimum.z
  );
  // Build once; all surface-point queries below are immutable and thread-safe.
  const CartesianSpatialIndex spatialIndex(gridData);

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
      const double outsideX = std::max(
        {domainMinimum.x - target.x, 0.0, target.x - domainMaximum.x});
      const double outsideY = std::max(
        {domainMinimum.y - target.y, 0.0, target.y - domainMaximum.y});
      const double outsideZ = std::max(
        {domainMinimum.z - target.z, 0.0, target.z - domainMaximum.z});
      const double outsideDistance = std::hypot(outsideX, outsideY, outsideZ);
      if (outsideDistance > options.domainRelativeTolerance * domainDiagonal) {
        throw std::out_of_range(
          "target lies outside the sampled Cartesian domain tolerance");
      }

      InterpolationDiagnostics diagnostics{
        .requestedMethod = options.method,
        .requestedDegree = options.method == CartesianToSurfaceMethod::NearestNeighbour ? 0u : 1u,
        .degreeUsed = options.method == CartesianToSurfaceMethod::NearestNeighbour ? 0u : 1u,
        .neighbourCount = options.method == CartesianToSurfaceMethod::NearestNeighbour
          ? 1u : options.numberOfNeighbours,
        .target = target,
        .rank = options.method == CartesianToSurfaceMethod::NearestNeighbour ? 1u : 0u,
        .conditionEstimate = 1.0,
        .orderReduced = false
      };

      bool useNearestNeighbour =
        options.method == CartesianToSurfaceMethod::NearestNeighbour;
      if (!useNearestNeighbour
          && (options.numberOfNeighbours < 4
              || options.numberOfNeighbours > gridData.coordinates.size())) {
        if (!options.allowOrderReduction) {
          throw std::invalid_argument(
            "degree-1 reconstruction has insufficient neighbours and order reduction is disabled");
        }
        useNearestNeighbour = true;
        diagnostics.degreeUsed = 0;
        diagnostics.neighbourCount = 1;
        diagnostics.rank = 1;
        diagnostics.orderReduced = true;
      }

      std::vector<std::size_t> neighbours = FindNearestGridCells(
        spatialIndex, target, useNearestNeighbour ? 1 : options.numberOfNeighbours);

      const auto validateStencilValues = [&gridData](std::span<const std::size_t> stencil) {
        for (std::size_t cellIndex : stencil) {
          const auto row = gridData.cell(cellIndex).row();
          for (double value : row) {
            if (!std::isfinite(value)) {
              throw std::invalid_argument(
                "selected interpolation stencil contains a non-finite compact value");
            }
          }
        }
      };
      validateStencilValues(neighbours);

      std::array<double, NumberOfCompactFields> compactValues{};
      std::array<Vector3, NumberOfCompactFields> compactGradients{};
      bool derivativesAvailable = false;

      // 4. Sample one cell directly, or reconstruct all fields with degree-1 MLS.
      if (useNearestNeighbour) {
        const auto row = gridData.cell(neighbours.front()).row();
        std::copy(row.begin(), row.end(), compactValues.begin());
      } else {
        const LinearReconstructionOptions reconstructionOptions{
          .method = options.method == CartesianToSurfaceMethod::WeightedDegree1
            ? LinearReconstructionMethod::WeightedMovingLeastSquares
            : LinearReconstructionMethod::Unweighted,
          .numberOfNeighbours = neighbours.size(),
          .weightSupport = options.weightSupport,
          .rankToleranceMultiplier = options.rankToleranceMultiplier
        };
        try {
          const CompactLinearFit fit = FitCompactLinearPolynomial(
            gridData, neighbours, target, reconstructionOptions);
          compactValues = fit.values;
          compactGradients = fit.gradients;
          derivativesAvailable = true;
          diagnostics.rank = fit.rank;
          diagnostics.conditionEstimate = fit.conditionEstimate;
        } catch (const std::invalid_argument& error) {
          const std::string message = error.what();
          const bool rankFailure = message.find("rank deficient") != std::string::npos;
          if (!rankFailure || !options.allowOrderReduction) {
            throw;
          }
          neighbours = FindNearestGridCells(spatialIndex, target, 1);
          validateStencilValues(neighbours);
          const auto row = gridData.cell(neighbours.front()).row();
          std::copy(row.begin(), row.end(), compactValues.begin());
          diagnostics.degreeUsed = 0;
          diagnostics.neighbourCount = 1;
          diagnostics.rank = 1;
          diagnostics.conditionEstimate = 1.0;
          diagnostics.orderReduced = true;
        }
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
      point.interpolationDiagnostics = diagnostics;
      const AnalyticCoordinateGeometry coordinateGeometry =
        ComputeAnalyticCoordinateGeometry(target, referencePoint);
      point.geometry.T = coordinateGeometry.T;
      point.geometry.coordinateGradient = coordinateGeometry.gradient;
      point.geometry.coordinateHessian = coordinateGeometry.hessian;
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
