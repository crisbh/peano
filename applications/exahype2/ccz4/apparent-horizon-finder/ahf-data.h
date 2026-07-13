#ifndef APPARENT_HORIZON_FINDER_AHF_DATA_H_
#define APPARENT_HORIZON_FINDER_AHF_DATA_H_

#include <array>
#include <cstddef>
#include <limits>
#include <vector>

namespace exahype2::ccz4::ahf {

/**
 * Cartesian reference coordinate for one apparent-horizon centre point.
 */
struct ReferencePoint {
  double x;
  double y;
  double z;
};

/**
 * One point on an apparent-horizon guess, expressed in spherical coordinates.
 */
struct AHGuessPoint {
  double r;
  double rho;
  double sigma;
};

/**
 * One apparent-horizon hypersurface guess.
 *
 * Multiple candidate guesses are represented by std::vector<AHGuess>.
 */
struct AHGuess {
  std::vector<AHGuessPoint> points;
};

/**
 * Cartesian coordinate for one sampled grid point.
 */
struct GridCoordinate {
  double x;
  double y;
  double z;
};

/**
 * Pointwise Cartesian grid data.
 *
 * Values are stored as a flat point-major array:
 *   values[pointIndex * numberOfVariables + variableIndex]
 */
struct GridData {
  std::vector<GridCoordinate> coordinates;
  std::vector<double> values;
  std::size_t numberOfVariables = 0;

  std::size_t numberOfPoints() const {
    return coordinates.size();
  }

  void setNumberOfVariables(std::size_t newNumberOfVariables) {
    numberOfVariables = newNumberOfVariables;
    values.resize(coordinates.size() * numberOfVariables);
  }

  bool isConsistent() const {
    return numberOfVariables > 0 && values.size() == coordinates.size() * numberOfVariables;
  }

  double& value(std::size_t pointIndex, std::size_t variableIndex) {
    return values[pointIndex * numberOfVariables + variableIndex];
  }

  const double& value(std::size_t pointIndex, std::size_t variableIndex) const {
    return values[pointIndex * numberOfVariables + variableIndex];
  }
};

/**
 * User-provided mapping from AHF quantities to variable indices in GridData.
 *
 * Symmetric tensors use the component order xx, xy, xz, yy, yz, zz. The
 * mapping contains indices only; the corresponding values remain in
 * GridData::values.
 */
struct VariableMapping {
  using Index = std::size_t;

  static constexpr Index NotStored = std::numeric_limits<Index>::max();

  enum class ConformalFactorEncoding {
    Phi,
    LogPhi
  };

  std::array<Index, 6> conformalMetric{
    NotStored, NotStored, NotStored, NotStored, NotStored, NotStored
  };
  std::array<Index, 6> conformalExtrinsicCurvature{
    NotStored, NotStored, NotStored, NotStored, NotStored, NotStored
  };
  Index traceK = NotStored;
  Index conformalFactor = NotStored;
  ConformalFactorEncoding conformalFactorEncoding = ConformalFactorEncoding::Phi;

  bool isValid(std::size_t numberOfVariables) const {
    for (Index index : conformalMetric) {
      if (index >= numberOfVariables) {
        return false;
      }
    }
    for (Index index : conformalExtrinsicCurvature) {
      if (index >= numberOfVariables) {
        return false;
      }
    }
    return traceK < numberOfVariables && conformalFactor < numberOfVariables;
  }
};

}  // namespace exahype2::ccz4::ahf

#endif  // APPARENT_HORIZON_FINDER_AHF_DATA_H_
