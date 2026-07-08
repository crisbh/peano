#ifndef APPARENT_HORIZON_FINDER_AHF_DATA_H_
#define APPARENT_HORIZON_FINDER_AHF_DATA_H_

#include <cstddef>
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

}  // namespace exahype2::ccz4::ahf

#endif  // APPARENT_HORIZON_FINDER_AHF_DATA_H_
