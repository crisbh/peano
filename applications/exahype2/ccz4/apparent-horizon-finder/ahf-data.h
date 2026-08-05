#ifndef APPARENT_HORIZON_FINDER_AHF_DATA_H_
#define APPARENT_HORIZON_FINDER_AHF_DATA_H_

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
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
 * User-provided angular coordinate at which the shell should be sampled.
 */
struct AngularSample {
  double rho = 0.0;
  double sigma = 0.0;
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
 * Fixed column layout of every compact CCZ4 row stored in GridData.
 *
 * This enum is the only field-index contract used after file input. Its values
 * are independent of the component numbers in a source VTU tuple.
 */
enum class CompactField : std::size_t {
  GammaTildeXX = 0,
  GammaTildeXY = 1,
  GammaTildeXZ = 2,
  GammaTildeYY = 3,
  GammaTildeYZ = 4,
  GammaTildeZZ = 5,
  ATildeXX = 6,
  ATildeXY = 7,
  ATildeXZ = 8,
  ATildeYY = 9,
  ATildeYZ = 10,
  ATildeZZ = 11,
  TraceK = 12,
  Phi = 13,
  Count = 14
};

inline constexpr std::size_t NumberOfCompactCCZ4Fields =
  static_cast<std::size_t>(CompactField::Count);

inline constexpr std::array<const char*, NumberOfCompactCCZ4Fields> CompactFieldNames{
  "gammaTilde_xx", "gammaTilde_xy", "gammaTilde_xz",
  "gammaTilde_yy", "gammaTilde_yz", "gammaTilde_zz",
  "ATilde_xx", "ATilde_xy", "ATilde_xz",
  "ATilde_yy", "ATilde_yz", "ATilde_zz",
  "K", "phi"
};

/**
 * Reader-side mapping from compact CCZ4 fields to components of a source tuple.
 * Source indices are used only while constructing a compact GridData row.
 */
struct SourceVariableMapping {
  using Index = std::size_t;
  static constexpr Index NotStored = std::numeric_limits<Index>::max();

  std::array<Index, 6> conformalMetric{
    NotStored, NotStored, NotStored, NotStored, NotStored, NotStored
  };
  std::array<Index, 6> conformalExtrinsicCurvature{
    NotStored, NotStored, NotStored, NotStored, NotStored, NotStored
  };
  Index traceK = NotStored;
  Index conformalFactor = NotStored;

  std::array<Index, NumberOfCompactCCZ4Fields> componentsInCompactOrder() const {
    return {
      conformalMetric[0], conformalMetric[1], conformalMetric[2],
      conformalMetric[3], conformalMetric[4], conformalMetric[5],
      conformalExtrinsicCurvature[0], conformalExtrinsicCurvature[1],
      conformalExtrinsicCurvature[2], conformalExtrinsicCurvature[3],
      conformalExtrinsicCurvature[4], conformalExtrinsicCurvature[5],
      traceK, conformalFactor
    };
  }

  bool isComplete() const {
    const auto components = componentsInCompactOrder();
    return std::all_of(components.begin(), components.end(), [](Index index) {
      return index != NotStored;
    });
  }

  void validateRequiredFields() const {
    const auto components = componentsInCompactOrder();
    for (std::size_t field = 0; field < components.size(); ++field) {
      if (components[field] == NotStored) {
        throw std::invalid_argument(
          "SourceVariableMapping: required field "
          + std::string(CompactFieldNames[field]) + " is NotStored");
      }
    }

    for (std::size_t first = 0; first < components.size(); ++first) {
      for (std::size_t second = first + 1; second < components.size(); ++second) {
        if (components[first] == components[second]) {
          throw std::invalid_argument(
            "SourceVariableMapping: required fields "
            + std::string(CompactFieldNames[first]) + " and "
            + CompactFieldNames[second] + " both map to source component "
            + std::to_string(components[first]));
        }
      }
    }
  }

  void validateAgainstSourceTuple(std::size_t numberOfSourceComponents) const {
    validateRequiredFields();
    const auto components = componentsInCompactOrder();
    for (std::size_t field = 0; field < components.size(); ++field) {
      if (components[field] >= numberOfSourceComponents) {
        throw std::invalid_argument(
          "SourceVariableMapping: required field "
          + std::string(CompactFieldNames[field]) + " selects source component "
          + std::to_string(components[field]) + ", but the source tuple has only "
          + std::to_string(numberOfSourceComponents) + " components");
      }
    }
  }

  bool isValid(std::size_t numberOfSourceComponents) const {
    const auto components = componentsInCompactOrder();
    for (std::size_t field = 0; field < components.size(); ++field) {
      if (components[field] >= numberOfSourceComponents) {
        return false;
      }
      for (std::size_t previous = 0; previous < field; ++previous) {
        if (components[field] == components[previous]) {
          return false;
        }
      }
    }
    return true;
  }
};

/**
 * Apply the sole source-index to compact-index transformation. Both symmetric
 * tensors use xx, xy, xz, yy, yz, zz order; K and direct phi occupy columns 12
 * and 13 respectively.
 */
inline std::array<double, NumberOfCompactCCZ4Fields> CompactCCZ4SourceTuple(
  std::span<const double> sourceTuple,
  const SourceVariableMapping& sourceMapping
) {
  sourceMapping.validateAgainstSourceTuple(sourceTuple.size());

  const auto sourceComponents = sourceMapping.componentsInCompactOrder();

  std::array<double, NumberOfCompactCCZ4Fields> compact{};
  for (std::size_t column = 0; column < compact.size(); ++column) {
    compact[column] = sourceTuple[sourceComponents[column]];
  }
  return compact;
}

/**
 * Pointwise Cartesian grid data.
 *
 * Values are stored as a flat point-major array:
 *   values[pointIndex * numberOfVariables + variableIndex]
 */
struct GridData {
  std::vector<GridCoordinate> coordinates;
  std::vector<double> values;

  static constexpr std::size_t numberOfVariables = NumberOfCompactCCZ4Fields;

  std::size_t numberOfPoints() const {
    return coordinates.size();
  }

  void resize(std::size_t numberOfPoints) {
    coordinates.resize(numberOfPoints);
    values.resize(numberOfPoints * numberOfVariables);
  }

  bool isConsistent() const {
    return values.size() == coordinates.size() * numberOfVariables;
  }

  double& value(const GridCoordinate& coordinate, std::size_t variableIndex) {
    return values[valueIndex(coordinate, variableIndex)];
  }

  const double& value(const GridCoordinate& coordinate, std::size_t variableIndex) const {
    return values[valueIndex(coordinate, variableIndex)];
  }

  double& value(const GridCoordinate& coordinate, CompactField field) {
    return value(coordinate, static_cast<std::size_t>(field));
  }

  const double& value(const GridCoordinate& coordinate, CompactField field) const {
    return value(coordinate, static_cast<std::size_t>(field));
  }

private:
  std::size_t valueIndex(const GridCoordinate& coordinate, std::size_t variableIndex) const {
    if (variableIndex >= numberOfVariables) {
      throw std::out_of_range("GridData variable index is out of range");
    }

    const auto point = std::find_if(
      coordinates.begin(),
      coordinates.end(),
      [&coordinate](const GridCoordinate& candidate) {
        return coordinatesEqual(candidate.x, coordinate.x)
            && coordinatesEqual(candidate.y, coordinate.y)
            && coordinatesEqual(candidate.z, coordinate.z);
      }
    );
    if (point == coordinates.end()) {
      throw std::out_of_range("GridData coordinate was not found");
    }

    const auto pointIndex = static_cast<std::size_t>(point - coordinates.begin());
    return pointIndex * numberOfVariables + variableIndex;
  }

  static bool coordinatesEqual(double first, double second) {
    constexpr double AbsoluteTolerance = 1.0e-12;
    constexpr double RelativeTolerance = 1.0e-10;
    return std::abs(first - second)
        <= AbsoluteTolerance
         + RelativeTolerance * std::max(std::abs(first), std::abs(second));
  }
};

using Vector3 = std::array<double, 3>;
using Tensor3 = std::array<std::array<double, 3>, 3>;

/**
 * Physical and geometrical quantities shared by shell and surface samples.
 *
 * The first index of coordinateGradient and coordinateHessian selects the
 * surface coordinate u: 0 for rho and 1 for sigma.
 */
struct InterpolatedGeometryData {
  Tensor3 T{};
  std::array<Vector3, 2> coordinateGradient{};  // X^u_i
  std::array<Tensor3, 2> coordinateHessian{};   // X^u_ij
  Tensor3 spatialMetric{};                     // gamma_ij
  std::array<Tensor3, 3> spatialMetricDerivative{};  // partial_i gamma_jk
  Tensor3 extrinsicCurvature{};                // K_ij
  double traceK = 0.0;
};

/**
 * Interpolated quantities at one spherical shell sample (r, rho, sigma).
 */
struct InterpolatedShellPoint {
  double r = 0.0;
  double rho = 0.0;
  double sigma = 0.0;
  InterpolatedGeometryData geometry{};
};

/**
 * Interpolated data on the shell-like subdomain around an AH guess.
 *
 * Each entry in points owns both its spherical coordinate and all quantities
 * interpolated at that coordinate. Points belonging to one angular sample
 * should be stored consecutively and ordered by increasing r so that the
 * shell-to-surface step can perform one-dimensional radial interpolation.
 */
struct InterpolatedShell {
  double innerRadius = 0.0;
  double outerRadius = 0.0;
  std::vector<InterpolatedShellPoint> points;

  std::size_t numberOfPoints() const {
    return points.size();
  }

  InterpolatedShellPoint& value(const AHGuessPoint& coordinate) {
    return *findPoint(coordinate);
  }

  const InterpolatedShellPoint& value(const AHGuessPoint& coordinate) const {
    return *findPoint(coordinate);
  }

  bool isConsistent() const {
    if (!std::isfinite(innerRadius) || !std::isfinite(outerRadius)
        || innerRadius > outerRadius || points.empty()) {
      return false;
    }

    for (const InterpolatedShellPoint& point : points) {
      if (!std::isfinite(point.r) || !std::isfinite(point.rho)
          || !std::isfinite(point.sigma) || point.r < innerRadius
          || point.r > outerRadius) {
        return false;
      }
    }
    return true;
  }

private:
  std::vector<InterpolatedShellPoint>::iterator findPoint(const AHGuessPoint& coordinate) {
    const auto point = std::find_if(
      points.begin(),
      points.end(),
      [&coordinate](const InterpolatedShellPoint& candidate) {
        return coordinatesEqual(candidate.r, coordinate.r)
            && coordinatesEqual(candidate.rho, coordinate.rho)
            && coordinatesEqual(candidate.sigma, coordinate.sigma);
      }
    );
    if (point == points.end()) {
      throw std::out_of_range("InterpolatedShell coordinate was not found");
    }
    return point;
  }

  std::vector<InterpolatedShellPoint>::const_iterator findPoint(const AHGuessPoint& coordinate) const {
    const auto point = std::find_if(
      points.begin(),
      points.end(),
      [&coordinate](const InterpolatedShellPoint& candidate) {
        return coordinatesEqual(candidate.r, coordinate.r)
            && coordinatesEqual(candidate.rho, coordinate.rho)
            && coordinatesEqual(candidate.sigma, coordinate.sigma);
      }
    );
    if (point == points.end()) {
      throw std::out_of_range("InterpolatedShell coordinate was not found");
    }
    return point;
  }

  static bool coordinatesEqual(double first, double second) {
    constexpr double AbsoluteTolerance = 1.0e-12;
    constexpr double RelativeTolerance = 1.0e-10;
    return std::abs(first - second)
        <= AbsoluteTolerance
         + RelativeTolerance * std::max(std::abs(first), std::abs(second));
  }
};

/**
 * Interpolated quantities at one point on the current AH surface.
 */
struct InterpolatedSurfacePoint {
  double r = 0.0;
  double rho = 0.0;
  double sigma = 0.0;

  Vector3 surfaceNormal{};             // s_i
  Tensor3 surfaceNormalDerivative{};   // partial_i s_j
  InterpolatedGeometryData geometry{};
};

/**
 * Interpolated data on the current AH surface r = h(rho, sigma).
 *
 * Surface point i corresponds to AHGuess::points[i].
 */
struct InterpolatedSurface {
  std::vector<InterpolatedSurfacePoint> points;

  std::size_t numberOfPoints() const {
    return points.size();
  }

  InterpolatedSurfacePoint& value(const AHGuessPoint& coordinate) {
    return *findPoint(coordinate);
  }

  const InterpolatedSurfacePoint& value(const AHGuessPoint& coordinate) const {
    return *findPoint(coordinate);
  }

  bool isConsistent() const {
    if (points.empty()) {
      return false;
    }

    for (const InterpolatedSurfacePoint& point : points) {
      if (!std::isfinite(point.r) || !std::isfinite(point.rho)
          || !std::isfinite(point.sigma) || point.r < 0.0) {
        return false;
      }
    }
    return true;
  }

  bool matches(const AHGuess& ahGuess) const {
    if (points.size() != ahGuess.points.size()) {
      return false;
    }

    for (std::size_t i = 0; i < points.size(); ++i) {
      if (points[i].r != ahGuess.points[i].r
          || points[i].rho != ahGuess.points[i].rho
          || points[i].sigma != ahGuess.points[i].sigma) {
        return false;
      }
    }
    return true;
  }

private:
  std::vector<InterpolatedSurfacePoint>::iterator findPoint(const AHGuessPoint& coordinate) {
    const auto point = std::find_if(
      points.begin(),
      points.end(),
      [&coordinate](const InterpolatedSurfacePoint& candidate) {
        return coordinatesEqual(candidate.r, coordinate.r)
            && coordinatesEqual(candidate.rho, coordinate.rho)
            && coordinatesEqual(candidate.sigma, coordinate.sigma);
      }
    );
    if (point == points.end()) {
      throw std::out_of_range("InterpolatedSurface coordinate was not found");
    }
    return point;
  }

  std::vector<InterpolatedSurfacePoint>::const_iterator findPoint(const AHGuessPoint& coordinate) const {
    const auto point = std::find_if(
      points.begin(),
      points.end(),
      [&coordinate](const InterpolatedSurfacePoint& candidate) {
        return coordinatesEqual(candidate.r, coordinate.r)
            && coordinatesEqual(candidate.rho, coordinate.rho)
            && coordinatesEqual(candidate.sigma, coordinate.sigma);
      }
    );
    if (point == points.end()) {
      throw std::out_of_range("InterpolatedSurface coordinate was not found");
    }
    return point;
  }

  static bool coordinatesEqual(double first, double second) {
    constexpr double AbsoluteTolerance = 1.0e-12;
    constexpr double RelativeTolerance = 1.0e-10;
    return std::abs(first - second)
        <= AbsoluteTolerance
         + RelativeTolerance * std::max(std::abs(first), std::abs(second));
  }
};

}  // namespace exahype2::ccz4::ahf

#endif  // APPARENT_HORIZON_FINDER_AHF_DATA_H_
