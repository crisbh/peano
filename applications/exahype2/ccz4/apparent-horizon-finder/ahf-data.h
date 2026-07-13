#ifndef APPARENT_HORIZON_FINDER_AHF_DATA_H_
#define APPARENT_HORIZON_FINDER_AHF_DATA_H_

#include <array>
#include <cmath>
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
};

}  // namespace exahype2::ccz4::ahf

#endif  // APPARENT_HORIZON_FINDER_AHF_DATA_H_
