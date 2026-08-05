#include "interpolation.h"

#include <algorithm>
#include <numeric>

namespace exahype2::ccz4::ahf {

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

void InterpolateToSurfaceNearestNeighbour(
  const GridData& gridData,
  const AHGuess& ahGuess,
  const ReferencePoint& referencePoint,
  InterpolatedSurface& interpolatedSurface
) {
  std::vector<InterpolatedSurfacePoint> points;
  points.reserve(ahGuess.points.size());

  for (const AHGuessPoint& guessPoint : ahGuess.points) {
    const GridCoordinate target = SurfacePointToCartesian(guessPoint, referencePoint);
    const std::size_t cellIndex = FindNearestGridCells(gridData, target, 1).front();
    const CompactCCZ4CellView cell = gridData.cell(cellIndex);

    std::array<double, NumberOfSymmetricTensorComponents> gammaTildeComponents{};
    std::array<double, NumberOfSymmetricTensorComponents> ATildeComponents{};
    for (std::size_t component = 0;
         component < NumberOfSymmetricTensorComponents;
         ++component) {
      gammaTildeComponents[component] = cell.conformalMetricComponent(component);
      ATildeComponents[component] =
        cell.conformalExtrinsicCurvatureComponent(component);
    }

    const Tensor3 gammaTilde = UnpackSymmetricTensor(gammaTildeComponents);
    const Tensor3 ATilde = UnpackSymmetricTensor(ATildeComponents);
    const double traceK = cell.traceK();
    const double phi = cell.phi();

    InterpolatedSurfacePoint point{
      .r = guessPoint.r,
      .rho = guessPoint.rho,
      .sigma = guessPoint.sigma
    };
    point.geometry.spatialMetric = ComputePhysicalSpatialMetric(gammaTilde, phi);
    point.geometry.extrinsicCurvature =
      ComputePhysicalExtrinsicCurvature(ATilde, traceK, gammaTilde, phi);
    point.geometry.traceK = traceK;
    // Nearest-neighbour sampling supplies no meaningful Cartesian derivatives.
    point.geometry.spatialMetricDerivative = {};
    point.geometry.spatialMetricDerivativeAvailable = false;
    points.push_back(point);
  }

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
