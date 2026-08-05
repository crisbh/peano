#ifndef APPARENT_HORIZON_FINDER_INTERPOLATION_H_
#define APPARENT_HORIZON_FINDER_INTERPOLATION_H_

#include "ahf-data.h"

#include <span>
#include <vector>

namespace exahype2::ccz4::ahf {

/**
 * Convert one AH point (r, longitude rho, latitude sigma) into its global
 * Cartesian target coordinate relative to @p referencePoint.
 *
 * @throws std::invalid_argument for non-finite inputs or non-positive radius.
 */
GridCoordinate SurfacePointToCartesian(
  const AHGuessPoint& surfacePoint,
  const ReferencePoint& referencePoint
);

/**
 * Return the indices of the @p neighbourCount nearest Cartesian cell centres.
 * Results are ordered by increasing squared Euclidean distance, with cell index
 * as the deterministic tie-breaker. This brute-force baseline is deliberately
 * isolated behind this API so its implementation can later use a spatial index.
 *
 * @throws std::invalid_argument for an empty grid, an invalid neighbour count,
 * or a non-finite target coordinate.
 */
std::vector<std::size_t> FindNearestGridCells(
  const GridData& gridData,
  const GridCoordinate& target,
  std::size_t neighbourCount
);

struct LinearPolynomialFit {
  double value = 0.0;
  Vector3 derivative{};
  double lengthScale = 0.0;
};

/**
 * Unweighted degree-1 scalar least-squares reconstruction using the basis
 * [1,X,Y,Z] and a column-pivoted QR solve. The length scale ell is the maximum
 * Euclidean distance of any sample from the target, which bounds every scaled
 * sample radius by one and is insensitive to sample ordering.
 *
 * @throws std::invalid_argument for malformed/non-finite input, fewer than four
 * samples, or a geometrically rank-deficient sample distribution.
 */
LinearPolynomialFit FitLinearPolynomial(
  std::span<const GridCoordinate> sampleCoordinates,
  std::span<const double> sampleValues,
  const GridCoordinate& target
);

struct CompactLinearFit {
  std::array<double, NumberOfCompactFields> values{};
  std::array<Vector3, NumberOfCompactFields> gradients{};
  std::size_t numberOfNeighbours = 0;
  std::size_t polynomialDegree = 1;
  std::size_t rank = 0;
  double lengthScale = 0.0;
};

/**
 * Reconstruct all fourteen compact fields with one shared degree-1 design
 * matrix and one column-pivoted QR factorisation. Result fields retain compact
 * order, including K at 12 and direct phi at 13.
 */
CompactLinearFit FitCompactLinearPolynomial(
  const GridData& gridData,
  std::span<const std::size_t> neighbourIndices,
  const GridCoordinate& target
);

/**
 * Sample each AH surface point from its nearest Cartesian cell centre and
 * reconstruct the physical spatial metric and extrinsic curvature.
 *
 * This low-order baseline does not compute metric derivatives, coordinate
 * derivatives, or surface normals. Their storage is zero-initialised and
 * spatialMetricDerivativeAvailable is false for every output point.
 */
void InterpolateToSurfaceNearestNeighbour(
  const GridData& gridData,
  const AHGuess& ahGuess,
  const ReferencePoint& referencePoint,
  InterpolatedSurface& interpolatedSurface
);

/**
 * Interpolate the Cartesian grid data to a shell-like subdomain around the current AH guess.
 */
void InterpolateFromCartesianToShell(const GridData& gridData,
                                     const AHGuess& ahGuess,
                                     const ReferencePoint& referencePoint,
                                     InterpolatedShell& interpolatedShell,
                                     std::span<const AngularSample> angularSamples = {});

/**
 * Interpolate the shell data to the apparent-horizon surface defined by the current AH guess.
 */
void InterpolateFromShellToSurface(const InterpolatedShell& interpolatedShell,
                                   const AHGuess& ahGuess,
                                   const ReferencePoint& referencePoint,
                                   InterpolatedSurface& interpolatedSurface);

}  // namespace exahype2::ccz4::ahf

#endif  // APPARENT_HORIZON_FINDER_INTERPOLATION_H_
