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

struct AnalyticCoordinateGeometry {
  Tensor3 T{};
  /** [u][i], with u=0 for rho, u=1 for sigma, and i=x,y,z. */
  std::array<Vector3, 2> gradient{};
  /** [u][i][j] = partial_i partial_j u, with u=0 rho and u=1 sigma. */
  std::array<Tensor3, 2> hessian{};
};

/** Compute T_ij = r^2 delta_ij - x_i x_j from relative coordinates. */
Tensor3 ComputeRadialProjectionTensor(const GridCoordinate& relativeCoordinate);

/**
 * Compute the appendix's analytic rho/sigma coordinate geometry. Global input
 * is translated by @p referencePoint before evaluating every formula.
 * Longitude is undefined near |cos(sigma)|=0, so poles are rejected rather
 * than assigned an artificial treatment.
 */
AnalyticCoordinateGeometry ComputeAnalyticCoordinateGeometry(
  const GridCoordinate& globalCoordinate,
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

enum class LinearReconstructionMethod {
  Unweighted,
  WeightedMovingLeastSquares
};

struct LinearReconstructionOptions {
  LinearReconstructionMethod method = LinearReconstructionMethod::Unweighted;
  /** Zero uses every supplied sample; otherwise use the first requested samples. */
  std::size_t numberOfNeighbours = 0;
  /** Gaussian support alpha in exp(-(distance/(alpha*ell))^2). */
  double weightSupport = 1.0;
  /**
   * Multiplier applied to epsilon*n*largest-column-norm for rank detection.
   * Controls how strictly the QR solver decides that the neighbour geometry
   * is rank deficient.
   */
  double rankToleranceMultiplier = 64.0;
};

/** Evaluate the isolated unweighted/Gaussian row weight. */
double LinearReconstructionWeight(
  double distance,
  double lengthScale,
  const LinearReconstructionOptions& options
);

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
  const GridCoordinate& target,
  const LinearReconstructionOptions& options = {}
);

struct CompactLinearFit {
  std::array<double, NumberOfCompactFields> values{};
  std::array<Vector3, NumberOfCompactFields> gradients{};
  std::size_t numberOfNeighbours = 0;
  std::size_t polynomialDegree = 1;
  std::size_t rank = 0;
  double lengthScale = 0.0;
  /** Approximate max/min absolute diagonal ratio of the pivoted R factor. */
  double conditionEstimate = 1.0;
};

/**
 * Reconstruct all fourteen compact fields with one shared degree-1 design
 * matrix and one column-pivoted QR factorisation. Result fields retain compact
 * order, including K at 12 and direct phi at 13.
 */
CompactLinearFit FitCompactLinearPolynomial(
  const GridData& gridData,
  std::span<const std::size_t> neighbourIndices,
  const GridCoordinate& target,
  const LinearReconstructionOptions& options = {}
);

struct CartesianToSurfaceOptions {
  CartesianToSurfaceMethod method = CartesianToSurfaceMethod::WeightedDegree1;
  std::size_t numberOfNeighbours = 8;
  double weightSupport = 1.0;
  double rankToleranceMultiplier = 64.0;
  bool allowOrderReduction = false;
  /** Allowed distance outside the grid AABB as a fraction of its diagonal. */
  double domainRelativeTolerance = 0.05;
};

/**
 * Interpolate compact Cartesian fields directly to every AH surface point.
 * Degree-1 modes populate physical metric derivatives in [k][i][j] order,
 * where k is the Cartesian differentiation direction. Surface-shape and
 * angular derivatives are intentionally not computed.
 */
void InterpolateCartesianToSurface(
  const GridData& gridData,
  const AHGuess& ahGuess,
  const ReferencePoint& referencePoint,
  InterpolatedSurface& interpolatedSurface,
  const CartesianToSurfaceOptions& options = {}
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
