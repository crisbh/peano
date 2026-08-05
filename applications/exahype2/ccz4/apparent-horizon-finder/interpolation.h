#ifndef APPARENT_HORIZON_FINDER_INTERPOLATION_H_
#define APPARENT_HORIZON_FINDER_INTERPOLATION_H_

#include "ahf-data.h"

#include <span>

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
