#ifndef APPARENT_HORIZON_FINDER_INTERPOLATION_H_
#define APPARENT_HORIZON_FINDER_INTERPOLATION_H_

#include "ahf-data.h"

#include <vector>

namespace exahype2::ccz4::ahf {

struct InterpolatedShell;
struct InterpolatedSurface;
struct AHGuess;

/**
 * Interpolate the Cartesian grid data to a shell-like subdomain around the current AH guess.
 */
void InterpolateFromCartesianToShell(const GridData& gridData,
                                     const std::vector<AHGuess>& ahGuess,
                                     const ReferencePoint& referencePoint,
                                     const std::vector<double>& angularSamples,
                                     const VariableMapping& variableMapping,
                                     InterpolatedShell& interpolatedShell);

/**
 * Interpolate the shell data to the apparent-horizon surface defined by the current AH guess.
 */
void InterpolateFromShellToSurface(const InterpolatedShell& interpolatedShell,
                                   const std::vector<AHGuess>& ahGuess,
                                   const ReferencePoint& referencePoint,
                                   InterpolatedSurface& interpolatedSurface);

}  // namespace exahype2::ccz4::ahf

#endif  // APPARENT_HORIZON_FINDER_INTERPOLATION_H_
