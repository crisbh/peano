#ifndef APPARENT_HORIZON_FINDER_ITERATION_H_
#define APPARENT_HORIZON_FINDER_ITERATION_H_

#include "ahf-data.h"

#include <vector>

namespace exahype2::ccz4::ahf {

struct AHGuess;

/**
 * Solve the apparent-horizon PDE problem iteratively.
 */
void SolvePDEforAH(const AHGuess& initialAHGuess,
                   const ReferencePoint& referencePoint,
                   const std::vector<double>& masses,
                   AHGuess& solvedAHGuess);

/**
 * Determine whether a common-horizon search should be triggered for a binary-black-hole setup.
 */
void CheckCommonHorizonOfBBH(const AHGuess& initialAHGuess,
                             const ReferencePoint& referencePoint,
                             const std::vector<double>& masses,
                             bool& enableCommonHorizon,
                             ReferencePoint& commonReferencePoint,
                             double& commonMass,
                             AHGuess& commonAHGuesses);

}  // namespace exahype2::ccz4::ahf

#endif  // APPARENT_HORIZON_FINDER_ITERATION_H_
