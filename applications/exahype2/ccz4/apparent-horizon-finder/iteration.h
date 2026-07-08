#ifndef APPARENT_HORIZON_FINDER_ITERATION_H_
#define APPARENT_HORIZON_FINDER_ITERATION_H_

#include "ahf-data.h"

#include <vector>

namespace exahype2::ccz4::ahf {

struct AHGuess;

/**
 * Solve the apparent-horizon PDE problem iteratively.
 */
void SolvePDEforAH(const std::vector<AHGuess>& initialAHGuesses,
                   const std::vector<ReferencePoint>& referencePoints,
                   const std::vector<double>& masses,
                   std::vector<AHGuess>& solvedAHGuesses);

/**
 * Determine whether a common-horizon search should be triggered for a binary-black-hole setup.
 */
void CheckCommonHorizonOfBBH(const std::vector<AHGuess>& initialAHGuesses,
                             const std::vector<ReferencePoint>& referencePoints,
                             const std::vector<double>& masses,
                             bool& enableCommonHorizon,
                             ReferencePoint& commonReferencePoint,
                             double& commonMass,
                             std::vector<AHGuess>& commonAHGuesses);

}  // namespace exahype2::ccz4::ahf

#endif  // APPARENT_HORIZON_FINDER_ITERATION_H_
