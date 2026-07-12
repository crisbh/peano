#include "iteration.h"

namespace exahype2::ccz4::ahf {

void SolvePDEforAH(const AHGuess& initialAHGuess,
                   const ReferencePoint& referencePoint,
                   const std::vector<double>& masses,
                   AHGuess& solvedAHGuess) {
  (void)initialAHGuess;
  (void)referencePoint;
  (void)masses;
  (void)solvedAHGuess;
  // Call InterpolateFromCartesianToShell() [though this function is called in every iteration,
  //      it would not do anything except the very first iteration and the two cases state above],
  //      InterpolateFromShellToSurface() and AssembleGlobalMatrix(), then determine the new AH
  //      guess function ($ \{h^{(k+1)}(\rho _I, \sigma _I)\}$) (say, use Newton-Raphson). Call
  //      the functions again. Repeat until convergence using the mass-scaled threshold.
  // Notice we are on a 2D surface, so no boundary condition is required (natural periodic BC).
  // We search the AH one by one (i.e., try the first set of a reference point, black hole mass,
  //      and initial AH guess; then move to the next set).
  // We may also consider Picard iteration, as the physics quantities near the AH should change
  //      smoothly enough.
  // – Input: Initial AH Guess functions ($ \{h_{i}^{(0)}(\rho _I, \sigma _I)\}$), Reference
  //      points($\{s_i\}$), The black hole mass associated with corresponding reference points.
  // – Output: Solution AH function ($ \{h_i^{(f)}(\rho _I, \sigma _I)\}$).
  // This function can accept multiple reference points and the corresponding initial guesses.
  //     The threshold of iteration stop would depend on the input black hole mass (scaled by that).
  // Failing control: If one initial guess fails, we trigger a two-step scheme:
  // 1. The resolution on the surface may not be good enough, so we should refine our sample points
  //    and call the interpolation function to initialize the data structure again. We would always
  //    do this refine-retry first when an AH search fails. However, we set an upper limit of
  //    resolution here, and if we still fail when the limit is reached, we stop and switch to step 2.
  // 2. This initial guess is wrong (including reference point choice and/or initial guess on $h$),
  //    so we move to the next set of initial guess.
}

void CheckCommonHorizonOfBBH(const AHGuess& initialAHGuess,
                             const ReferencePoint& referencePoint,
                             const std::vector<double>& masses,
                             bool& enableCommonHorizon,
                             ReferencePoint& commonReferencePoint,
                             double& commonMass,
                             AHGuess& commonAHGuesses) {
  (void)initialAHGuess;
  (void)referencePoint;
  (void)masses;
  (void)enableCommonHorizon;
  (void)commonReferencePoint;
  (void)commonMass;
  (void)commonAHGuesses;
  //used in SolvePDEforAH() for pre-conditioning for binary black hole scenarios. It calculates
  //     the coordinate distance between the two reference points and compares it to the sum of
  //     the apparent horizon sizes ($h_{1,2}^{(0)}$) in the corresponding direction. If the latter
  //     is larger, we tell SolvePDEforAH() to solve AH first using the weight center of the two
  //     black holes as a reference point, and provide it with associated mass and new initial
  //     guess functions.
  // – Input: Initial AH Guess functions ($ \{h_{1,2}^{(0)}(\rho _I, \sigma _I)\}$), Reference
  //      points($\{s_1,s_2\}$), The black hole mass associated with corresponding reference
  //      points($\{M_1, M_2\}$).
  // – Output: Flag to enable common horizon search, new reference point($s_{com}$), mass
  //      scale($M_{com}$), and Initial AH Guess functions ($ \{h_{com}^{(0)}(\rho _I, \sigma _I)\}$)
}

}  // namespace exahype2::ccz4::ahf
