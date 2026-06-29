#include "interpolation.h"

namespace exahype2::ccz4::ahf {

void InterpolateFromCartesianToShell(const GridData& gridData,
                                     const std::vector<AHGuess>& ahGuess,
                                     const ReferencePoint& referencePoint,
                                     const std::vector<double>& angularSamples,
                                     const VariableMapping& variableMapping,
                                     InterpolatedShell& interpolatedShell) {
  (void)gridData;
  (void)ahGuess;
  (void)referencePoint;
  (void)angularSamples;
  (void)variableMapping;
  (void)interpolatedShell;
  // Interpolate and calculate all necessary quantities from the 3D Cartesian grid to the 3D spherical
  //          coordinates using the given angular coordinates sample ($\{\rho_I,\sigma_I\}$) within a radius range
  //          (thickness) determined by AH guess input, i.e., we prepare a shell-like subdomain around
  //          the AH guess with spherical coordinates. The angular coordinate sample is given as input.
  // The radial coordinate sample can be determined by this function. Notice we need to calculate derivatives 
  //          here, so the subdomain may require halos if we do the process in parallel.
  // – Input: Full grid dataset ($Q\{x_I,y_I,z_I\}$), AH Guess function ($ \{h^{(k)}(\rho _I, \sigma _I)\}$),
  //          Reference point($s_0$). angular coordinate samples, Variable mapping (see below).
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
                                   const std::vector<AHGuess>& ahGuess,
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
