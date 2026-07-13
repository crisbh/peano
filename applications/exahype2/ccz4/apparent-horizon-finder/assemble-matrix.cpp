#include "assemble-matrix.h"

namespace exahype2::ccz4::ahf {

void AssembleGlobalMatrix(const InterpolatedSurface& interpolatedSurface,
                          const AHGuess& ahGuess,
                          std::vector<std::vector<double>>& matrixCoefficients,
                          std::vector<double>& residuals) {
  (void)interpolatedSurface;
  (void)ahGuess;
  (void)matrixCoefficients;
  (void)residuals;
  // Use the quantities on the surface to calculate the expression of the local matrix,
  //     and the corresponding residuals of $\Theta$ over the sample points, then create
  //     the global sparse matrix.
  // – Input: the data structure InterpolatedSurface, AH Guess function ($ \{h^{(k)}(\rho _I, \sigma _I)\}$).
  // – Output: Global Matrix $A$ and Residuals $f$.
  // (We assume an expression of $A\left(h^{(k)}\right)h^{(k)}=f$ for Picard;
  //     $J\!\left[h^{(k)}\right]\delta h^{(k)} = -\Theta\!\left[h^{(k)}\right]$, followed by
  //     $h^{(k+1)} = h^{(k)} + \lambda \delta h^{(k)}$ for N-R.
}
}  // namespace exahype2::ccz4::ahf
