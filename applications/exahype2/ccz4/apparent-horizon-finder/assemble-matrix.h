#ifndef APPARENT_HORIZON_FINDER_ASSEMBLE_MATRIX_H_
#define APPARENT_HORIZON_FINDER_ASSEMBLE_MATRIX_H_

#include <vector>

namespace exahype2::ccz4::ahf {

struct InterpolatedSurface;
struct AHGuess;

/**
 * Assemble the global system matrix and residuals from the interpolated surface data.
 */
void AssembleGlobalMatrix(const InterpolatedSurface& interpolatedSurface,
                          const AHGuess& ahGuess,
                          std::vector<std::vector<double>>& matrixCoefficients,
                          std::vector<double>& residuals);

}  // namespace exahype2::ccz4::ahf

#endif  // APPARENT_HORIZON_FINDER_ASSEMBLE_MATRIX_H_
