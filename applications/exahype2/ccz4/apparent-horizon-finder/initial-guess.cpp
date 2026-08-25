#include "initial-guess.h"

#include <cmath>
#include <stdexcept>

namespace exahype2::ccz4::ahf {

std::vector<double> MakeSphericalInitialCoefficients(
  const SpectralBasis& basis,
  double radius
) {
  if (!std::isfinite(radius) || radius <= 0.0) {
    throw std::invalid_argument(
      "MakeSphericalInitialCoefficients: radius must be finite and strictly positive"
    );
  }

  const std::vector<double> field(basis.nPoints, radius);
  std::vector<double> coefficients = Project(basis, field);
  if (coefficients.size() != basis.nModes) {
    throw std::runtime_error(
      "MakeSphericalInitialCoefficients: projected coefficient count does not match basis.nModes"
    );
  }
  return coefficients;
}

}  // namespace exahype2::ccz4::ahf
