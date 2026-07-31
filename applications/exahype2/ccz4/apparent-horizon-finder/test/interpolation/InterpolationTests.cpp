#include <iostream>

#include "../..//apparent-horizon-finder.h"

namespace exahype2::ccz4::ahf::tests::interpolation {

static_assert(requires(const GridData& gridData,
                       const AHGuess& ahGuess,
                       const ReferencePoint& referencePoint,
                       InterpolatedShell& interpolatedShell,
                       const std::vector<AngularSample>& angularSamples) {
  InterpolateFromCartesianToShell(
    gridData,
    ahGuess,
    referencePoint,
    interpolatedShell
  );
  InterpolateFromCartesianToShell(
    gridData,
    ahGuess,
    referencePoint,
    interpolatedShell,
    angularSamples
  );
});

void testInterpolateFromCartesianToShell() {
  std::cout << "Placeholder test for InterpolateFromCartesianToShell" << std::endl;
}

void testInterpolateFromShellToSurface() {
  std::cout << "Placeholder test for InterpolateFromShellToSurface" << std::endl;
}

}  // namespace exahype2::ccz4::ahf::tests::interpolation

int main() {
  using namespace exahype2::ccz4::ahf::tests::interpolation;

  std::cout << "Apparent Horizon Finder interpolation tests begin..." << std::endl;

  testInterpolateFromCartesianToShell();
  testInterpolateFromShellToSurface();

  std::cout << "All interpolation tests completed." << std::endl;
  return 0;
}
