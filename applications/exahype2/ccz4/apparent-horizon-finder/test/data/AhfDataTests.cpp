#include <cassert>
#include <iostream>

#include "../../ahf-data.h"

namespace exahype2::ccz4::ahf::tests::data {

void testReferencePointStruct() {
  const ReferencePoint referencePoint{1.0, 0.25, 0.5};

  assert(referencePoint.x == 1.0);
  assert(referencePoint.y == 0.25);
  assert(referencePoint.z == 0.5);

  std::cout << "ReferencePoint stores x = " << referencePoint.x << ", y = " << referencePoint.y << ", and z = " << referencePoint.z << std::endl;
}

void testAHGuessStruct() {
  const std::vector<AHGuess> initialAHGuesses = {
    {{{1.0, 0.0, 0.0}, {1.1, 0.5, 1.0}}},
    {{{2.0, 0.0, 0.0}, {2.2, 0.5, 1.0}}}
  };

  assert(initialAHGuesses.size() == 2);
  assert(initialAHGuesses[0].points.size() == 2);
  assert(initialAHGuesses[0].points[1].r == 1.1);
  assert(initialAHGuesses[0].points[1].rho == 0.5);
  assert(initialAHGuesses[0].points[1].sigma == 1.0);
  assert(initialAHGuesses[1].points[0].r == 2.0);

  std::cout << "AHGuess stores lists of (r, rho, sigma) points, and initialAHGuesses stores multiple guesses" << std::endl;
}

void testGridDataStruct() {
  GridData gridData;
  gridData.coordinates = {
    {0.0, 0.1, 0.2},
    {1.0, 1.1, 1.2}
  };
  gridData.setNumberOfVariables(3);

  gridData.value(0, 0) = 10.0;
  gridData.value(0, 1) = 11.0;
  gridData.value(0, 2) = 12.0;
  gridData.value(1, 0) = 20.0;
  gridData.value(1, 1) = 21.0;
  gridData.value(1, 2) = 22.0;

  assert(gridData.numberOfPoints() == 2);
  assert(gridData.numberOfVariables == 3);
  assert(gridData.isConsistent());
  assert(gridData.coordinates[1].x == 1.0);
  assert(gridData.coordinates[1].y == 1.1);
  assert(gridData.coordinates[1].z == 1.2);
  assert(gridData.value(0, 2) == 12.0);
  assert(gridData.value(1, 1) == 21.0);

  const GridData& constGridData = gridData;
  assert(constGridData.value(1, 2) == 22.0);

  std::cout << "GridData stores coordinates and point-major scalar values" << std::endl;
}

}  // namespace exahype2::ccz4::ahf::tests::data

int main() {
  using namespace exahype2::ccz4::ahf::tests::data;

  std::cout << "Apparent Horizon Finder data tests begin..." << std::endl;

  testReferencePointStruct();
  testAHGuessStruct();
  testGridDataStruct();

  std::cout << "All data tests completed." << std::endl;
  return 0;
}
