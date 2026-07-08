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
  testGridDataStruct();

  std::cout << "All data tests completed." << std::endl;
  return 0;
}
