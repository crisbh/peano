#include <iostream>

#include "../..//apparent-horizon-finder.h"

namespace exahype2::ccz4::ahf::tests::iteration {

void testSolvePDEforAH() {
  std::cout << "Placeholder test for SolvePDEforAH" << std::endl;
}

void testCheckCommonHorizonOfBBH() {
  std::cout << "Placeholder test for CheckCommonHorizonOfBBH" << std::endl;
}

}  // namespace exahype2::ccz4::ahf::tests::iteration

int main() {
  using namespace exahype2::ccz4::ahf::tests::iteration;

  std::cout << "Apparent Horizon Finder iteration tests begin..." << std::endl;

  testSolvePDEforAH();
  testCheckCommonHorizonOfBBH();

  std::cout << "All iteration tests completed." << std::endl;
  return 0;
}
