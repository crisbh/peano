#include <iostream>

#include "../..//apparent-horizon-finder.h"

namespace exahype2::ccz4::ahf::tests::assemble_matrix {

void testAssembleGlobalMatrix() {
  std::cout << "Placeholder test for AssembleGlobalMatrix" << std::endl;
}

}  // namespace exahype2::ccz4::ahf::tests::assemble_matrix

int main() {
  using namespace exahype2::ccz4::ahf::tests::assemble_matrix;

  std::cout << "Apparent Horizon Finder assemble-matrix tests begin..." << std::endl;

  testAssembleGlobalMatrix();

  std::cout << "All assemble-matrix tests completed." << std::endl;
  return 0;
}
