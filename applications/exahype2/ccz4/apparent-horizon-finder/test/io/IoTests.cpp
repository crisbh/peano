#include <iostream>

#include "../../file-input.h"
#include "../../output.h"
#include "../../peano-input.h"

namespace exahype2::ccz4::ahf::tests::io {

void testReadGridFromFile() {
  std::cout << "Placeholder test for file-input ReadGridFromFile()" << std::endl;
}

void testReadGridFromPeano() {
  std::cout << "Placeholder test for peano-input ReadGridFromPeano()" << std::endl;
}

void testWriteAHToFile() {
  std::cout << "Placeholder test for output WriteAHToFile()" << std::endl;
}

void testWriteAHToVTK() {
  std::cout << "Placeholder test for output WriteAHToVTK()" << std::endl;
}

}  // namespace exahype2::ccz4::ahf::tests::io

int main() {
  using namespace exahype2::ccz4::ahf::tests::io;

  std::cout << "Apparent Horizon Finder I/O tests begin..." << std::endl;

  testReadGridFromFile();
  testReadGridFromPeano();
  testWriteAHToFile();
  testWriteAHToVTK();

  std::cout << "All I/O tests completed." << std::endl;
  return 0;
}
