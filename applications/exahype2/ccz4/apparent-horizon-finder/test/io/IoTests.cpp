#include <iostream>

#include "../..//apparent-horizon-finder.h"

namespace exahype2::ccz4::ahf::tests::io {

void testReadGridFromFile() {
  std::cout << "Placeholder test for ReadGridFromFile()" << std::endl;
}

void testReadGridFromPeano() {
  std::cout << "Placeholder test for ReadGridFromPeano()" << std::endl;
}

void testWriteAHToFile() {
  std::cout << "Placeholder test for WriteAHToFile()" << std::endl;
}

void testWriteAHToVTK() {
  std::cout << "Placeholder test for WriteAHToVTK()" << std::endl;
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
