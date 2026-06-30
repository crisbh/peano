#ifndef APPARENT_HORIZON_FINDER_FILE_INPUT_H_
#define APPARENT_HORIZON_FINDER_FILE_INPUT_H_

#include <string>
#include <vector>

namespace exahype2::ccz4::ahf {

/**
 * Load grid dataset from VTK/HDF5 files.
 */
void ReadGridFromFile(const std::string& inputFile);

}  // namespace exahype2::ccz4::ahf

#endif  // APPARENT_HORIZON_FINDER_FILE_INPUT_H_
