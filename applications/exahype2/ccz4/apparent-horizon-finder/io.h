#ifndef APPARENT_HORIZON_FINDER_IO_H_
#define APPARENT_HORIZON_FINDER_IO_H_

#include <string>
#include <vector>

namespace exahype2::ccz4::ahf {

struct PeanoGridData;
struct AHGuess;
struct ReferencePoint;

/**
 * Load grid dataset from VTK/HDF5 files.
 */
void ReadGridFromFile(const std::string& inputFile);

/**
 * Load grid dataset from the Peano simulation on-the-fly.
 */
void ReadGridFromPeano(const PeanoGridData& peanoGridData);

/**
 * Write the apparent-horizon hypersurface to an ASCII file.
 */
void WriteAHToFile(const std::vector<AHGuess>& initialAHGuesses,
                   const std::vector<ReferencePoint>& referencePoints,
                   const std::vector<double>& masses,
                   const std::string& outputFile);

/**
 * Convert the ASCII apparent-horizon hypersurface output into a VTK file.
 */
void WriteAHToVTK(const std::string& inputFile, const std::string& outputFile);

}  // namespace exahype2::ccz4::ahf

#endif  // APPARENT_HORIZON_FINDER_IO_H_
