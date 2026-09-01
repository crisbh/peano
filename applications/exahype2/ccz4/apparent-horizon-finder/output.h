#ifndef APPARENT_HORIZON_FINDER_OUTPUT_H_
#define APPARENT_HORIZON_FINDER_OUTPUT_H_

#include "ah-finder-result.h"

#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace exahype2::ccz4::ahf {

struct AHGuess;

/** Resolution of the output mesh only; it does not affect the AH solve. */
struct AHVisualizationOptions {
  int numberOfRhoPoints = 64;
  int numberOfSigmaPoints = 33;
};

/** Dense closed surface reconstructed directly from the spectral coefficients. */
struct AHVisualizationSurface {
  std::vector<std::array<double, 3>> points;
  std::vector<double> radius;
  std::vector<std::array<std::size_t, 3>> triangles;
};

/** Write all solver and reconstructed-surface diagnostics as deterministic text. */
void WriteAHSummary(const AHFinderResult& result, const std::string& filename);

/**
 * Evaluate the retained spectral representation on an independent dense grid.
 * The grid has one vertex at each pole and periodic azimuthal rings between.
 */
AHVisualizationSurface BuildAHVisualizationSurface(
  const AHFinderResult& result,
  const AHVisualizationOptions& options = {}
);

/** Write a converged dense horizon surface as VTK XML PolyData (.vtp). */
void WriteAHSurfaceToVTP(
  const AHFinderResult& result,
  const std::string& filename,
  const AHVisualizationOptions& options = {}
);

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

#endif  // APPARENT_HORIZON_FINDER_OUTPUT_H_
