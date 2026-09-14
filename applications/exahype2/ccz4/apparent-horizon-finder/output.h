#ifndef APPARENT_HORIZON_FINDER_OUTPUT_H_
#define APPARENT_HORIZON_FINDER_OUTPUT_H_

#include "ah-finder-result.h"
#include "ah-finder-options.h"
#include "ah-search.h"

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

/** Filenames and visualisation resolution for the single-horizon façade. */
struct AHOutputOptions {
  std::string summaryFilename = "apparent-horizon.txt";
  std::string surfaceFilename = "apparent-horizon.vtp";
  AHVisualizationOptions visualization{};
};

/** Validate output filenames and visualisation resolution before solving. */
void ValidateAHOutputOptions(const AHOutputOptions& options);

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
 * Complete single-horizon user workflow.
 *
 * The grid is referenced by one persistent CartesianPointInterpolator for the
 * solve and final reconstruction. A summary is always written after a solver
 * result is obtained. The normal VTP horizon is written only for a converged,
 * reconstructable surface.
 */
AHFinderResult WriteAHToFile(
  const GridData& gridData,
  const AHFinderOptions& finderOptions,
  const AHOutputOptions& outputOptions = {}
);

/**
 * Write ordered diagnostics for every executed attempt and, on success, the
 * selected horizon diagnostics.
 */
void WriteAHSearchSummary(
  const AHSearchResult& result,
  const std::string& filename
);

/**
 * Multi-attempt user workflow. Attempts are paired entries tried in vector
 * order; the first converged, reconstructable horizon stops the search.
 */
AHSearchResult WriteAHToFile(
  const GridData& gridData,
  const AHSearchOptions& searchOptions,
  const AHOutputOptions& outputOptions = {}
);

/**
 * Convert the ASCII apparent-horizon hypersurface output into a VTK file.
 */
void WriteAHToVTK(const std::string& inputFile, const std::string& outputFile);

}  // namespace exahype2::ccz4::ahf

#endif  // APPARENT_HORIZON_FINDER_OUTPUT_H_
