#ifndef APPARENT_HORIZON_FINDER_AH_SEARCH_H_
#define APPARENT_HORIZON_FINDER_AH_SEARCH_H_

#include "ah-finder-result.h"
#include "interpolation.h"
#include "spectral-solver.h"

#include <cstddef>
#include <limits>
#include <optional>
#include <vector>

namespace exahype2::ccz4::ahf {

/** One explicit, ordered spherical initial-guess/reference-point pair. */
struct AHSearchAttempt {
  ReferencePoint referencePoint{0.0, 0.0, 0.0};
  double initialRadius = 1.0;
};

/** Shared numerical configuration and an arbitrary ordered attempt list. */
struct AHSearchOptions {
  std::vector<AHSearchAttempt> attempts;
  int lMax = 2;
  CartesianToSurfaceOptions interpolation{};
  SpectralNewtonOptions solver{};
};

/** Complete diagnostics from one executed attempt. */
struct AHSearchAttemptResult {
  ReferencePoint referencePoint{0.0, 0.0, 0.0};
  double initialRadius = 1.0;
  AHFinderResult result;
  /** True only when Newton converged and final reconstruction succeeded. */
  bool succeeded = false;
};

/** Aggregate ordered-search outcome. */
struct AHSearchResult {
  static constexpr std::size_t NoSuccessfulAttempt =
    std::numeric_limits<std::size_t>::max();

  bool converged = false;
  std::size_t successfulAttempt = NoSuccessfulAttempt;
  std::optional<AHFinderResult> horizon;
  std::vector<AHSearchAttemptResult> attempts;
};

/** Reject an empty list or any structurally invalid shared/attempt option. */
void ValidateAHSearchOptions(const AHSearchOptions& options);

/**
 * Try explicit attempt pairs in vector order and stop at the first converged,
 * reconstructable horizon. One persistent Cartesian interpolator is reused.
 */
AHSearchResult FindApparentHorizon(
  const GridData& gridData,
  const AHSearchOptions& options
);

}  // namespace exahype2::ccz4::ahf

#endif  // APPARENT_HORIZON_FINDER_AH_SEARCH_H_
