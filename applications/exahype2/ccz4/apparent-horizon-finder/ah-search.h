#ifndef APPARENT_HORIZON_FINDER_AH_SEARCH_H_
#define APPARENT_HORIZON_FINDER_AH_SEARCH_H_

#include "ah-finder-result.h"
#include "interpolation.h"
#include "spectral-solver.h"

#include <cstddef>
#include <limits>
#include <optional>
#include <string>
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

/** Physical role assigned by the user to one independently searched surface. */
enum class AHHorizonKind {
  Individual,
  Common
};

/**
 * One physical horizon and its ordered alternative guesses. The target ID is
 * the stable association between attempts and the intended black hole/common
 * horizon; no geometric post-hoc matching is performed.
 */
struct AHHorizonTarget {
  std::string id;
  AHHorizonKind kind = AHHorizonKind::Individual;
  std::vector<AHSearchAttempt> attempts;
};

/** Shared numerical configuration for any number of physical targets. */
struct AHMultiHorizonOptions {
  std::vector<AHHorizonTarget> targets;
  int lMax = 2;
  CartesianToSurfaceOptions interpolation{};
  SpectralNewtonOptions solver{};
};

/** Ordered attempts and selected result for one physical target. */
struct AHHorizonTargetResult {
  std::string id;
  AHHorizonKind kind = AHHorizonKind::Individual;
  bool converged = false;
  std::size_t successfulAttempt = AHSearchResult::NoSuccessfulAttempt;
  std::optional<AHFinderResult> horizon;
  std::vector<AHSearchAttemptResult> attempts;
};

/** Aggregate result without collapsing partial success into one boolean. */
struct AHMultiHorizonResult {
  bool anyHorizonFound = false;
  bool allTargetsFound = false;
  std::vector<AHHorizonTargetResult> targets;
};

void ValidateAHMultiHorizonOptions(const AHMultiHorizonOptions& options);

/**
 * Search every physical target in vector order. The first successful attempt
 * stops retries for that target only; remaining targets are still searched.
 * One spectral basis and one persistent Cartesian interpolator are shared.
 */
AHMultiHorizonResult FindApparentHorizons(
  const GridData& gridData,
  const AHMultiHorizonOptions& options
);

}  // namespace exahype2::ccz4::ahf

#endif  // APPARENT_HORIZON_FINDER_AH_SEARCH_H_
