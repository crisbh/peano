#ifndef APPARENT_HORIZON_FINDER_COMMON_HORIZON_CANDIDATES_H_
#define APPARENT_HORIZON_FINDER_COMMON_HORIZON_CANDIDATES_H_

#include "ah-finder-result.h"
#include "ah-search.h"

#include <optional>
#include <string>
#include <vector>

namespace exahype2::ccz4::ahf {

/** User identity and retry configuration for one physical black hole. */
struct AHPhysicalBlackHole {
  std::string id;
  double mass = 1.0;
  std::vector<AHSearchAttempt> attempts;
  /** Physical puncture/location estimate, distinct from numerical AH origins. */
  std::optional<ReferencePoint> physicalPosition;
};

/** A successfully reconstructed individual horizon associated with its ID. */
struct AHIndividualHorizonObservation {
  std::string id;
  double mass = 1.0;
  AHFinderResult horizon;
};

/** Directional separation diagnostic for one unordered pair. */
struct AHHorizonPairGap {
  std::string firstId;
  std::string secondId;
  double centreDistance = 0.0;
  double firstFacingRadius = 0.0;
  double secondFacingRadius = 0.0;
  double surfaceGap = 0.0;
  bool proximityTriggered = false;
};

struct AHCommonHorizonCandidateOptions {
  /** Trigger when surfaceGap <= searchMargin; coordinate-length units. */
  double searchMargin = 0.0;
  bool includeTriggeredPairs = true;
  bool includeConnectedComponents = true;
  /** Always add the group containing every observed black hole. */
  bool unconditionalAllBlackHoles = false;
  /** Always add these ID subsets, independently of proximity. */
  std::vector<std::vector<std::string>> requestedSubsets;
};

/** Why a canonical group survived deterministic deduplication. */
struct AHCommonHorizonCandidateReason {
  bool triggeredPair = false;
  bool connectedComponent = false;
  bool explicitlyRequested = false;
  bool unconditionalAllBlackHoles = false;
};

/** Geometry and initial spherical guess for one prospective common AH solve. */
struct AHCommonHorizonCandidate {
  /** Sorted stable IDs; also the canonical deduplication key. */
  std::vector<std::string> memberIds;
  AHCommonHorizonCandidateReason reason;
  ReferencePoint referencePoint{0.0, 0.0, 0.0};
  double totalMass = 0.0;
  double enclosingInitialRadius = 0.0;
};

struct AHCommonHorizonCandidatePlan {
  std::vector<AHHorizonPairGap> pairGaps;
  std::vector<AHCommonHorizonCandidate> candidates;
};

void ValidateAHPhysicalBlackHoles(
  const std::vector<AHPhysicalBlackHole>& blackHoles
);

/** Convert user identities into the existing individual-target representation. */
std::vector<AHHorizonTarget> MakeIndividualHorizonTargets(
  const std::vector<AHPhysicalBlackHole>& blackHoles
);

void ValidateAHCommonHorizonCandidateOptions(
  const AHCommonHorizonCandidateOptions& options
);

/**
 * Plan prospective common searches from already successful individual AHs.
 * Proximity only creates a candidate; it is not evidence that a common AH
 * exists, and this function deliberately performs no Newton solve.
 */
AHCommonHorizonCandidatePlan PlanCommonHorizonCandidates(
  const std::vector<AHIndividualHorizonObservation>& observations,
  const AHCommonHorizonCandidateOptions& options = {}
);

}  // namespace exahype2::ccz4::ahf

#endif
