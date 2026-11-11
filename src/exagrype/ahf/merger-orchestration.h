#ifndef APPARENT_HORIZON_FINDER_MERGER_ORCHESTRATION_H_
#define APPARENT_HORIZON_FINDER_MERGER_ORCHESTRATION_H_

#include "ah-search.h"

#include <array>
#include <string>
#include <vector>

namespace exahype2::ccz4::ahf {

/** One known member of a binary system and its ordered individual-AH retries. */
struct AHBinaryComponent {
  AHHorizonTarget target;
  double mass = 1.0;
};

/** Opt-in binary-merger planning plus the shared numerical search options. */
struct AHBinaryMergerOptions {
  std::array<AHBinaryComponent, 2> components;
  std::string commonHorizonId = "common";
  /** Multipliers applied to the smallest enclosing common coordinate sphere. */
  std::vector<double> commonRadiusFactors{1.0, 1.5, 2.0};
  int lMax = 2;
  CartesianToSurfaceOptions interpolation{};
  SpectralNewtonOptions solver{};
};

/** Documented trigger/preconditioning diagnostics and the generated plan. */
struct AHBinaryMergerPlan {
  bool commonSearchEnabled = false;
  double separation = 0.0;
  std::array<double, 2> facingRadii{0.0, 0.0};
  ReferencePoint commonReferencePoint{0.0, 0.0, 0.0};
  double commonMass = 0.0;
  double enclosingInitialRadius = 0.0;
  AHMultiHorizonOptions search;
};

struct AHBinaryMergerResult {
  AHBinaryMergerPlan plan;
  AHMultiHorizonResult horizons;
};

void ValidateAHBinaryMergerOptions(const AHBinaryMergerOptions& options);

/** Build the deterministic common-first search plan without accessing grid data. */
AHBinaryMergerPlan BuildBinaryMergerPlan(
  const AHBinaryMergerOptions& options
);

/** Project-consistent name for the documented binary common-horizon check. */
AHBinaryMergerPlan CheckCommonHorizonOfBBH(
  const AHBinaryMergerOptions& options
);

/** Build the plan and execute it with the ordinary multi-horizon machinery. */
AHBinaryMergerResult FindBinaryMergerHorizons(
  const GridData& gridData,
  const AHBinaryMergerOptions& options
);

}  // namespace exahype2::ccz4::ahf

#endif
