#include "ah-search.h"

#include "ah-finder-options.h"
#include "harmonics.h"
#include "iteration.h"

#include <stdexcept>
#include <unordered_set>
#include <utility>

namespace exahype2::ccz4::ahf {
namespace {

AHFinderOptions SingleAttemptOptions(
  const AHSearchOptions& options,
  const AHSearchAttempt& attempt
) {
  return {
    .referencePoint = attempt.referencePoint,
    .initialRadius = attempt.initialRadius,
    .lMax = options.lMax,
    .interpolation = options.interpolation,
    .solver = options.solver
  };
}

AHSearchResult ExecuteAttempts(
  const SpectralBasis& basis,
  const CartesianPointInterpolator& interpolator,
  const std::vector<AHSearchAttempt>& attempts,
  const SpectralNewtonOptions& solverOptions
) {
  AHSearchResult search;
  search.attempts.reserve(attempts.size());
  for (std::size_t index = 0; index < attempts.size(); ++index) {
    const AHSearchAttempt& attempt = attempts[index];
    const SpectralNewtonResult solver = SolvePDEforAH(
      basis, interpolator, attempt.referencePoint, attempt.initialRadius,
      solverOptions);
    AHFinderResult result = BuildAHFinderResult(
      basis, interpolator, attempt.referencePoint, solver);
    const bool succeeded = result.solver.converged && result.surfaceAvailable;
    search.attempts.push_back({
      .referencePoint = attempt.referencePoint,
      .initialRadius = attempt.initialRadius,
      .result = result,
      .succeeded = succeeded
    });
    if (succeeded) {
      search.converged = true;
      search.successfulAttempt = index;
      search.horizon = std::move(result);
      return search;
    }
  }
  return search;
}

}  // namespace

void ValidateAHSearchOptions(const AHSearchOptions& options) {
  if (options.attempts.empty()) {
    throw std::invalid_argument(
      "ValidateAHSearchOptions: at least one search attempt is required");
  }
  for (const AHSearchAttempt& attempt : options.attempts) {
    ValidateAHFinderOptions(SingleAttemptOptions(options, attempt));
  }
}

AHSearchResult FindApparentHorizon(const GridData& gridData,
                                   const AHSearchOptions& options) {
  ValidateAHSearchOptions(options);

  const SpectralBasis basis = BuildSpectralBasis(options.lMax);
  const CartesianPointInterpolator interpolator(gridData, options.interpolation);

  return ExecuteAttempts(basis, interpolator, options.attempts, options.solver);
}

void ValidateAHMultiHorizonOptions(const AHMultiHorizonOptions& options) {
  if (options.targets.empty()) {
    throw std::invalid_argument(
      "ValidateAHMultiHorizonOptions: at least one horizon target is required");
  }
  std::unordered_set<std::string> ids;
  for (const AHHorizonTarget& target : options.targets) {
    if (target.id.empty()) {
      throw std::invalid_argument(
        "ValidateAHMultiHorizonOptions: target IDs must not be empty");
    }
    if (!ids.insert(target.id).second) {
      throw std::invalid_argument(
        "ValidateAHMultiHorizonOptions: target IDs must be unique");
    }
    ValidateAHSearchOptions({
      .attempts = target.attempts,
      .lMax = options.lMax,
      .interpolation = options.interpolation,
      .solver = options.solver
    });
  }
}

AHMultiHorizonResult FindApparentHorizons(
  const GridData& gridData,
  const AHMultiHorizonOptions& options
) {
  ValidateAHMultiHorizonOptions(options);
  const SpectralBasis basis = BuildSpectralBasis(options.lMax);
  const CartesianPointInterpolator interpolator(gridData, options.interpolation);

  AHMultiHorizonResult result;
  result.allTargetsFound = true;
  result.targets.reserve(options.targets.size());
  for (const AHHorizonTarget& target : options.targets) {
    AHSearchResult search = ExecuteAttempts(
      basis, interpolator, target.attempts, options.solver);
    result.anyHorizonFound = result.anyHorizonFound || search.converged;
    result.allTargetsFound = result.allTargetsFound && search.converged;
    result.targets.push_back({
      .id = target.id,
      .kind = target.kind,
      .converged = search.converged,
      .successfulAttempt = search.successfulAttempt,
      .horizon = std::move(search.horizon),
      .attempts = std::move(search.attempts)
    });
  }
  return result;
}

}  // namespace exahype2::ccz4::ahf
