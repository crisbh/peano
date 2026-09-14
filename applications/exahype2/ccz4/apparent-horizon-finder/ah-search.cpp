#include "ah-search.h"

#include "ah-finder-options.h"
#include "harmonics.h"
#include "iteration.h"

#include <stdexcept>
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

  AHSearchResult search;
  search.attempts.reserve(options.attempts.size());
  for (std::size_t index = 0; index < options.attempts.size(); ++index) {
    const AHSearchAttempt& attempt = options.attempts[index];
    const SpectralNewtonResult solver = SolvePDEforAH(
      basis,
      interpolator,
      attempt.referencePoint,
      attempt.initialRadius,
      options.solver
    );
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

}  // namespace exahype2::ccz4::ahf
