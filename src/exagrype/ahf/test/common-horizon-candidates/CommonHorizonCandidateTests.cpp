#include "../../apparent-horizon-finder.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace exahype2::ccz4::ahf::tests {
namespace {

AHIndividualHorizonObservation Sphere(
  std::string id, double mass, ReferencePoint centre, double radius
) {
  const SpectralBasis basis = BuildSpectralBasis(2);
  AHFinderResult horizon;
  horizon.solver.converged = true;
  horizon.solver.coefficients = MakeSphericalInitialCoefficients(basis, radius);
  horizon.referencePoint = centre;
  horizon.lMax = 2;
  horizon.minRadius = radius;
  horizon.meanRadius = radius;
  horizon.maxRadius = radius;
  horizon.surfaceAvailable = true;
  return {std::move(id), mass, std::move(horizon)};
}

const AHCommonHorizonCandidate& FindCandidate(
  const AHCommonHorizonCandidatePlan& plan,
  const std::vector<std::string>& members
) {
  for (const auto& candidate : plan.candidates) {
    if (candidate.memberIds == members) return candidate;
  }
  throw std::runtime_error("candidate not found");
}

void TestPhysicalConfigurationConversion() {
  const std::vector<AHPhysicalBlackHole> blackHoles{
    {.id = "A", .mass = 1.0,
     .attempts = {{{0.0, 0.0, 0.0}, 1.0}}},
    {.id = "B", .mass = 2.0,
     .attempts = {{{2.0, 0.0, 0.0}, 0.5}}}
  };
  const auto targets = MakeIndividualHorizonTargets(blackHoles);
  assert(targets.size() == 2);
  assert(targets[0].id == "A" && targets[1].id == "B");
  assert(targets[0].kind == AHHorizonKind::Individual);
  assert(*targets[1].massScale == 2.0);

  bool threw = false;
  auto invalid = blackHoles;
  invalid[1].id = "A";
  try { (void)MakeIndividualHorizonTargets(invalid); }
  catch (const std::invalid_argument&) { threw = true; }
  assert(threw);
}

void TestMarginBeforeContact() {
  const std::vector<AHIndividualHorizonObservation> observations{
    Sphere("A", 1.0, {0.0, 0.0, 0.0}, 1.0),
    Sphere("B", 1.0, {2.2, 0.0, 0.0}, 1.0)
  };
  const auto withoutMargin = PlanCommonHorizonCandidates(observations);
  assert(withoutMargin.pairGaps.size() == 1);
  assert(std::abs(withoutMargin.pairGaps[0].surfaceGap - 0.2) < 1.0e-12);
  assert(!withoutMargin.pairGaps[0].proximityTriggered);
  assert(withoutMargin.candidates.empty());

  const auto withMargin = PlanCommonHorizonCandidates(
    observations, {.searchMargin = 0.25});
  assert(withMargin.pairGaps[0].proximityTriggered);
  assert(withMargin.candidates.size() == 1); // pair == component, deduplicated
  assert(withMargin.candidates[0].memberIds ==
         std::vector<std::string>({"A", "B"}));
  assert(withMargin.candidates[0].reason.triggeredPair);
  assert(withMargin.candidates[0].reason.connectedComponent);
}

void TestSingleObservationNeverCreatesCommonCandidate() {
  const auto plan = PlanCommonHorizonCandidates(
    {Sphere("single", 1.0, {0.0, 0.0, 0.0}, 1.0)},
    {.searchMargin = 100.0,
     .includeTriggeredPairs = true,
     .includeConnectedComponents = true,
     .unconditionalAllBlackHoles = true});
  assert(plan.pairGaps.empty());
  assert(plan.candidates.empty());
}

void TestPairsComponentsExplicitAndUnconditional() {
  const std::vector<AHIndividualHorizonObservation> observations{
    Sphere("C", 3.0, {3.0, 0.0, 0.0}, 1.0),
    Sphere("A", 1.0, {0.0, 0.0, 0.0}, 1.0),
    Sphere("B", 2.0, {1.5, 0.0, 0.0}, 1.0),
    Sphere("D", 4.0, {10.0, 0.0, 0.0}, 0.5)
  };
  const auto plan = PlanCommonHorizonCandidates(
    observations,
    {.searchMargin = 0.0,
     .includeTriggeredPairs = true,
     .includeConnectedComponents = true,
     .unconditionalAllBlackHoles = true,
     .requestedSubsets = {{"D", "A"}, {"B", "A"}}}
  );

  // Triggered chain A-B-C gives AB, BC and component ABC. Explicit AD and
  // unconditional ABCD are added. Explicit BA merges with triggered AB.
  assert(plan.candidates.size() == 5);
  const auto& ab = FindCandidate(plan, {"A", "B"});
  assert(ab.reason.triggeredPair && ab.reason.explicitlyRequested);
  const auto& abc = FindCandidate(plan, {"A", "B", "C"});
  assert(abc.reason.connectedComponent);
  assert(std::abs(abc.totalMass - 6.0) < 1.0e-14);
  assert(std::abs(abc.referencePoint.x - 2.0) < 1.0e-14);
  assert(std::abs(abc.enclosingInitialRadius - 3.0) < 1.0e-14);
  const auto& all = FindCandidate(plan, {"A", "B", "C", "D"});
  assert(all.reason.unconditionalAllBlackHoles);

  for (std::size_t i = 1; i < plan.candidates.size(); ++i) {
    assert(plan.candidates[i - 1].memberIds < plan.candidates[i].memberIds);
  }
}

void TestDirectionalDistortedSurface() {
  const SpectralBasis basis = BuildSpectralBasis(2);
  AHIndividualHorizonObservation first =
    Sphere("A", 1.0, {0.0, 0.0, 0.0}, 1.0);
  std::vector<double> radii(basis.nPoints);
  for (std::size_t p = 0; p < radii.size(); ++p) {
    radii[p] = 1.0 + 0.2 * std::cos(basis.rho[p]) *
      std::cos(basis.sigma[p]);
  }
  first.horizon.solver.coefficients =
    MakeDistortedInitialCoefficients(basis, radii);
  first.horizon.maxRadius = *std::max_element(radii.begin(), radii.end());
  const auto plan = PlanCommonHorizonCandidates({
    first, Sphere("B", 1.0, {2.1, 0.0, 0.0}, 1.0)
  });
  assert(plan.pairGaps[0].firstFacingRadius > 1.19);
  assert(plan.pairGaps[0].firstFacingRadius < 1.21);
  assert(plan.pairGaps[0].proximityTriggered);
}

void TestInvalidObservationsAndRequests() {
  auto observations = std::vector<AHIndividualHorizonObservation>{
    Sphere("A", 1.0, {0.0, 0.0, 0.0}, 1.0),
    Sphere("B", 1.0, {3.0, 0.0, 0.0}, 1.0)
  };
  observations[0].horizon.solver.converged = false;
  bool threw = false;
  try { (void)PlanCommonHorizonCandidates(observations); }
  catch (const std::invalid_argument&) { threw = true; }
  assert(threw);

  observations[0].horizon.solver.converged = true;
  threw = false;
  try {
    (void)PlanCommonHorizonCandidates(
      observations, {.requestedSubsets = {{"A", "unknown"}}});
  } catch (const std::invalid_argument&) { threw = true; }
  assert(threw);
}

}  // namespace
}  // namespace exahype2::ccz4::ahf::tests

int main() {
  using namespace exahype2::ccz4::ahf::tests;
  TestPhysicalConfigurationConversion();
  TestMarginBeforeContact();
  TestSingleObservationNeverCreatesCommonCandidate();
  TestPairsComponentsExplicitAndUnconditional();
  TestDirectionalDistortedSurface();
  TestInvalidObservationsAndRequests();
}
