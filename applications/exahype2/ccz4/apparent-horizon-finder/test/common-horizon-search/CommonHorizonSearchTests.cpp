#include "../../apparent-horizon-finder.h"

#include <array>
#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>

namespace exahype2::ccz4::ahf::tests {
namespace {

AHFinderResult Sphere(ReferencePoint centre, double radius) {
  const SpectralBasis basis = BuildSpectralBasis(2);
  AHFinderResult result;
  result.solver.converged = true;
  result.solver.coefficients = MakeSphericalInitialCoefficients(basis, radius);
  result.referencePoint = centre;
  result.lMax = 2;
  result.minRadius = result.meanRadius = result.maxRadius = radius;
  result.spectralResidualInfinityNorm = 1.0e-12;
  result.maxAbsExpansion = 2.0e-12;
  result.surfaceAvailable = true;
  return result;
}

AHIndividualHorizonObservation Observation(
  std::string id, double mass, ReferencePoint centre, double radius) {
  return {std::move(id), mass, Sphere(centre, radius)};
}

GridData MakeGrid() {
  GridData grid;
  for (int ix = -2; ix <= 2; ++ix) for (int iy = -2; iy <= 2; ++iy)
    for (int iz = -2; iz <= 2; ++iz) {
      grid.coordinates.push_back({double(ix), double(iy), double(iz)});
      const std::array<double, NumberOfCompactFields> row{
        1,0,0,1,0,1, 1,0,0,1,0,1, 3,1};
      grid.values.insert(grid.values.end(), row.begin(), row.end());
    }
  return grid;
}

CartesianToSurfaceOptions Interpolation() {
  return {.method = CartesianToSurfaceMethod::UnweightedPolynomial,
          .numberOfNeighbours = 64, .polynomialDegree = 1,
          .allowOrderReduction = false};
}

void TestCandidateConversion() {
  AHCommonHorizonCandidate candidate;
  candidate.memberIds = {"C", "A", "B"};
  candidate.reason = {.triggeredPair = true, .explicitlyRequested = true};
  candidate.referencePoint = {1, 2, 3};
  candidate.totalMass = 6;
  candidate.enclosingInitialRadius = 4;
  const auto target = MakeCommonHorizonTarget(
    candidate, {.initialRadiusInflation = 1.25,
                .retryRadiusMultipliers = {1, 1.2, 2}});
  assert(target.id == "common[A+B+C]");
  assert(target.kind == AHHorizonKind::Common);
  assert(*target.massScale == 6);
  assert(target.attempts.size() == 3);
  assert(target.attempts[0].initialRadius == 5);
  assert(target.attempts[1].initialRadius == 6);
  assert(target.attempts[2].initialRadius == 10);
  for (const auto& attempt : target.attempts) {
    assert(attempt.referencePoint.x == 1 && attempt.referencePoint.y == 2);
  }
}

void TestValidationAndEnclosure() {
  const auto a = Observation("A", 1, {-0.5, 0, 0}, 0.4);
  const auto b = Observation("B", 1, { 0.5, 0, 0}, 0.4);
  AHFinderResult common = Sphere({0, 0, 0}, 1.1);
  auto validation = ValidateCommonHorizon(common, {a, b});
  assert(validation.valid);
  assert(validation.numericallyValid);
  assert(validation.status ==
         AHCommonHorizonValidationStatus::FullyValidated);
  assert(validation.memberHorizonEnclosure ==
         AHEnclosureValidationStatus::Passed);
  assert(validation.enclosesMembers);
  assert(validation.minimumEnclosureClearance > 0.19);

  common = Sphere({0, 0, 0}, 0.7);
  validation = ValidateCommonHorizon(common, {a, b});
  assert(!validation.valid && !validation.enclosesMembers);
  assert(validation.status == AHCommonHorizonValidationStatus::Rejected);

  common = Sphere({0, 0, 0}, 1.1);
  common.maxAbsExpansion = 1.0e-2;
  validation = ValidateCommonHorizon(common, {a, b});
  assert(!validation.valid && !validation.pointwiseExpansionValid);

  common = Sphere({0, 0, 0}, 1.1);
  common.solver.coefficients[0] = std::numeric_limits<double>::infinity();
  validation = ValidateCommonHorizon(common, {a, b});
  assert(!validation.valid && !validation.coefficientsFinite);
}

void TestProvisionalAndPhysicalEnclosure() {
  const auto a = Observation("A", 1, {-0.5, 0, 0}, 0.4);
  AHFinderResult common = Sphere({0, 0, 0}, 1.1);
  auto validation = ValidateCommonHorizon(
    common, {a}, {{-0.5, 0, 0}, {0.5, 0, 0}}, 2);
  assert(validation.numericallyValid);
  assert(!validation.valid);
  assert(validation.status ==
         AHCommonHorizonValidationStatus::ProvisionallyValidated);
  assert(validation.memberHorizonEnclosure ==
         AHEnclosureValidationStatus::Unavailable);
  assert(validation.physicalPositionEnclosure ==
         AHEnclosureValidationStatus::Passed);

  validation = ValidateCommonHorizon(
    common, {a}, {{-0.5, 0, 0}, {2.0, 0, 0}}, 2);
  assert(validation.numericallyValid);
  assert(validation.status == AHCommonHorizonValidationStatus::Rejected);
  assert(validation.physicalPositionEnclosure ==
         AHEnclosureValidationStatus::Failed);
}

AHCommonHorizonSearchOptions OneHoleOptions(ReferencePoint referencePoint) {
  AHCommonHorizonSearchOptions options;
  options.blackHoles = {{
    .id = "single", .mass = 1.0,
    .attempts = {{referencePoint, 1.3}},
    .physicalPosition = ReferencePoint{0, 0, 0}
  }};
  options.lMax = 0;
  options.interpolation = Interpolation();
  // N=1 must remain individual-only even when this global option is enabled.
  options.candidates.unconditionalAllBlackHoles = true;
  return options;
}

void TestZeroAndOneBlackHole() {
  AHCommonHorizonSearchOptions empty;
  bool threw = false;
  try { (void)FindApparentHorizons(MakeGrid(), empty); }
  catch (const std::invalid_argument& error) {
    threw = std::string(error.what()).find("at least one black hole") !=
      std::string::npos;
  }
  assert(threw);

  const auto success = FindApparentHorizons(
    MakeGrid(), OneHoleOptions({0, 0, 0}));
  assert(success.individuals.targets.size() == 1);
  assert(success.individuals.targets[0].converged);
  assert(success.candidatePlan.pairGaps.empty());
  assert(success.candidatePlan.candidates.empty());
  assert(success.commonHorizons.empty());

  const auto failure = FindApparentHorizons(
    MakeGrid(), OneHoleOptions({100, 0, 0}));
  assert(failure.individuals.targets.size() == 1);
  assert(!failure.individuals.targets[0].converged);
  assert(!failure.individuals.targets[0].attempts.empty());
  assert(failure.candidatePlan.candidates.empty());
  assert(failure.commonHorizons.empty());

  const auto directory = std::filesystem::temp_directory_path();
  const std::string successPrefix =
    (directory / "ahf-single-success").string();
  const std::string successSummary = successPrefix + "-summary.txt";
  const auto writtenSuccess = WriteAHToFile(
    MakeGrid(), OneHoleOptions({0, 0, 0}),
    {.summaryFilename = successSummary,
     .horizonFilenamePrefix = successPrefix,
     .visualization = {8, 5}});
  assert(writtenSuccess.individuals.targets.size() == 1);
  assert(writtenSuccess.commonHorizons.empty());
  assert(std::filesystem::exists(successSummary));
  assert(std::filesystem::exists(
    successPrefix + "-individual-0-single.txt"));
  assert(std::filesystem::exists(
    successPrefix + "-individual-0-single.vtp"));

  const std::string failurePrefix =
    (directory / "ahf-single-failure").string();
  const std::string failureSummary = failurePrefix + "-summary.txt";
  const auto writtenFailure = WriteAHToFile(
    MakeGrid(), OneHoleOptions({100, 0, 0}),
    {.summaryFilename = failureSummary,
     .horizonFilenamePrefix = failurePrefix,
     .visualization = {8, 5}});
  assert(writtenFailure.individuals.targets.size() == 1);
  assert(!writtenFailure.individuals.targets[0].converged);
  assert(writtenFailure.commonHorizons.empty());
  assert(std::filesystem::exists(failureSummary));
  assert(std::filesystem::exists(
    failurePrefix + "-individual-0-single.txt"));
  assert(!std::filesystem::exists(
    failurePrefix + "-individual-0-single.vtp"));

  for (const auto& path : {
         successSummary,
         successPrefix + "-individual-0-single.txt",
         successPrefix + "-individual-0-single.vtp",
         failureSummary,
         failurePrefix + "-individual-0-single.txt"}) {
    std::filesystem::remove(path);
  }
}

void TestExplicitFallbackAfterIndividualFailure() {
  AHCommonHorizonSearchOptions options;
  options.blackHoles = {
    {.id = "failed", .mass = 1,
     .attempts = {{{100, 0, 0}, 1.0}},
     .physicalPosition = ReferencePoint{-0.1, 0, 0}},
    {.id = "found", .mass = 2,
     .attempts = {{{0.1, 0, 0}, 1.3}},
     .physicalPosition = ReferencePoint{0.1, 0, 0}}
  };
  options.candidates.requestedSubsets = {{"found", "failed"}};
  options.commonTargets.retryRadiusMultipliers = {1};
  options.lMax = 0;
  options.interpolation = Interpolation();
  options.validation.spectralResidualTolerance = 1.0;
  options.validation.pointwiseExpansionTolerance = 1.0;
  options.validation.enclosureTolerance = 0.5;
  const auto result = FindCommonApparentHorizons(MakeGrid(), options);
  assert(result.individuals.targets.size() == 2);
  assert(!result.individuals.targets[0].converged);
  assert(result.individuals.targets[1].converged);
  assert(result.commonHorizons.size() == 1);
  assert(result.commonHorizons[0].id == "common[failed+found]");
  assert(result.commonHorizons[0].reason.explicitlyRequested);
  assert(result.commonHorizons[0].search.converged);
  assert(result.commonHorizons[0].validation.numericallyValid);
  assert(result.commonHorizons[0].validation.status ==
         AHCommonHorizonValidationStatus::ProvisionallyValidated);
  assert(!result.commonHorizons[0].validation.valid);
  assert(result.commonHorizons[0].validation.physicalPositionEnclosure ==
         AHEnclosureValidationStatus::Passed);
  assert(!result.commonHorizons[0].search.attempts.empty());
  assert(std::abs(result.candidatePlan.candidates[0].referencePoint.x -
                  (0.1 / 3.0)) < 1.0e-14);

  const auto directory = std::filesystem::temp_directory_path();
  const std::string prefix = (directory / "ahf-provisional").string();
  const std::string summary = prefix + "-summary.txt";
  (void)WriteAHToFile(
    MakeGrid(), options,
    {.summaryFilename = summary, .horizonFilenamePrefix = prefix,
     .visualization = {8, 5}});
  const std::string commonBase =
    prefix + "-common-0-common_failed_found_";
  assert(std::filesystem::exists(commonBase + ".txt"));
  assert(!std::filesystem::exists(commonBase + ".vtp"));
  std::ifstream summaryInput(summary);
  std::ostringstream summaryText;
  summaryText << summaryInput.rdbuf();
  assert(summaryText.str().find(
    "validation.status = provisionally-validated") != std::string::npos);
  std::filesystem::remove(summary);
  std::filesystem::remove(commonBase + ".txt");
  std::filesystem::remove(prefix + "-individual-0-failed.txt");
  std::filesystem::remove(prefix + "-individual-1-found.txt");
  std::filesystem::remove(prefix + "-individual-1-found.vtp");
}

void TestFallbackRequiresPhysicalPositions() {
  AHCommonHorizonSearchOptions options;
  options.blackHoles = {
    {.id = "failed", .mass = 1,
     .attempts = {{{100, 0, 0}, 1.0}}},
    {.id = "found", .mass = 1,
     .attempts = {{{0.1, 0, 0}, 1.3}},
     .physicalPosition = ReferencePoint{0.1, 0, 0}}
  };
  options.candidates.requestedSubsets = {{"failed", "found"}};
  options.lMax = 0;
  options.interpolation = Interpolation();
  bool threw = false;
  try { (void)FindCommonApparentHorizons(MakeGrid(), options); }
  catch (const std::invalid_argument& error) {
    threw = std::string(error.what()).find("physicalPosition is unavailable") !=
      std::string::npos;
  }
  assert(threw);
}

void TestTwoIndividualsProduceOneCommonSearch() {
  AHCommonHorizonSearchOptions options;
  options.blackHoles = {
    {.id = "left", .mass = 1,
     .attempts = {{{-0.1, 0, 0}, 1.3}}},
    {.id = "right", .mass = 1,
     .attempts = {{{0.1, 0, 0}, 1.3}}}
  };
  options.candidates.searchMargin = 3.0;
  options.commonTargets.retryRadiusMultipliers = {1};
  options.lMax = 0;
  options.interpolation = Interpolation();
  options.validation.spectralResidualTolerance = 1.0;
  options.validation.pointwiseExpansionTolerance = 1.0;
  options.validation.enclosureTolerance = 0.5;
  const auto result = FindCommonApparentHorizons(MakeGrid(), options);
  assert(result.individuals.allTargetsFound);
  assert(result.commonHorizons.size() == 1);
  assert(result.commonHorizons[0].id == "common[left+right]");
  assert(result.commonHorizons[0].search.converged);
  assert(result.commonHorizons[0].validation.valid);

  const auto directory = std::filesystem::temp_directory_path();
  const std::string prefix = (directory / "ahf-catalogue-horizon").string();
  const std::string summary = (directory / "ahf-catalogue-summary.txt").string();
  const auto written = WriteAHToFile(
    MakeGrid(), options,
    {.summaryFilename = summary, .horizonFilenamePrefix = prefix,
     .visualization = {8, 5}});
  assert(written.commonHorizons[0].validation.valid);
  assert(std::filesystem::exists(summary));
  assert(std::filesystem::exists(prefix + "-individual-0-left.txt"));
  assert(std::filesystem::exists(prefix + "-individual-1-right.vtp"));
  assert(std::filesystem::exists(prefix + "-common-0-common_left_right_.txt"));
  assert(std::filesystem::exists(prefix + "-common-0-common_left_right_.vtp"));
  std::filesystem::remove(summary);
  for (const auto& suffix : {
         "-individual-0-left.txt", "-individual-0-left.vtp",
         "-individual-1-right.txt", "-individual-1-right.vtp",
         "-common-0-common_left_right_.txt", "-common-0-common_left_right_.vtp"}) {
    std::filesystem::remove(prefix + suffix);
  }

  auto rejectedOptions = options;
  rejectedOptions.blackHoles[0].physicalPosition = ReferencePoint{-100, 0, 0};
  rejectedOptions.blackHoles[1].physicalPosition = ReferencePoint{100, 0, 0};
  const std::string rejectedPrefix =
    (directory / "ahf-rejected-horizon").string();
  const std::string rejectedSummary =
    (directory / "ahf-rejected-summary.txt").string();
  const auto rejected = WriteAHToFile(
    MakeGrid(), rejectedOptions,
    {.summaryFilename = rejectedSummary,
     .horizonFilenamePrefix = rejectedPrefix,
     .visualization = {8, 5}});
  assert(rejected.commonHorizons[0].validation.status ==
         AHCommonHorizonValidationStatus::Rejected);
  const std::string rejectedCommon =
    rejectedPrefix + "-common-0-common_left_right_";
  assert(std::filesystem::exists(rejectedCommon + ".txt"));
  assert(!std::filesystem::exists(rejectedCommon + ".vtp"));
  std::filesystem::remove(rejectedSummary);
  for (const auto& suffix : {
         "-individual-0-left.txt", "-individual-0-left.vtp",
         "-individual-1-right.txt", "-individual-1-right.vtp",
         "-common-0-common_left_right_.txt"}) {
    std::filesystem::remove(rejectedPrefix + suffix);
  }
}

void TestThreeHoleDeterminismAndSummary() {
  AHHorizonCatalogue catalogue;
  catalogue.individuals.targets = {
    {.id = "A", .kind = AHHorizonKind::Individual, .converged = false},
    {.id = "B", .kind = AHHorizonKind::Individual, .converged = true}
  };
  AHCommonHorizonCatalogueEntry ab;
  ab.id = "common[A+B]";
  ab.memberIds = {"A", "B"};
  ab.reason.triggeredPair = true;
  ab.validation.reasons = {"Newton did not converge"};
  AHCommonHorizonCatalogueEntry abc;
  abc.id = "common[A+B+C]";
  abc.memberIds = {"A", "B", "C"};
  abc.reason.connectedComponent = true;
  catalogue.commonHorizons = {ab, abc};
  const auto path = std::filesystem::temp_directory_path() /
    "ahf-common-catalogue-test.txt";
  WriteAHHorizonCatalogueSummary(catalogue, path.string());
  std::ifstream input(path);
  std::ostringstream text;
  text << input.rdbuf();
  assert(text.str().find("common[0].id = \"common[A+B]\"") != std::string::npos);
  assert(text.str().find("Newton did not converge") != std::string::npos);
  std::filesystem::remove(path);
}

}  // namespace
}  // namespace exahype2::ccz4::ahf::tests

int main() {
  using namespace exahype2::ccz4::ahf::tests;
  TestCandidateConversion();
  TestValidationAndEnclosure();
  TestProvisionalAndPhysicalEnclosure();
  TestZeroAndOneBlackHole();
  TestExplicitFallbackAfterIndividualFailure();
  TestFallbackRequiresPhysicalPositions();
  TestTwoIndividualsProduceOneCommonSearch();
  TestThreeHoleDeterminismAndSummary();
}
