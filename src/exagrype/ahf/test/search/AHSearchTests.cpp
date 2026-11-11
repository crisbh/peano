#include "../../apparent-horizon-finder.h"

#include <array>
#include <cassert>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

namespace exahype2::ccz4::ahf::tests {
namespace {

GridData MakeGrid() {
  GridData grid;
  for (int ix = -2; ix <= 2; ++ix) {
    for (int iy = -2; iy <= 2; ++iy) {
      for (int iz = -2; iz <= 2; ++iz) {
        grid.coordinates.push_back({double(ix), double(iy), double(iz)});
        const std::array<double, NumberOfCompactFields> row{
          1.0, 0.0, 0.0, 1.0, 0.0, 1.0,
          1.0, 0.0, 0.0, 1.0, 0.0, 1.0,
          3.0, 1.0
        };
        grid.values.insert(grid.values.end(), row.begin(), row.end());
      }
    }
  }
  return grid;
}

AHSearchOptions SearchOptions(std::vector<AHSearchAttempt> attempts) {
  AHSearchOptions options;
  options.attempts = std::move(attempts);
  options.lMax = 0;
  options.interpolation = {
    .method = CartesianToSurfaceMethod::UnweightedPolynomial,
    .numberOfNeighbours = 64,
    .polynomialDegree = 1,
    .allowOrderReduction = false
  };
  return options;
}

std::string ReadText(const std::filesystem::path& path) {
  std::ifstream input(path);
  std::ostringstream text;
  text << input.rdbuf();
  return text.str();
}

void TestFirstAttemptStopsSearch() {
  const AHSearchResult result = FindApparentHorizon(
    MakeGrid(),
    SearchOptions({{{0.0, 0.0, 0.0}, 1.3}, {{100.0, 0.0, 0.0}, 1.0}})
  );
  assert(result.converged);
  assert(result.successfulAttempt == 0);
  assert(result.attempts.size() == 1);
  assert(result.horizon.has_value());
  assert(result.horizon->solver.converged);
}

void TestRetryOrderingAndArbitraryCount() {
  std::vector<AHSearchAttempt> attempts;
  for (int i = 0; i < 5; ++i) {
    attempts.push_back({{100.0 + i, 0.0, 0.0}, 1.0});
  }
  attempts.push_back({{0.0, 0.0, 0.0}, 1.3});
  const AHSearchResult result = FindApparentHorizon(
    MakeGrid(), SearchOptions(attempts));
  assert(result.converged);
  assert(result.successfulAttempt == 5);
  assert(result.attempts.size() == 6);
  for (std::size_t i = 0; i < 5; ++i) {
    assert(result.attempts[i].referencePoint.x == 100.0 + double(i));
    assert(!result.attempts[i].succeeded);
    assert(!result.attempts[i].result.solver.message.empty());
  }
  assert(result.attempts[5].succeeded);
}

void TestAllAttemptsFail() {
  const AHSearchResult result = FindApparentHorizon(
    MakeGrid(),
    SearchOptions({{{100.0, 0.0, 0.0}, 1.0},
                   {{101.0, 0.0, 0.0}, 1.1},
                   {{102.0, 0.0, 0.0}, 1.2}})
  );
  assert(!result.converged);
  assert(result.successfulAttempt == AHSearchResult::NoSuccessfulAttempt);
  assert(!result.horizon.has_value());
  assert(result.attempts.size() == 3);
  for (const AHSearchAttemptResult& attempt : result.attempts) {
    assert(!attempt.succeeded);
    assert(!attempt.result.solver.message.empty());
  }
}

void TestValidation() {
  bool threw = false;
  try { (void)FindApparentHorizon(MakeGrid(), SearchOptions({})); }
  catch (const std::invalid_argument&) { threw = true; }
  assert(threw);

  threw = false;
  try {
    (void)FindApparentHorizon(
      MakeGrid(), SearchOptions({{{0.0, 0.0, 0.0}, 0.0}}));
  } catch (const std::invalid_argument&) { threw = true; }
  assert(threw);
}

void TestSearchOutput() {
  const GridData grid = MakeGrid();
  const auto directory = std::filesystem::temp_directory_path();
  const auto summary = directory / "ahf-search-summary.txt";
  const auto surface = directory / "ahf-search-surface.vtp";
  std::filesystem::remove(summary);
  std::filesystem::remove(surface);

  AHSearchResult failed = WriteAHToFile(
    grid,
    SearchOptions({{{100.0, 0.0, 0.0}, 1.0},
                   {{101.0, 0.0, 0.0}, 1.1}}),
    {.summaryFilename = summary.string(), .surfaceFilename = surface.string()}
  );
  assert(!failed.converged);
  assert(std::filesystem::is_regular_file(summary));
  assert(!std::filesystem::exists(surface));
  const std::string failedText = ReadText(summary);
  assert(failedText.find("attempts.count = 2") != std::string::npos);
  assert(failedText.find("attempts[0].referencePoint = 100 0 0") !=
         std::string::npos);
  assert(failedText.find("attempts[1].referencePoint = 101 0 0") !=
         std::string::npos);

  std::filesystem::remove(summary);
  AHSearchResult success = WriteAHToFile(
    grid,
    SearchOptions({{{100.0, 0.0, 0.0}, 1.0},
                   {{0.0, 0.0, 0.0}, 1.3}}),
    {.summaryFilename = summary.string(),
     .surfaceFilename = surface.string(),
     .visualization = {9, 6}}
  );
  assert(success.converged);
  assert(success.successfulAttempt == 1);
  assert(std::filesystem::is_regular_file(summary));
  assert(std::filesystem::is_regular_file(surface));
  const std::string successText = ReadText(summary);
  assert(successText.find("AHSearchSummaryVersion = 1") != std::string::npos);
  assert(successText.find("successfulAttempt = 1") != std::string::npos);
  assert(successText.find("attempts.count = 2") != std::string::npos);
  assert(successText.find("attempts[0].succeeded = false") != std::string::npos);
  assert(successText.find("attempts[1].succeeded = true") != std::string::npos);
  assert(successText.find("horizon.meanRadius = ") != std::string::npos);
  assert(ReadText(surface).find("NumberOfPoints=\"38\"") != std::string::npos);
  std::filesystem::remove(summary);
  std::filesystem::remove(surface);
}

AHMultiHorizonOptions MultiOptions(std::vector<AHHorizonTarget> targets) {
  AHMultiHorizonOptions options;
  options.targets = std::move(targets);
  options.lMax = 0;
  options.interpolation = {
    .method = CartesianToSurfaceMethod::UnweightedPolynomial,
    .numberOfNeighbours = 64,
    .polynomialDegree = 1,
    .allowOrderReduction = false
  };
  return options;
}

void TestPhysicalTargetGrouping() {
  const AHMultiHorizonResult result = FindApparentHorizons(
    MakeGrid(),
    MultiOptions({
      {
        .id = "BH-plus",
        .kind = AHHorizonKind::Individual,
        .attempts = {{{100.0, 0.0, 0.0}, 1.0},
                     {{0.25, 0.0, 0.0}, 1.3},
                     {{0.0, 0.0, 0.0}, 1.3}}
      },
      {
        .id = "BH-minus",
        .kind = AHHorizonKind::Individual,
        .attempts = {{{-0.25, 0.0, 0.0}, 1.3}}
      }
    })
  );
  assert(result.anyHorizonFound);
  assert(result.allTargetsFound);
  assert(result.targets.size() == 2);
  assert(result.targets[0].id == "BH-plus");
  assert(result.targets[0].successfulAttempt == 1);
  assert(result.targets[0].attempts.size() == 2); // third retry was not run
  assert(!result.targets[0].attempts[0].succeeded);
  assert(result.targets[0].attempts[1].succeeded);
  assert(result.targets[1].id == "BH-minus");
  assert(result.targets[1].converged); // first target success did not stop it
  assert(result.targets[1].attempts.size() == 1);
}

void TestPhysicalTargetPartialSuccess() {
  const AHMultiHorizonResult result = FindApparentHorizons(
    MakeGrid(),
    MultiOptions({
      {.id = "missing", .attempts = {{{100.0, 0.0, 0.0}, 1.0}}},
      {.id = "common", .kind = AHHorizonKind::Common,
       .attempts = {{{0.0, 0.0, 0.0}, 1.3}}}
    })
  );
  assert(result.anyHorizonFound);
  assert(!result.allTargetsFound);
  assert(!result.targets[0].converged);
  assert(!result.targets[0].horizon.has_value());
  assert(result.targets[1].converged);
  assert(result.targets[1].horizon.has_value());
  assert(result.targets[1].kind == AHHorizonKind::Common);
}

void TestPhysicalTargetValidation() {
  bool threw = false;
  try { (void)FindApparentHorizons(MakeGrid(), MultiOptions({})); }
  catch (const std::invalid_argument&) { threw = true; }
  assert(threw);

  threw = false;
  try {
    (void)FindApparentHorizons(
      MakeGrid(),
      MultiOptions({{.id = "same", .attempts = {{{0, 0, 0}, 1.0}}},
                    {.id = "same", .attempts = {{{0, 0, 0}, 1.0}}}}));
  } catch (const std::invalid_argument&) { threw = true; }
  assert(threw);
}

void TestSpectralAndPreviousSurfaceAttempts() {
  const SpectralBasis basis = BuildSpectralBasis(0);
  const std::vector<double> coefficients =
    MakeSphericalInitialCoefficients(basis, 1.3);
  AHSearchOptions options = SearchOptions({{
    .referencePoint = {0.0, 0.0, 0.0},
    .initialRadius = 0.0,
    .initialCoefficients = coefficients
  }});
  const AHSearchResult first = FindApparentHorizon(MakeGrid(), options);
  assert(first.converged);
  assert(first.attempts.size() == 1);
  assert(first.attempts[0].initialCoefficients == coefficients);

  const AHSearchAttempt continuation =
    MakePreviousSurfaceAttempt(*first.horizon);
  assert(continuation.referencePoint.x == first.horizon->referencePoint.x);
  assert(continuation.initialCoefficients == first.horizon->solver.coefficients);
  const AHSearchResult second = FindApparentHorizon(
    MakeGrid(), SearchOptions({continuation}));
  assert(second.converged);

  options.attempts[0].initialCoefficients.push_back(0.0);
  bool rejected = false;
  try { (void)FindApparentHorizon(MakeGrid(), options); }
  catch (const std::invalid_argument&) { rejected = true; }
  assert(rejected);

  AHSearchOptions distortedOptions = SearchOptions({});
  distortedOptions.lMax = 2;
  const SpectralBasis distortedBasis = BuildSpectralBasis(2);
  std::vector<double> radii(distortedBasis.nPoints);
  for (std::size_t p = 0; p < radii.size(); ++p) {
    radii[p] = 1.3 + 0.05 * EvalMode(
      2, 0, distortedBasis.rho[p], distortedBasis.sigma[p]).Y;
  }
  distortedOptions.attempts = {{
    .referencePoint = {0.0, 0.0, 0.0},
    .initialCoefficients =
      MakeDistortedInitialCoefficients(distortedBasis, radii)
  }};
  const AHSearchResult distorted =
    FindApparentHorizon(MakeGrid(), distortedOptions);
  assert(distorted.converged);
  assert(distorted.horizon->solver.coefficients.size() ==
         distortedBasis.nModes);
}

void TestPhysicalTargetOutput() {
  const GridData grid = MakeGrid();
  const auto directory = std::filesystem::temp_directory_path();
  const auto aggregate = directory / "ahf-multi-horizon-summary.txt";
  const auto prefix = directory / "ahf-multi-horizon";
  const auto missingSummary = directory / "ahf-multi-horizon-0-missing.txt";
  const auto missingSurface = directory / "ahf-multi-horizon-0-missing.vtp";
  const auto plusSummary = directory / "ahf-multi-horizon-1-BH_plus.txt";
  const auto plusSurface = directory / "ahf-multi-horizon-1-BH_plus.vtp";
  const auto commonSummary = directory / "ahf-multi-horizon-2-common.txt";
  const auto commonSurface = directory / "ahf-multi-horizon-2-common.vtp";
  for (const auto& path : {aggregate, missingSummary, missingSurface,
                           plusSummary, plusSurface, commonSummary,
                           commonSurface}) {
    std::filesystem::remove(path);
  }

  const AHMultiHorizonResult result = WriteAHToFile(
    grid,
    MultiOptions({
      {.id = "missing", .attempts = {{{100.0, 0.0, 0.0}, 1.0}}},
      {.id = "BH plus", .attempts = {{{0.25, 0.0, 0.0}, 1.3}}},
      {.id = "common", .kind = AHHorizonKind::Common,
       .attempts = {{{0.0, 0.0, 0.0}, 1.3}}}
    }),
    {.summaryFilename = aggregate.string(),
     .horizonFilenamePrefix = prefix.string(),
     .visualization = {9, 6}}
  );

  assert(result.anyHorizonFound);
  assert(!result.allTargetsFound);
  assert(std::filesystem::is_regular_file(aggregate));
  assert(!std::filesystem::exists(missingSummary));
  assert(!std::filesystem::exists(missingSurface));
  assert(std::filesystem::is_regular_file(plusSummary));
  assert(std::filesystem::is_regular_file(plusSurface));
  assert(std::filesystem::is_regular_file(commonSummary));
  assert(std::filesystem::is_regular_file(commonSurface));

  const std::string aggregateText = ReadText(aggregate);
  assert(aggregateText.find("AHMultiHorizonSummaryVersion = 1") !=
         std::string::npos);
  assert(aggregateText.find("anyHorizonFound = true") != std::string::npos);
  assert(aggregateText.find("allTargetsFound = false") != std::string::npos);
  assert(aggregateText.find("targets.count = 3") != std::string::npos);
  assert(aggregateText.find("targets[0].id = \"missing\"") !=
         std::string::npos);
  assert(aggregateText.find("targets[0].solver.message") == std::string::npos);
  assert(aggregateText.find("targets[0].attempts[0].solver.message") !=
         std::string::npos);
  assert(aggregateText.find("targets[1].id = \"BH plus\"") !=
         std::string::npos);
  assert(aggregateText.find("targets[2].kind = common") != std::string::npos);
  assert(ReadText(plusSurface).find("NumberOfPoints=\"38\"") !=
         std::string::npos);

  for (const auto& path : {aggregate, plusSummary, plusSurface, commonSummary,
                           commonSurface}) {
    std::filesystem::remove(path);
  }
}

}  // namespace
}  // namespace exahype2::ccz4::ahf::tests

int main() {
  exahype2::ccz4::ahf::tests::TestFirstAttemptStopsSearch();
  exahype2::ccz4::ahf::tests::TestRetryOrderingAndArbitraryCount();
  exahype2::ccz4::ahf::tests::TestAllAttemptsFail();
  exahype2::ccz4::ahf::tests::TestValidation();
  exahype2::ccz4::ahf::tests::TestSearchOutput();
  exahype2::ccz4::ahf::tests::TestPhysicalTargetGrouping();
  exahype2::ccz4::ahf::tests::TestPhysicalTargetPartialSuccess();
  exahype2::ccz4::ahf::tests::TestPhysicalTargetValidation();
  exahype2::ccz4::ahf::tests::TestSpectralAndPreviousSurfaceAttempts();
  exahype2::ccz4::ahf::tests::TestPhysicalTargetOutput();
}
