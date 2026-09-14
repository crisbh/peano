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
    assert(!attempt.result.solver.residualHistory.empty());
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

}  // namespace
}  // namespace exahype2::ccz4::ahf::tests

int main() {
  exahype2::ccz4::ahf::tests::TestFirstAttemptStopsSearch();
  exahype2::ccz4::ahf::tests::TestRetryOrderingAndArbitraryCount();
  exahype2::ccz4::ahf::tests::TestAllAttemptsFail();
  exahype2::ccz4::ahf::tests::TestValidation();
  exahype2::ccz4::ahf::tests::TestSearchOutput();
}
