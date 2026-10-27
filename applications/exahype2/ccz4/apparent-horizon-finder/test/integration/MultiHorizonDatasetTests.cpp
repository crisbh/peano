#include "../../apparent-horizon-finder.h"
#include "../../file-input.h"

#include <cassert>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>

namespace exahype2::ccz4::ahf::tests::integration {
namespace {

std::string DatasetPath(const std::string& filename) {
  const std::string source = __FILE__;
  const std::string directory = source.substr(0, source.find_last_of("/\\"));
  return directory + "/../../ccz4Results/" + filename;
}

SourceVariableMapping SourceMapping() {
  return {
    .conformalMetric = {0, 1, 2, 3, 4, 5},
    .conformalExtrinsicCurvature = {6, 7, 8, 9, 10, 11},
    .traceK = 53,
    .conformalFactor = 54
  };
}

AHMultiHorizonOptions BaseOptions(double supportRadius) {
  AHMultiHorizonOptions options;
  options.lMax = 6;
  options.interpolation = {
    .method = CartesianToSurfaceMethod::RadiusWeightedPolynomial,
    .polynomialDegree = 3,
    .supportRadius = supportRadius,
    .rankToleranceMultiplier = 64.0,
    .allowOrderReduction = false
  };
  options.solver = {
    .tolerance = 1.0e-8,
    .massScale = 1.0,
    .maxIterations = 40,
    .finiteDifferenceStep = 1.0e-6,
    .maxLineSearchSteps = 20,
    .maxRadiusFactor = 6.0
  };
  return options;
}

AHBinaryMergerOptions MergerOptions(double separation, double supportRadius) {
  const AHMultiHorizonOptions base = BaseOptions(supportRadius);
  AHBinaryMergerOptions options;
  options.components = {{
    {{.id = "individual-plus",
      .attempts = {{{separation / 2.0, 0.0, 0.0}, 0.2},
                   {{separation / 2.0, 0.0, 0.0}, 0.12},
                   {{separation / 2.0, 0.0, 0.0}, 0.3}}}, 0.5},
    {{.id = "individual-minus",
      .attempts = {{{-separation / 2.0, 0.0, 0.0}, 0.2},
                   {{-separation / 2.0, 0.0, 0.0}, 0.12},
                   {{-separation / 2.0, 0.0, 0.0}, 0.3}}}, 0.5}
  }};
  options.commonRadiusFactors = {1.0, 1.5, 2.0};
  options.lMax = base.lMax;
  options.interpolation = base.interpolation;
  options.solver = base.solver;
  return options;
}

AHMultiHorizonOutputOptions DatasetOutput(const std::string& name) {
  return {
    .summaryFilename = name + "-summary.txt",
    .horizonFilenamePrefix = name,
    .visualization = {32, 17}
  };
}

bool KeepTestOutputs() {
  const char* value = std::getenv("AHF_KEEP_TEST_OUTPUTS");
  return value != nullptr && std::string(value) == "1";
}

void RemoveOutput(const AHMultiHorizonOutputOptions& output,
                  const AHMultiHorizonResult& result) {
  if (KeepTestOutputs()) {
    std::cout << "[AHF] Keeping multi-horizon test outputs in "
              << std::filesystem::current_path() << "\n";
    return;
  }
  std::filesystem::remove(output.summaryFilename);
  for (std::size_t i = 0; i < result.targets.size(); ++i) {
    std::string id = result.targets[i].id;
    for (char& character : id) {
      if (!std::isalnum(static_cast<unsigned char>(character)) &&
          character != '-' && character != '_' && character != '.') {
        character = '_';
      }
    }
    const std::string base = output.horizonFilenamePrefix + "-" +
      std::to_string(i) + "-" + id;
    std::filesystem::remove(base + ".txt");
    std::filesystem::remove(base + ".vtp");
  }
}

void TestSeparatedIndividualHorizons() {
  std::cout << "[AHF] Loading separated D=2 two-puncture dataset...\n";
  GridData grid;
  ReadGridFromFile(
    DatasetPath("two-punctures-D-2-t0.vtu"), grid, SourceMapping());
  assert(grid.numberOfPoints() == 373248);

  const AHBinaryMergerOptions options = MergerOptions(2.0, 0.2);

  const AHMultiHorizonOutputOptions output =
    DatasetOutput("ahf-two-punctures-separated");
  const AHBinaryMergerResult merger = WriteAHToFile(grid, options, output);
  const AHMultiHorizonResult& result = merger.horizons;
  assert(!merger.plan.commonSearchEnabled);
  assert(result.allTargetsFound);
  assert(result.targets.size() == 2);
  for (std::size_t i = 0; i < result.targets.size(); ++i) {
    const AHHorizonTargetResult& target = result.targets[i];
    assert(target.id == options.components[i].target.id);
    assert(target.converged);
    assert(target.successfulAttempt == 0);
    assert(target.attempts.size() == 1);
    assert(target.horizon.has_value());
    assert(target.horizon->surfaceAvailable);
    assert(target.horizon->meanRadius > 0.19);
    assert(target.horizon->meanRadius < 0.22);
    const std::string base = output.horizonFilenamePrefix + "-" +
      std::to_string(i) + "-" + target.id;
    assert(std::filesystem::is_regular_file(base + ".txt"));
    assert(std::filesystem::is_regular_file(base + ".vtp"));
  }
  assert(std::filesystem::is_regular_file(output.summaryFilename));
  RemoveOutput(output, result);
}

void TestCloseCommonHorizon() {
  std::cout << "[AHF] Loading close D=0.25 two-puncture dataset...\n";
  GridData grid;
  ReadGridFromFile(
    DatasetPath("two-punctures-D-0.25-t0.vtu"), grid, SourceMapping());
  assert(grid.numberOfPoints() == 373248);

  const AHBinaryMergerOptions options = MergerOptions(0.25, 0.15);

  const AHMultiHorizonOutputOptions output =
    DatasetOutput("ahf-two-punctures-close");
  const AHBinaryMergerResult merger = WriteAHToFile(grid, options, output);
  const AHMultiHorizonResult& result = merger.horizons;
  assert(merger.plan.commonSearchEnabled);
  assert(std::abs(merger.plan.separation - 0.25) < 1.0e-14);
  assert(std::abs(merger.plan.commonReferencePoint.x) < 1.0e-14);
  assert(result.allTargetsFound);
  assert(result.targets.size() == 3);
  const AHHorizonTargetResult& common = result.targets.front();
  assert(common.kind == AHHorizonKind::Common);
  assert(common.successfulAttempt == 0);
  assert(common.attempts.size() == 1);
  assert(common.horizon.has_value());
  assert(common.horizon->surfaceAvailable);
  assert(std::abs(common.horizon->meanRadius - 0.46142331149) < 1.0e-3);
  assert(std::filesystem::is_regular_file(output.summaryFilename));
  assert(std::filesystem::is_regular_file(
    output.horizonFilenamePrefix + "-0-common.txt"));
  assert(std::filesystem::is_regular_file(
    output.horizonFilenamePrefix + "-0-common.vtp"));
  RemoveOutput(output, result);
}

}  // namespace
}  // namespace exahype2::ccz4::ahf::tests::integration

int main() {
  exahype2::ccz4::ahf::tests::integration::TestSeparatedIndividualHorizons();
  exahype2::ccz4::ahf::tests::integration::TestCloseCommonHorizon();
}
