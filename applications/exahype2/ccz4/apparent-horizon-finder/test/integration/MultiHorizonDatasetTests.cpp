#include "../../apparent-horizon-finder.h"
#include "../../file-input.h"

#include <cassert>
#include <cmath>
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

void TestSeparatedIndividualHorizons() {
  std::cout << "[AHF] Loading separated D=2 two-puncture dataset...\n";
  GridData grid;
  ReadGridFromFile(
    DatasetPath("two-punctures-D-2-t0.vtu"), grid, SourceMapping());
  assert(grid.numberOfPoints() == 373248);

  AHMultiHorizonOptions options = BaseOptions(0.2);
  options.targets = {
    {
      .id = "individual-plus",
      .kind = AHHorizonKind::Individual,
      .attempts = {{{1.0, 0.0, 0.0}, 0.2},
                   {{1.0, 0.0, 0.0}, 0.12},
                   {{1.0, 0.0, 0.0}, 0.3}}
    },
    {
      .id = "individual-minus",
      .kind = AHHorizonKind::Individual,
      .attempts = {{{-1.0, 0.0, 0.0}, 0.2},
                   {{-1.0, 0.0, 0.0}, 0.12},
                   {{-1.0, 0.0, 0.0}, 0.3}}
    }
  };

  const AHMultiHorizonResult result = FindApparentHorizons(grid, options);
  assert(result.allTargetsFound);
  assert(result.targets.size() == 2);
  for (std::size_t i = 0; i < result.targets.size(); ++i) {
    const AHHorizonTargetResult& target = result.targets[i];
    assert(target.id == options.targets[i].id);
    assert(target.converged);
    assert(target.successfulAttempt == 0);
    assert(target.attempts.size() == 1);
    assert(target.horizon.has_value());
    assert(target.horizon->surfaceAvailable);
    assert(target.horizon->meanRadius > 0.19);
    assert(target.horizon->meanRadius < 0.22);
  }
}

void TestCloseCommonHorizon() {
  std::cout << "[AHF] Loading close D=0.25 two-puncture dataset...\n";
  GridData grid;
  ReadGridFromFile(
    DatasetPath("two-punctures-D-0.25-t0.vtu"), grid, SourceMapping());
  assert(grid.numberOfPoints() == 373248);

  AHMultiHorizonOptions options = BaseOptions(0.15);
  options.targets = {{
    .id = "common",
    .kind = AHHorizonKind::Common,
    .attempts = {{{0.0, 0.0, 0.0}, 0.6},
                 {{0.0, 0.0, 0.0}, 0.9},
                 {{0.0, 0.0, 0.0}, 1.125}}
  }};

  const AHMultiHorizonResult result = FindApparentHorizons(grid, options);
  assert(result.allTargetsFound);
  assert(result.targets.size() == 1);
  const AHHorizonTargetResult& common = result.targets.front();
  assert(common.kind == AHHorizonKind::Common);
  assert(common.successfulAttempt == 0);
  assert(common.attempts.size() == 1);
  assert(common.horizon.has_value());
  assert(common.horizon->surfaceAvailable);
  assert(std::abs(common.horizon->meanRadius - 0.46142331149) < 1.0e-3);
}

}  // namespace
}  // namespace exahype2::ccz4::ahf::tests::integration

int main() {
  exahype2::ccz4::ahf::tests::integration::TestSeparatedIndividualHorizons();
  exahype2::ccz4::ahf::tests::integration::TestCloseCommonHorizon();
}
