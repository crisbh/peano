#include "../../apparent-horizon-finder.h"
#include "../../file-input.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <numeric>
#include <string>

namespace exahype2::ccz4::ahf::tests::integration {
namespace {

std::string DatasetPath() {
  const std::string source = __FILE__;
  const std::string directory = source.substr(0, source.find_last_of("/\\"));
  return directory + "/../../ccz4Results/sample_t=0.vtu";
}

SourceVariableMapping SourceMapping() {
  return {
    .conformalMetric = {0, 1, 2, 3, 4, 5},
    .conformalExtrinsicCurvature = {6, 7, 8, 9, 10, 11},
    .traceK = 53,
    .conformalFactor = 54
  };
}

void TestSampleT0UserWorkflow() {
  GridData grid;
  ReadGridFromFile(DatasetPath(), grid, SourceMapping());

  const CartesianPointInterpolator interpolator(
    grid,
    {
      .method = CartesianToSurfaceMethod::RadiusWeightedPolynomial,
      .polynomialDegree = 1,
      .supportRadius = 0.1,
      .rankToleranceMultiplier = 64.0,
      .allowOrderReduction = false
    }
  );
  const SpectralBasis basis = BuildSpectralBasis(2);
  const ReferencePoint referencePoint{0.0, 0.0, 0.0};
  const SpectralNewtonOptions options{
    .tolerance = 1.0e-8,
    .massScale = 1.0,
    .maxIterations = 12,
    .finiteDifferenceStep = 1.0e-6,
    .maxLineSearchSteps = 20,
    .maxRadiusFactor = 6.0
  };

  const SpectralNewtonResult result = SolvePDEforAH(
    basis, interpolator, referencePoint, 0.8, options
  );
  assert(result.converged);

  const SynthesizedField horizon = Synthesize(basis, result.coefficients);
  const double meanRadius = std::accumulate(
    horizon.h.begin(), horizon.h.end(), 0.0
  ) / static_cast<double>(horizon.h.size());

  constexpr double expectedRadius = 0.5;
  constexpr double radiusTolerance = 1.0e-2;
  assert(std::isfinite(meanRadius));
  assert(std::abs(meanRadius - expectedRadius) < radiusTolerance);

  std::cout << "sample_t=0.vtu user workflow: converged in "
            << result.iterations << " iterations; mean AH radius = "
            << meanRadius << ", expected = " << expectedRadius << '\n';
}

}  // namespace
}  // namespace exahype2::ccz4::ahf::tests::integration

int main() {
  exahype2::ccz4::ahf::tests::integration::TestSampleT0UserWorkflow();
}
