#include "../../apparent-horizon-finder.h"
#include "../../file-input.h"

#include <cassert>
#include <chrono>
#include <cmath>
#include <iostream>
#include <string>

namespace exahype2::ccz4::ahf::tests::integration {
namespace {

using Clock = std::chrono::steady_clock;

double Seconds(Clock::time_point begin, Clock::time_point end) {
  return std::chrono::duration<double>(end - begin).count();
}

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
  const auto totalStart = Clock::now();

  std::cout << "[AHF] Building the spectral basis and initial guess..."
            << std::endl;
  const auto setupStart = Clock::now();
  const SpectralBasis basis = BuildSpectralBasis(2);
  const ReferencePoint referencePoint{0.0, 0.0, 0.0};
  constexpr double initialAHRadiusGuess = 0.9;
  const auto setupFinish = Clock::now();
  std::cout << "[AHF] Basis and initial-guess setup took "
            << Seconds(setupStart, setupFinish) << " s." << std::endl;

  std::cout << "[AHF] Loading Cartesian data from sample_t=0.vtu..."
            << std::endl;
  const auto loadStart = Clock::now();
  GridData grid;
  ReadGridFromFile(DatasetPath(), grid, SourceMapping());
  const auto loadFinish = Clock::now();
  std::cout << "[AHF] Loaded " << grid.numberOfPoints()
            << " Cartesian cells in " << Seconds(loadStart, loadFinish)
            << " s." << std::endl;

  std::cout << "[AHF] Constructing the persistent Cartesian interpolator..."
            << std::endl;
  const auto interpolatorStart = Clock::now();
  const CartesianPointInterpolator interpolator(
    grid,
    {
      .method = CartesianToSurfaceMethod::RadiusWeightedPolynomial,
      .polynomialDegree = 3,
      .supportRadius = 0.1,
      .rankToleranceMultiplier = 64.0,
      .allowOrderReduction = false
    }
  );
  const auto interpolatorFinish = Clock::now();
  std::cout << "[AHF] Interpolator construction took "
            << Seconds(interpolatorStart, interpolatorFinish) << " s."
            << std::endl;

  const SpectralNewtonOptions options{
    .tolerance = 1.0e-8,
    .massScale = 1.0,
    .maxIterations = 12,
    .finiteDifferenceStep = 1.0e-6,
    .maxLineSearchSteps = 20,
    .maxRadiusFactor = 6.0
  };

  std::cout << "[AHF] Iterating the apparent-horizon equation..."
            << std::endl;
  const auto solveStart = Clock::now();
  const SpectralNewtonResult result = SolvePDEforAH(
    basis, interpolator, referencePoint, initialAHRadiusGuess, options
  );
  const auto solveFinish = Clock::now();
  assert(result.converged);
  std::cout << "[AHF] Iteration converged after " << result.iterations
            << " Newton iterations in " << Seconds(solveStart, solveFinish)
            << " s." << std::endl;

  std::cout << "[AHF] Reconstructing and validating the final horizon..."
            << std::endl;
  const auto reconstructionStart = Clock::now();
  const AHFinderResult horizon = BuildAHFinderResult(
    basis, interpolator, referencePoint, result
  );

  constexpr double expectedRadius = 0.5;
  constexpr double radiusTolerance = 1.0e-2;

  std::cout << "[AHF] Finished: mean AH radius = " << horizon.meanRadius
            << ", expected = " << expectedRadius << std::endl;

  const auto reconstructionFinish = Clock::now();
  std::cout << "[AHF] Final reconstruction and validation took "
            << Seconds(reconstructionStart, reconstructionFinish) << " s."
            << std::endl
            << "[AHF] Total wall-clock time: "
            << Seconds(totalStart, reconstructionFinish) << " s."
            << std::endl;

  assert(horizon.solver.converged);
  assert(horizon.surfaceAvailable);
  assert(horizon.referencePoint.x == referencePoint.x);
  assert(horizon.referencePoint.y == referencePoint.y);
  assert(horizon.referencePoint.z == referencePoint.z);
  assert(std::isfinite(horizon.minRadius));
  assert(std::isfinite(horizon.meanRadius));
  assert(std::isfinite(horizon.maxRadius));
  assert(horizon.minRadius <= horizon.meanRadius);
  assert(horizon.meanRadius <= horizon.maxRadius);
  assert(std::isfinite(horizon.maxAbsExpansion));
  assert(std::isfinite(horizon.spectralResidualInfinityNorm));
  assert(!horizon.solver.coefficients.empty());
  assert(std::abs(horizon.meanRadius - expectedRadius) < radiusTolerance);

}

}  // namespace
}  // namespace exahype2::ccz4::ahf::tests::integration

int main() {
  exahype2::ccz4::ahf::tests::integration::TestSampleT0UserWorkflow();
}
