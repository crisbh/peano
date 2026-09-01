#include "../../apparent-horizon-finder.h"

#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

namespace exahype2::ccz4::ahf::tests {
namespace {

AHFinderResult SphereResult(double radius, int lMax = 2, bool converged = true) {
  const SpectralBasis basis = BuildSpectralBasis(lMax);
  AHFinderResult result;
  result.solver.coefficients = SphereCoefficients(basis, radius);
  result.solver.converged = converged;
  result.solver.iterations = 3;
  result.solver.residualHistory = {1.0, 0.25, 1.0e-12};
  result.solver.message = converged ? "converged" : "line search failed";
  result.referencePoint = {0.25, -0.5, 0.75};
  result.lMax = lMax;
  result.minRadius = radius;
  result.meanRadius = radius;
  result.maxRadius = radius;
  result.spectralResidualInfinityNorm = 1.0e-12;
  result.maxAbsExpansion = 2.0e-12;
  result.surfaceAvailable = true;
  return result;
}

std::string ReadText(const std::filesystem::path& path) {
  std::ifstream input(path);
  std::ostringstream text;
  text << input.rdbuf();
  return text.str();
}

void TestSummary() {
  const auto directory = std::filesystem::temp_directory_path();
  const auto first = directory / "ahf-output-summary-1.txt";
  const auto second = directory / "ahf-output-summary-2.txt";
  const AHFinderResult result = SphereResult(2.5);
  WriteAHSummary(result, first.string());
  WriteAHSummary(result, second.string());
  const std::string text = ReadText(first);
  assert(text == ReadText(second));
  assert(text.find("solver.converged = true") != std::string::npos);
  assert(text.find("referencePoint = 0.25 -0.5 0.75") != std::string::npos);
  assert(text.find("meanRadius = 2.5") != std::string::npos);
  assert(text.find("residualHistory[2]") != std::string::npos);
  assert(text.find("coefficients[8]") != std::string::npos);

  AHFinderResult failure = result;
  failure.solver.converged = false;
  failure.surfaceAvailable = false;
  failure.solver.message = "controlled failure";
  WriteAHSummary(failure, second.string());
  const std::string failureText = ReadText(second);
  assert(failureText.find("solver.converged = false") != std::string::npos);
  assert(failureText.find("surfaceAvailable = false") != std::string::npos);
  assert(failureText.find("minRadius = unavailable") != std::string::npos);
  assert(failureText.find("controlled failure") != std::string::npos);
  std::filesystem::remove(first);
  std::filesystem::remove(second);
}

void TestDenseSphereGeometry() {
  constexpr int nRho = 12;
  constexpr int nSigma = 7;
  const AHFinderResult result = SphereResult(2.5, 1);
  const AHVisualizationSurface surface =
    BuildAHVisualizationSurface(result, {nRho, nSigma});
  assert(surface.points.size() == 2 + (nSigma - 2) * nRho);
  assert(surface.triangles.size() == 2 * (nSigma - 2) * nRho);
  assert(surface.radius.size() == surface.points.size());
  for (std::size_t i = 0; i < surface.points.size(); ++i) {
    assert(std::abs(surface.radius[i] - 2.5) < 1.0e-12);
    const auto& p = surface.points[i];
    assert(std::isfinite(p[0]) && std::isfinite(p[1]) && std::isfinite(p[2]));
    const double distance = std::sqrt(
      std::pow(p[0] - result.referencePoint.x, 2) +
      std::pow(p[1] - result.referencePoint.y, 2) +
      std::pow(p[2] - result.referencePoint.z, 2));
    assert(std::abs(distance - 2.5) < 1.0e-12);
  }

  std::map<std::pair<std::size_t, std::size_t>, int> edgeUses;
  for (const auto& triangle : surface.triangles) {
    for (std::size_t vertex : triangle) assert(vertex < surface.points.size());
    assert(triangle[0] != triangle[1] && triangle[1] != triangle[2] &&
           triangle[0] != triangle[2]);
    for (int edge = 0; edge < 3; ++edge) {
      std::size_t a = triangle[edge];
      std::size_t b = triangle[(edge + 1) % 3];
      if (a > b) std::swap(a, b);
      ++edgeUses[{a, b}];
    }
  }
  for (const auto& [edge, count] : edgeUses) {
    (void)edge;
    assert(count == 2); // closed surface, including the periodic rho seam
  }

  const auto denser = BuildAHVisualizationSurface(result, {19, 9});
  assert(denser.points.size() != surface.points.size());
  assert(result.lMax == 1); // visualization resolution did not alter the solve

  bool threw = false;
  try { (void)BuildAHVisualizationSurface(result, {2, 7}); }
  catch (const std::invalid_argument&) { threw = true; }
  assert(threw);
}

void TestVTP() {
  const auto path = std::filesystem::temp_directory_path() / "ahf-output-sphere.vtp";
  const AHFinderResult result = SphereResult(1.75);
  WriteAHSurfaceToVTP(result, path.string(), {10, 6});
  assert(std::filesystem::is_regular_file(path));
  const std::string xml = ReadText(path);
  assert(xml.find("<VTKFile type=\"PolyData\"") != std::string::npos);
  assert(xml.find("NumberOfPoints=\"42\"") != std::string::npos);
  assert(xml.find("NumberOfPolys=\"80\"") != std::string::npos);
  assert(xml.find("Name=\"radius\"") != std::string::npos);
  assert(xml.find("Name=\"connectivity\"") != std::string::npos);
  assert(xml.find("Name=\"offsets\"") != std::string::npos);

  AHFinderResult nonconverged = result;
  nonconverged.solver.converged = false;
  bool threw = false;
  try { WriteAHSurfaceToVTP(nonconverged, path.string()); }
  catch (const std::invalid_argument&) { threw = true; }
  assert(threw);
  std::filesystem::remove(path);
}

}  // namespace
}  // namespace exahype2::ccz4::ahf::tests

int main() {
  exahype2::ccz4::ahf::tests::TestSummary();
  exahype2::ccz4::ahf::tests::TestDenseSphereGeometry();
  exahype2::ccz4::ahf::tests::TestVTP();
}
