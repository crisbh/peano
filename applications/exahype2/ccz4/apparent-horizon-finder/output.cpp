#include "output.h"

#include "harmonics.h"
#include "iteration.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numbers>
#include <stdexcept>

namespace exahype2::ccz4::ahf {
namespace {

void ValidateVisualizationOptions(const AHVisualizationOptions& options) {
  if (options.numberOfRhoPoints < 3) {
    throw std::invalid_argument(
      "AH visualization requires at least 3 rho points");
  }
  if (options.numberOfSigmaPoints < 3) {
    throw std::invalid_argument(
      "AH visualization requires at least 3 sigma points including the poles");
  }
}

double EvaluateRadius(const AHFinderResult& result, double rho, double sigma,
                      bool atPole) {
  double radius = 0.0;
  const auto modes = EnumerateModes(result.lMax);
  if (modes.size() != result.solver.coefficients.size()) {
    throw std::invalid_argument(
      "AH visualization coefficient count does not match lMax");
  }
  for (std::size_t j = 0; j < modes.size(); ++j) {
    const auto [l, m] = modes[j];
    double y = 0.0;
    if (atPole) {
      // Only m=0 survives at a pole. Avoid derivative recurrences in EvalMode,
      // which intentionally assume collocation points strictly inside the poles.
      if (m == 0) {
        const double sign = sigma > 0.0 || l % 2 == 0 ? 1.0 : -1.0;
        y = sign * std::sqrt((2.0 * l + 1.0) /
                             (4.0 * std::numbers::pi));
      }
    } else {
      y = EvalMode(l, m, rho, sigma).Y;
    }
    radius += result.solver.coefficients[j] * y;
  }
  if (!std::isfinite(radius) || radius <= 0.0) {
    throw std::runtime_error(
      "AH visualization produced a non-positive or non-finite radius");
  }
  return radius;
}

void AppendPoint(AHVisualizationSurface& surface, const ReferencePoint& origin,
                 double radius, double rho, double sigma) {
  const double cosSigma = std::cos(sigma);
  surface.points.push_back({
    origin.x + radius * cosSigma * std::cos(rho),
    origin.y + radius * cosSigma * std::sin(rho),
    origin.z + radius * std::sin(sigma)
  });
  surface.radius.push_back(radius);
}

std::ofstream OpenOutput(const std::string& filename) {
  std::ofstream stream(filename);
  if (!stream) {
    throw std::runtime_error("Unable to open AH output file: " + filename);
  }
  return stream;
}

}  // namespace

void ValidateAHOutputOptions(const AHOutputOptions& options) {
  if (options.summaryFilename.empty()) {
    throw std::invalid_argument(
      "ValidateAHOutputOptions: summaryFilename must not be empty");
  }
  if (options.surfaceFilename.empty()) {
    throw std::invalid_argument(
      "ValidateAHOutputOptions: surfaceFilename must not be empty");
  }
  if (options.summaryFilename == options.surfaceFilename) {
    throw std::invalid_argument(
      "ValidateAHOutputOptions: summary and surface filenames must differ");
  }
  ValidateVisualizationOptions(options.visualization);
}

void WriteAHSummary(const AHFinderResult& result, const std::string& filename) {
  std::ofstream output = OpenOutput(filename);
  output << std::setprecision(std::numeric_limits<double>::max_digits10)
         << std::boolalpha
         << "AHFinderSummaryVersion = 1\n"
         << "solver.converged = " << result.solver.converged << '\n'
         << "surfaceAvailable = " << result.surfaceAvailable << '\n'
         << "solver.iterations = " << result.solver.iterations << '\n'
         << "solver.message = " << std::quoted(result.solver.message) << '\n'
         << "referencePoint = " << result.referencePoint.x << ' '
         << result.referencePoint.y << ' ' << result.referencePoint.z << '\n'
         << "lMax = " << result.lMax << '\n';

  if (result.surfaceAvailable) {
    output << "minRadius = " << result.minRadius << '\n'
           << "meanRadius = " << result.meanRadius << '\n'
           << "maxRadius = " << result.maxRadius << '\n'
           << "spectralResidualInfinityNorm = "
           << result.spectralResidualInfinityNorm << '\n'
           << "maxAbsExpansion = " << result.maxAbsExpansion << '\n';
  } else {
    output << "minRadius = unavailable\n"
           << "meanRadius = unavailable\n"
           << "maxRadius = unavailable\n"
           << "spectralResidualInfinityNorm = unavailable\n"
           << "maxAbsExpansion = unavailable\n";
  }

  output << "residualHistory.count = " << result.solver.residualHistory.size()
         << '\n';
  for (std::size_t i = 0; i < result.solver.residualHistory.size(); ++i) {
    output << "residualHistory[" << i << "] = "
           << result.solver.residualHistory[i] << '\n';
  }
  output << "coefficients.count = " << result.solver.coefficients.size() << '\n';
  for (std::size_t i = 0; i < result.solver.coefficients.size(); ++i) {
    output << "coefficients[" << i << "] = "
           << result.solver.coefficients[i] << '\n';
  }
  if (!output) {
    throw std::runtime_error("Failed while writing AH summary: " + filename);
  }
}

void WriteAHSearchSummary(const AHSearchResult& result,
                          const std::string& filename) {
  std::ofstream output = OpenOutput(filename);
  output << std::setprecision(std::numeric_limits<double>::max_digits10)
         << std::boolalpha
         << "AHSearchSummaryVersion = 1\n"
         << "converged = " << result.converged << '\n';
  if (result.converged) {
    output << "successfulAttempt = " << result.successfulAttempt << '\n';
  } else {
    output << "successfulAttempt = unavailable\n";
  }
  output << "attempts.count = " << result.attempts.size() << '\n';
  for (std::size_t i = 0; i < result.attempts.size(); ++i) {
    const AHSearchAttemptResult& attempt = result.attempts[i];
    const std::string prefix = "attempts[" + std::to_string(i) + "].";
    output << prefix << "referencePoint = " << attempt.referencePoint.x << ' '
           << attempt.referencePoint.y << ' ' << attempt.referencePoint.z << '\n'
           << prefix << "initialRadius = " << attempt.initialRadius << '\n'
           << prefix << "succeeded = " << attempt.succeeded << '\n'
           << prefix << "solver.converged = "
           << attempt.result.solver.converged << '\n'
           << prefix << "surfaceAvailable = "
           << attempt.result.surfaceAvailable << '\n'
           << prefix << "solver.iterations = "
           << attempt.result.solver.iterations << '\n'
           << prefix << "solver.message = "
           << std::quoted(attempt.result.solver.message) << '\n'
           << prefix << "residualHistory.count = "
           << attempt.result.solver.residualHistory.size() << '\n';
    for (std::size_t j = 0;
         j < attempt.result.solver.residualHistory.size(); ++j) {
      output << prefix << "residualHistory[" << j << "] = "
             << attempt.result.solver.residualHistory[j] << '\n';
    }
  }
  if (result.converged && result.horizon.has_value()) {
    const AHFinderResult& horizon = *result.horizon;
    output << "horizon.referencePoint = " << horizon.referencePoint.x << ' '
           << horizon.referencePoint.y << ' ' << horizon.referencePoint.z << '\n'
           << "horizon.lMax = " << horizon.lMax << '\n'
           << "horizon.minRadius = " << horizon.minRadius << '\n'
           << "horizon.meanRadius = " << horizon.meanRadius << '\n'
           << "horizon.maxRadius = " << horizon.maxRadius << '\n'
           << "horizon.spectralResidualInfinityNorm = "
           << horizon.spectralResidualInfinityNorm << '\n'
           << "horizon.maxAbsExpansion = " << horizon.maxAbsExpansion << '\n'
           << "horizon.coefficients.count = "
           << horizon.solver.coefficients.size() << '\n';
    for (std::size_t i = 0; i < horizon.solver.coefficients.size(); ++i) {
      output << "horizon.coefficients[" << i << "] = "
             << horizon.solver.coefficients[i] << '\n';
    }
  }
  if (!output) {
    throw std::runtime_error(
      "Failed while writing AH search summary: " + filename);
  }
}

AHVisualizationSurface BuildAHVisualizationSurface(
  const AHFinderResult& result, const AHVisualizationOptions& options) {
  ValidateVisualizationOptions(options);
  if (!result.surfaceAvailable) {
    throw std::invalid_argument(
      "Cannot build AH visualization: final surface is unavailable");
  }

  AHVisualizationSurface surface;
  const std::size_t nRho = static_cast<std::size_t>(options.numberOfRhoPoints);
  const std::size_t nSigma =
    static_cast<std::size_t>(options.numberOfSigmaPoints);
  const std::size_t ringCount = nSigma - 2;
  surface.points.reserve(2 + ringCount * nRho);
  surface.radius.reserve(2 + ringCount * nRho);
  surface.triangles.reserve(2 * nRho * ringCount);

  const double southSigma = -std::numbers::pi / 2.0;
  const double northSigma = std::numbers::pi / 2.0;
  AppendPoint(surface, result.referencePoint,
              EvaluateRadius(result, 0.0, southSigma, true), 0.0, southSigma);

  for (std::size_t i = 1; i + 1 < nSigma; ++i) {
    const double sigma = southSigma + std::numbers::pi *
      static_cast<double>(i) / static_cast<double>(nSigma - 1);
    for (std::size_t j = 0; j < nRho; ++j) {
      const double rho = 2.0 * std::numbers::pi * static_cast<double>(j) /
                         static_cast<double>(nRho);
      AppendPoint(surface, result.referencePoint,
                  EvaluateRadius(result, rho, sigma, false), rho, sigma);
    }
  }
  const std::size_t north = surface.points.size();
  AppendPoint(surface, result.referencePoint,
              EvaluateRadius(result, 0.0, northSigma, true), 0.0, northSigma);

  const std::size_t firstRing = 1;
  for (std::size_t j = 0; j < nRho; ++j) {
    const std::size_t next = (j + 1) % nRho;
    surface.triangles.push_back({0, firstRing + next, firstRing + j});
  }
  for (std::size_t i = 0; i + 1 < ringCount; ++i) {
    const std::size_t lower = firstRing + i * nRho;
    const std::size_t upper = lower + nRho;
    for (std::size_t j = 0; j < nRho; ++j) {
      const std::size_t next = (j + 1) % nRho;
      surface.triangles.push_back({lower + j, lower + next, upper + j});
      surface.triangles.push_back({lower + next, upper + next, upper + j});
    }
  }
  const std::size_t lastRing = firstRing + (ringCount - 1) * nRho;
  for (std::size_t j = 0; j < nRho; ++j) {
    const std::size_t next = (j + 1) % nRho;
    surface.triangles.push_back({lastRing + j, lastRing + next, north});
  }
  return surface;
}

void WriteAHSurfaceToVTP(const AHFinderResult& result,
                         const std::string& filename,
                         const AHVisualizationOptions& options) {
  if (!result.solver.converged) {
    throw std::invalid_argument(
      "Cannot write normal AH visualization: solver did not converge");
  }
  const AHVisualizationSurface surface =
    BuildAHVisualizationSurface(result, options);
  std::ofstream output = OpenOutput(filename);
  output << std::setprecision(std::numeric_limits<double>::max_digits10)
         << "<?xml version=\"1.0\"?>\n"
         << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
         << "  <PolyData>\n"
         << "    <Piece NumberOfPoints=\"" << surface.points.size()
         << "\" NumberOfVerts=\"0\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\""
         << surface.triangles.size() << "\">\n"
         << "      <PointData Scalars=\"radius\">\n"
         << "        <DataArray type=\"Float64\" Name=\"radius\" format=\"ascii\">\n          ";
  for (double radius : surface.radius) output << radius << ' ';
  output << "\n        </DataArray>\n      </PointData>\n"
         << "      <CellData/>\n      <Points>\n"
         << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n          ";
  for (const auto& point : surface.points) {
    output << point[0] << ' ' << point[1] << ' ' << point[2] << ' ';
  }
  output << "\n        </DataArray>\n      </Points>\n      <Polys>\n"
         << "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n          ";
  for (const auto& triangle : surface.triangles) {
    output << triangle[0] << ' ' << triangle[1] << ' ' << triangle[2] << ' ';
  }
  output << "\n        </DataArray>\n"
         << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n          ";
  for (std::size_t i = 1; i <= surface.triangles.size(); ++i) output << 3 * i << ' ';
  output << "\n        </DataArray>\n      </Polys>\n"
         << "    </Piece>\n  </PolyData>\n</VTKFile>\n";
  if (!output) {
    throw std::runtime_error("Failed while writing AH VTP surface: " + filename);
  }
}

AHFinderResult WriteAHToFile(const GridData& gridData,
                             const AHFinderOptions& finderOptions,
                             const AHOutputOptions& outputOptions) {
  ValidateAHFinderOptions(finderOptions);
  ValidateAHOutputOptions(outputOptions);

  const SpectralBasis basis = BuildSpectralBasis(finderOptions.lMax);
  const CartesianPointInterpolator interpolator(
    gridData, finderOptions.interpolation);
  const SpectralNewtonResult solverResult = SolvePDEforAH(
    basis,
    interpolator,
    finderOptions.referencePoint,
    finderOptions.initialRadius,
    finderOptions.solver
  );
  AHFinderResult result = BuildAHFinderResult(
    basis, interpolator, finderOptions.referencePoint, solverResult);

  WriteAHSummary(result, outputOptions.summaryFilename);
  if (result.solver.converged && result.surfaceAvailable) {
    WriteAHSurfaceToVTP(
      result, outputOptions.surfaceFilename, outputOptions.visualization);
  }
  return result;
}

AHSearchResult WriteAHToFile(const GridData& gridData,
                             const AHSearchOptions& searchOptions,
                             const AHOutputOptions& outputOptions) {
  ValidateAHSearchOptions(searchOptions);
  ValidateAHOutputOptions(outputOptions);

  AHSearchResult result = FindApparentHorizon(gridData, searchOptions);
  WriteAHSearchSummary(result, outputOptions.summaryFilename);
  if (result.converged && result.horizon.has_value()) {
    WriteAHSurfaceToVTP(
      *result.horizon,
      outputOptions.surfaceFilename,
      outputOptions.visualization
    );
  }
  return result;
}

void WriteAHToVTK(const std::string& inputFile, const std::string& outputFile) {
  (void)inputFile;
  (void)outputFile;
  // Helper function to convert the ASKII AH file into a VTK file for visualisation.
  // – Input: ASCII file(s) contain all coordinate origins $\{s_i\}$ with a valid AH, and the 
  //          list of tuples of corresponding $ \{h_i^{(f)}(\rho, \sigma)\}$.
  // – Output: VTK file(s) for the AH surface matching the $ \{s_i, h_i^{(f)}(\rho, \sigma)\}$.
}

}  // namespace exahype2::ccz4::ahf
