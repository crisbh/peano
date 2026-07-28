#include <array>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

#include "../../ahf-data.h"
#include "../../file-input.h"
#include "../../output.h"
#include "../../peano-input.h"

namespace exahype2::ccz4::ahf::tests::io {

/**
 * Resolve a hard-coded sample file name relative to this source file so the test
 * works irrespective of the working directory ctest is launched from.
 */
std::string samplePath(const std::string& relative) {
  const std::string thisFile = __FILE__;                       // .../test/io/IoTests.cpp
  const std::string thisDir = thisFile.substr(0, thisFile.find_last_of("/\\"));
  return thisDir + "/../../" + relative;
}

bool fileExists(const std::string& path) {
  std::ifstream probe(path, std::ios::binary);
  return probe.good();
}

/**
 * Read the CCZ4 sample through ReadGridFromFile() and assert it loaded
 * correctly. The same data and ground truth hold whether it is read from the
 * single-piece .vtu or the repartitioned .pvtu series, so this checker is shared
 * between both modes.
 */
void checkSampleGrid(const std::string& inputFile) {
  // The CCZ4 sample stores 64 cell components; this is the sample's actual
  // mapping. It selects 14 of them, which become GridData's compact columns:
  // conformal metric 0-5 -> columns 0-5, conformal extrinsic curvature 6-11 ->
  // columns 6-11, traceK (source 53) -> column 12, conformal factor (source 54)
  // -> column 13.
  VariableMapping variableMapping;
  variableMapping.conformalMetric = {0, 1, 2, 3, 4, 5};
  variableMapping.conformalExtrinsicCurvature = {6, 7, 8, 9, 10, 11};
  variableMapping.traceK = 53;
  variableMapping.conformalFactor = 54;

  GridData gridData;
  ReadGridFromFile(inputFile, gridData, variableMapping);

  // sample_data_S: 3088448 cells; the mapping keeps 14 variables.
  constexpr std::size_t ExpectedNumberOfCells = 3088448;
  constexpr std::size_t ExpectedNumberOfVariables = 14;

  assert(gridData.numberOfPoints() == ExpectedNumberOfCells);
  assert(gridData.numberOfVariables == ExpectedNumberOfVariables);
  assert(gridData.isConsistent());
  assert(gridData.values.size() == ExpectedNumberOfCells * ExpectedNumberOfVariables);

  // Coordinates (cell centroids) and stored unknowns must be finite. Index the
  // value rows directly here; GridData::value() does a linear coordinate search,
  // which is too costly to call per variable over millions of cells.
  const GridCoordinate& firstCoordinate = gridData.coordinates.front();
  assert(std::isfinite(firstCoordinate.x));
  assert(std::isfinite(firstCoordinate.y));
  assert(std::isfinite(firstCoordinate.z));

  const GridCoordinate& lastCoordinate = gridData.coordinates.back();
  assert(std::isfinite(lastCoordinate.x));
  assert(std::isfinite(lastCoordinate.y));
  assert(std::isfinite(lastCoordinate.z));

  const std::size_t lastRow = (gridData.numberOfPoints() - 1) * gridData.numberOfVariables;
  for (std::size_t variable = 0; variable < gridData.numberOfVariables; ++variable) {
    assert(std::isfinite(gridData.values[variable]));            // first cell's row
    assert(std::isfinite(gridData.values[lastRow + variable]));  // last cell's row
  }

  // Exercise the coordinate-based enquiry API once (it is a linear search).
  assert(std::isfinite(gridData.value(firstCoordinate, 0)));

  // Ground truth for three known cells (centroids nearest to {1,0,0}, {0,3,3},
  // {5,5,6}), read directly from the sample. Setting each cell's coordinate must
  // return its stored unknowns, remapped into the compact layout: the checked
  // columns 0, 6, 12, 13 hold source components 0 (metric xx), 6 (extrinsic
  // curvature xx), 53 (traceK) and 54 (conformal factor). Only the conformal
  // factor differs across the three cells; the coordinate lookup pins down cell
  // identity. Repartitioning preserves cells, so these hold in both modes.
  struct KnownComponent {
    GridCoordinate coordinate;
    std::size_t column;
    double expected;
  };
  constexpr GridCoordinate cell1{
    0.98148149251937866, 0.018518491892685951, 0.018518491892685951};
  constexpr GridCoordinate cell2{
    -0.055555548518896103, 2.9444445371627808, 2.9444445371627808};
  constexpr GridCoordinate cell3{
    4.8333337306976318, 4.8333337306976318, 5.8333334922790527};

  const std::array<KnownComponent, 12> knownComponents{{
    {cell1, 0, 1.0},                   // source 0  (metric xx)
    {cell1, 6, 0.0},                   // source 6  (extrinsic curvature xx)
    {cell1, 12, 0.0},                  // source 53 (traceK)
    {cell1, 13, 0.43901168099999999},  // source 54 (conformal factor)
    {cell2, 0, 1.0},
    {cell2, 6, 0.0},
    {cell2, 12, 0.0},
    {cell2, 13, 0.79710270500000002},
    {cell3, 0, 1.0},
    {cell3, 6, 0.0},
    {cell3, 12, 0.0},
    {cell3, 13, 0.89736081099999998},
  }};

  for (const KnownComponent& known : knownComponents) {
    const double actual = gridData.value(known.coordinate, known.column);
    const double tolerance = 1.0e-10 * (1.0 + std::abs(known.expected));
    assert(std::abs(actual - known.expected) <= tolerance);
  }

  std::cout << "ReadGridFromFile() loaded " << gridData.numberOfPoints()
            << " cells with " << gridData.numberOfVariables
            << " variables each from " << inputFile << std::endl;
}

void testReadGridFromFile() {
  const std::string vtuFile = samplePath("sample_data_S.vtu");
  const std::string pvtuFile = samplePath("pieces/sample_data_S.pvtu");

  // Mode 1: single-piece .vtu, read in one shot.
  if (fileExists(vtuFile)) {
    checkSampleGrid(vtuFile);
  } else {
    std::cout << "Skipping .vtu mode: sample not found at " << vtuFile << std::endl;
  }

  // Mode 2: repartitioned .pvtu series, streamed piece by piece. Generate with:
  //   python3 repartition_vtu.py sample_data_S.vtu pieces/sample_data_S.pvtu
  if (fileExists(pvtuFile)) {
    checkSampleGrid(pvtuFile);
  } else {
    std::cout << "Skipping .pvtu mode: series not found at " << pvtuFile << std::endl;
  }
}

void testReadGridFromPeano() {
  std::cout << "Placeholder test for peano-input ReadGridFromPeano()" << std::endl;
}

void testWriteAHToFile() {
  std::cout << "Placeholder test for output WriteAHToFile()" << std::endl;
}

void testWriteAHToVTK() {
  std::cout << "Placeholder test for output WriteAHToVTK()" << std::endl;
}

}  // namespace exahype2::ccz4::ahf::tests::io

int main() {
  using namespace exahype2::ccz4::ahf::tests::io;

  std::cout << "Apparent Horizon Finder I/O tests begin..." << std::endl;

  testReadGridFromFile();
  testReadGridFromPeano();
  testWriteAHToFile();
  testWriteAHToVTK();

  std::cout << "All I/O tests completed." << std::endl;
  return 0;
}
