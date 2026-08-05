#include <cassert>
#include <iostream>
#include <stdexcept>

#include "../../ahf-data.h"
#include "../../file-input.h"

namespace exahype2::ccz4::ahf::tests::data {

void testReferencePointStruct() {
  const ReferencePoint referencePoint{1.0, 0.25, 0.5};

  assert(referencePoint.x == 1.0);
  assert(referencePoint.y == 0.25);
  assert(referencePoint.z == 0.5);

  std::cout << "ReferencePoint stores x = " << referencePoint.x << ", y = " << referencePoint.y << ", and z = " << referencePoint.z << std::endl;
}

void testAHGuessStruct() {
  const std::vector<AHGuess> initialAHGuesses = {
    {{{1.0, 0.0, 0.0}, {1.1, 0.5, 1.0}}},
    {{{2.0, 0.0, 0.0}, {2.2, 0.5, 1.0}}}
  };

  assert(initialAHGuesses.size() == 2);
  assert(initialAHGuesses[0].points.size() == 2);
  assert(initialAHGuesses[0].points[1].r == 1.1);
  assert(initialAHGuesses[0].points[1].rho == 0.5);
  assert(initialAHGuesses[0].points[1].sigma == 1.0);
  assert(initialAHGuesses[1].points[0].r == 2.0);

  std::cout << "AHGuess stores lists of (r, rho, sigma) points, and initialAHGuesses stores multiple guesses" << std::endl;
}

void testAngularSampleStruct() {
  const std::vector<AngularSample> angularSamples{
    {.rho = 0.25, .sigma = 0.5},
    {.rho = 0.75, .sigma = 1.0}
  };

  assert(angularSamples.size() == 2);
  assert(angularSamples[0].rho == 0.25);
  assert(angularSamples[0].sigma == 0.5);
  assert(angularSamples[1].rho == 0.75);
  assert(angularSamples[1].sigma == 1.0);

  std::cout << "AngularSample stores named (rho, sigma) coordinates" << std::endl;
}

void testGridDataStruct() {
  GridData gridData;
  gridData.coordinates = {
    {0.0, 0.1, 0.2},
    {1.0, 1.1, 1.2}
  };
  gridData.values.resize(gridData.coordinates.size() * GridData::numberOfVariables);

  gridData.value({0.0, 0.1, 0.2}, 0) = 10.0;
  gridData.value({0.0, 0.1, 0.2}, 1) = 11.0;
  gridData.value({0.0, 0.1, 0.2}, 2) = 12.0;
  gridData.value({1.0, 1.1, 1.2}, 0) = 20.0;
  gridData.value({1.0, 1.1, 1.2}, 1) = 21.0;
  gridData.value({1.0, 1.1, 1.2}, 2) = 22.0;

  assert(gridData.numberOfPoints() == 2);
  assert(gridData.numberOfVariables == 14);
  assert(gridData.isConsistent());
  assert(gridData.coordinates[1].x == 1.0);
  assert(gridData.coordinates[1].y == 1.1);
  assert(gridData.coordinates[1].z == 1.2);
  assert(gridData.value({0.0, 0.1, 0.2}, 2) == 12.0);
  assert(gridData.value({1.0, 1.1, 1.2}, 1) == 21.0);
  assert(gridData.value({1.0, 0.5 + 0.6, 1.2}, 1) == 21.0);

  GridData incompleteRow = gridData;
  incompleteRow.values.pop_back();
  assert(!incompleteRow.isConsistent());

  const GridData& constGridData = gridData;
  assert(constGridData.value({1.0, 1.1, 1.2}, 2) == 22.0);

  bool missingCoordinateRejected = false;
  try {
    (void)gridData.value({2.0, 2.1, 2.2}, 0);
  } catch (const std::out_of_range&) {
    missingCoordinateRejected = true;
  }
  assert(missingCoordinateRejected);

  bool invalidVariableRejected = false;
  try {
    (void)gridData.value({0.0, 0.1, 0.2}, GridData::numberOfVariables);
  } catch (const std::out_of_range&) {
    invalidVariableRejected = true;
  }
  assert(invalidVariableRejected);

  std::cout << "GridData stores coordinates and point-major scalar values" << std::endl;
}

void testCompactCellAccessor() {
  GridData gridData;
  gridData.coordinates = {{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}};
  gridData.values.resize(gridData.coordinates.size() * GridData::numberOfVariables);
  for (std::size_t column = 0; column < GridData::numberOfVariables; ++column) {
    gridData.values[column] = 100.0 + static_cast<double>(column);
    gridData.values[GridData::numberOfVariables + column] =
      200.0 + static_cast<double>(column);
  }

  const GridData& constGridData = gridData;
  const CompactCCZ4CellView first = constGridData.cell(0);
  const CompactCCZ4CellView second = constGridData.cell(1);

  // The view aliases GridData storage instead of copying the fourteen values.
  assert(first.row().data() == gridData.values.data());
  assert(second.row().data() == gridData.values.data() + GridData::numberOfVariables);
  assert(first.row().size() == NumberOfCompactCCZ4Fields);

  for (std::size_t component = 0; component < NumberOfSymmetricTensorComponents; ++component) {
    assert(first.conformalMetricComponent(component) == 100.0 + component);
    assert(first.conformalExtrinsicCurvatureComponent(component) == 106.0 + component);
    assert(second.conformalMetricComponent(component) == 200.0 + component);
    assert(second.conformalExtrinsicCurvatureComponent(component) == 206.0 + component);
  }

  assert(first.conformalMetricComponent(SymmetricTensorComponent::XX) == 100.0);
  assert(first.conformalMetricComponent(SymmetricTensorComponent::ZZ) == 105.0);
  assert(first.conformalExtrinsicCurvatureComponent(SymmetricTensorComponent::XX) == 106.0);
  assert(first.conformalExtrinsicCurvatureComponent(SymmetricTensorComponent::ZZ) == 111.0);
  assert(first.traceK() == 112.0);
  assert(first.phi() == 113.0);
  assert(second.traceK() == 212.0);
  assert(second.phi() == 213.0);

  // This specifically fails if the source VTU index 53 is used on a 14-value row.
  assert(first.traceK() == first.row()[static_cast<std::size_t>(CompactField::TraceK)]);
  assert(static_cast<std::size_t>(CompactField::TraceK) == 12);

  bool invalidCellRejected = false;
  try {
    (void)constGridData.cell(2);
  } catch (const std::out_of_range& error) {
    invalidCellRejected = true;
    assert(std::string(error.what()).find("cell index 2") != std::string::npos);
  }
  assert(invalidCellRejected);

  bool invalidMetricComponentRejected = false;
  try {
    (void)first.conformalMetricComponent(NumberOfSymmetricTensorComponents);
  } catch (const std::out_of_range& error) {
    invalidMetricComponentRejected = true;
    assert(std::string(error.what()).find("conformal metric component index 6")
           != std::string::npos);
  }
  assert(invalidMetricComponentRejected);

  bool invalidATildeComponentRejected = false;
  try {
    (void)first.conformalExtrinsicCurvatureComponent(NumberOfSymmetricTensorComponents);
  } catch (const std::out_of_range& error) {
    invalidATildeComponentRejected = true;
    assert(std::string(error.what()).find(
      "conformal extrinsic curvature component index 6") != std::string::npos);
  }
  assert(invalidATildeComponentRejected);

  GridData malformed = gridData;
  malformed.values.pop_back();
  bool malformedRejected = false;
  try {
    (void)malformed.cell(0);
  } catch (const std::logic_error& error) {
    malformedRejected = true;
    assert(std::string(error.what()).find("malformed compact data") != std::string::npos);
  }
  assert(malformedRejected);

  std::cout << "Compact cell accessor provides checked zero-copy named field access"
            << std::endl;
}

void testUnpackSymmetricTensor() {
  const std::array<double, NumberOfSymmetricTensorComponents> compact{
    1.0, 2.0, 3.0, 4.0, 5.0, 6.0
  };

  const Tensor3 tensor = UnpackSymmetricTensor(compact);
  const Tensor3 expected{{
    {{1.0, 2.0, 3.0}},
    {{2.0, 4.0, 5.0}},
    {{3.0, 5.0, 6.0}}
  }};

  assert(tensor == expected);
  for (std::size_t row = 0; row < 3; ++row) {
    for (std::size_t column = 0; column < 3; ++column) {
      assert(tensor[row][column] == tensor[column][row]);
    }
  }

  std::cout << "Symmetric tensor unpacking preserves component order and symmetry"
            << std::endl;
}

void testSourceVariableMappingStruct() {
  const SourceVariableMapping mapping{
    .conformalMetric = {0, 1, 2, 3, 4, 5},
    .conformalExtrinsicCurvature = {6, 7, 8, 9, 10, 11},
    .traceK = 53,
    .conformalFactor = 54
  };

  assert(mapping.conformalMetric[SourceVariableMapping::Index{0}] == 0);
  assert(mapping.conformalExtrinsicCurvature[SourceVariableMapping::Index{5}] == 11);
  assert(mapping.traceK == 53);
  assert(mapping.conformalFactor == 54);
  assert(mapping.isComplete());
  assert(mapping.isValid(59));
  assert(!mapping.isValid(54));

  const SourceVariableMapping incompleteMapping;
  assert(!incompleteMapping.isComplete());
  assert(!incompleteMapping.isValid(59));

  std::cout << "SourceVariableMapping stores validated source VTU component indices" << std::endl;
}

void testInterpolatedShellStruct() {
  InterpolatedShell shell{
    .innerRadius = 0.8,
    .outerRadius = 1.2,
    .points = {
      {.r = 0.8, .rho = 0.25, .sigma = 0.5},
      {.r = 1.0, .rho = 0.25, .sigma = 0.5},
      {.r = 1.2, .rho = 0.25, .sigma = 0.5}
    }
  };

  shell.points[1].geometry.T[0][1] = 2.0;
  shell.points[1].geometry.coordinateGradient[0][2] = 3.0;
  shell.points[1].geometry.coordinateHessian[1][2][0] = 4.0;
  shell.points[1].geometry.spatialMetric[2][2] = 5.0;
  shell.points[1].geometry.spatialMetricDerivative[0][1][2] = 6.0;
  shell.points[1].geometry.extrinsicCurvature[1][1] = 7.0;
  shell.points[1].geometry.traceK = 8.0;

  shell.value({1.0, 0.25, 0.5}).geometry.traceK = 9.0;

  assert(shell.numberOfPoints() == 3);
  assert(shell.points[1].r == 1.0);
  assert(shell.points[1].rho == 0.25);
  assert(shell.points[1].sigma == 0.5);
  assert(shell.points[1].geometry.T[0][1] == 2.0);
  assert(shell.points[1].geometry.coordinateGradient[0][2] == 3.0);
  assert(shell.points[1].geometry.coordinateHessian[1][2][0] == 4.0);
  assert(shell.points[1].geometry.spatialMetric[2][2] == 5.0);
  assert(shell.points[1].geometry.spatialMetricDerivative[0][1][2] == 6.0);
  assert(shell.points[1].geometry.extrinsicCurvature[1][1] == 7.0);
  assert(shell.points[1].geometry.traceK == 9.0);
  const InterpolatedShell& constShell = shell;
  assert(constShell.value({0.5 + 0.5, 0.1 + 0.15, 0.5}).geometry.traceK == 9.0);
  assert(shell.isConsistent());

  bool missingShellCoordinateRejected = false;
  try {
    (void)shell.value({1.0, 0.75, 0.5});
  } catch (const std::out_of_range&) {
    missingShellCoordinateRejected = true;
  }
  assert(missingShellCoordinateRejected);

  InterpolatedShell invalidShell = shell;
  invalidShell.points[0].r = 0.7;
  assert(!invalidShell.isConsistent());

  const InterpolatedShell emptyShell;
  assert(!emptyShell.isConsistent());

  std::cout << "InterpolatedShell stores pointwise (r, rho, sigma) data" << std::endl;
}

void testInterpolatedSurfaceStruct() {
  const AHGuess ahGuess{
    .points = {
      {.r = 1.0, .rho = 0.25, .sigma = 0.5},
      {.r = 1.1, .rho = 0.75, .sigma = 1.0}
    }
  };

  InterpolatedSurface surface{
    .points = {
      {.r = 1.0, .rho = 0.25, .sigma = 0.5},
      {.r = 1.1, .rho = 0.75, .sigma = 1.0}
    }
  };

  surface.points[0].surfaceNormal[2] = 2.0;
  surface.points[0].surfaceNormalDerivative[1][2] = 3.0;
  surface.points[0].geometry.spatialMetric[0][0] = 4.0;

  surface.value({1.0, 0.25, 0.5}).geometry.traceK = 5.0;

  assert(surface.numberOfPoints() == 2);
  assert(surface.points[0].surfaceNormal[2] == 2.0);
  assert(surface.points[0].surfaceNormalDerivative[1][2] == 3.0);
  assert(surface.points[0].geometry.spatialMetric[0][0] == 4.0);
  const InterpolatedSurface& constSurface = surface;
  assert(constSurface.value({0.5 + 0.5, 0.1 + 0.15, 0.5}).geometry.traceK == 5.0);
  assert(surface.isConsistent());
  assert(surface.matches(ahGuess));

  bool missingSurfaceCoordinateRejected = false;
  try {
    (void)surface.value({1.0, 0.75, 0.5});
  } catch (const std::out_of_range&) {
    missingSurfaceCoordinateRejected = true;
  }
  assert(missingSurfaceCoordinateRejected);

  AHGuess differentGuess = ahGuess;
  differentGuess.points[1].r = 1.2;
  assert(!surface.matches(differentGuess));

  InterpolatedSurface invalidSurface = surface;
  invalidSurface.points[0].r = -1.0;
  assert(!invalidSurface.isConsistent());

  const InterpolatedSurface emptySurface;
  assert(!emptySurface.isConsistent());

  std::cout << "InterpolatedSurface stores one point per AH guess point" << std::endl;
}

}  // namespace exahype2::ccz4::ahf::tests::data

int main() {
  using namespace exahype2::ccz4::ahf::tests::data;

  std::cout << "Apparent Horizon Finder data tests begin..." << std::endl;

  testReferencePointStruct();
  testAHGuessStruct();
  testAngularSampleStruct();
  testGridDataStruct();
  testCompactCellAccessor();
  testUnpackSymmetricTensor();
  testSourceVariableMappingStruct();
  testInterpolatedShellStruct();
  testInterpolatedSurfaceStruct();

  std::cout << "All data tests completed." << std::endl;
  return 0;
}
