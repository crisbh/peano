#include "file-input.h"

#include <vtkCellCenters.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLPUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridReader.h>

#include <stdexcept>
#include <string>
#include <vector>

namespace exahype2::ccz4::ahf {

namespace {

/**
 * check if @p path has the @extension suffix.
 */
bool hasExtension(const std::string& path, const std::string& extension) {
  return path.size() >= extension.size()
      && path.compare(path.size() - extension.size(), extension.size(), extension) == 0;
}

/**
 * Collect the source component indices selected by the mapping, in the fixed
 * field order documented in file-input.h. Fields left at NotStored are skipped;
 * each surviving index becomes one GridData column.
 */
std::vector<std::size_t> selectedComponents(const VariableMapping& variableMapping) {
  std::vector<std::size_t> selected;
  selected.reserve(14);

  const auto append = [&selected](VariableMapping::Index index) {
    if (index != VariableMapping::NotStored) {
      selected.push_back(index);
    }
  };

  for (VariableMapping::Index index : variableMapping.conformalMetric) {
    append(index);
  }
  for (VariableMapping::Index index : variableMapping.conformalExtrinsicCurvature) {
    append(index);
  }
  append(variableMapping.traceK);
  append(variableMapping.conformalFactor);

  return selected;
}

/**
 * Append every cell of @p grid to @p gridData: the cell centre as a coordinate
 * and the selected source components as one compact value row. gridData must
 * already have its number of variables set to selected.size().
 */
void appendGridCells(vtkUnstructuredGrid* grid,
                     const std::vector<std::size_t>& selected,
                     const std::string& context,
                     GridData& gridData) {
  const vtkIdType numberOfCells = (grid != nullptr) ? grid->GetNumberOfCells() : 0;
  if (numberOfCells == 0) {
    return;
  }

  // The CCZ4 solution is stored as cell data: one multi-component tuple per
  // finite volume. Prefer the array flagged as the active scalars, falling back
  // to the first cell-data array.
  vtkCellData* cellData = grid->GetCellData();
  vtkDataArray* unknowns = (cellData != nullptr) ? cellData->GetScalars() : nullptr;
  if (unknowns == nullptr && cellData != nullptr && cellData->GetNumberOfArrays() > 0) {
    unknowns = cellData->GetArray(0);
  }
  if (unknowns == nullptr) {
    throw std::runtime_error("ReadGridFromFile: no cell-data array found in " + context);
  }

  const int numberOfComponents = unknowns->GetNumberOfComponents();
  if (unknowns->GetNumberOfTuples() != numberOfCells) {
    throw std::runtime_error(
      "ReadGridFromFile: cell-data tuple count does not match the number of cells in "
      + context);
  }
  for (std::size_t sourceIndex : selected) {
    if (sourceIndex >= static_cast<std::size_t>(numberOfComponents)) {
      throw std::runtime_error(
        "ReadGridFromFile: the variable mapping selects a component index that the "
        "cell-data array does not have in " + context);
    }
  }

  // The sample coordinate of each finite volume is its cell centre. vtkCellCenters
  // evaluates the cell's parametric centre, which is correct for every cell type
  // (unlike a plain corner average, which only matches for linear box/simplex
  // cells). Output point i corresponds to cell i.
  vtkNew<vtkCellCenters> cellCenters;
  cellCenters->SetInputData(grid);
  cellCenters->Update();
  vtkPoints* centers =
    (cellCenters->GetOutput() != nullptr) ? cellCenters->GetOutput()->GetPoints() : nullptr;
  if (centers == nullptr || centers->GetNumberOfPoints() != numberOfCells) {
    throw std::runtime_error(
      "ReadGridFromFile: failed to compute a cell centre for every cell in " + context);
  }

  const std::size_t base = gridData.coordinates.size();
  const std::size_t numberOfVariables = gridData.numberOfVariables;
  gridData.coordinates.resize(base + static_cast<std::size_t>(numberOfCells));
  gridData.values.resize((base + static_cast<std::size_t>(numberOfCells)) * numberOfVariables);

  // Read one full tuple per cell and keep only the selected components, so the
  // large all-component array is never duplicated into GridData.
  std::vector<double> sourceTuple(static_cast<std::size_t>(numberOfComponents));
  for (vtkIdType cellIndex = 0; cellIndex < numberOfCells; ++cellIndex) {
    const std::size_t index = base + static_cast<std::size_t>(cellIndex);

    double centre[3];
    centers->GetPoint(cellIndex, centre);
    gridData.coordinates[index] = GridCoordinate{centre[0], centre[1], centre[2]};

    unknowns->GetTuple(cellIndex, sourceTuple.data());
    double* row = &gridData.values[index * numberOfVariables];
    for (std::size_t column = 0; column < selected.size(); ++column) {
      row[column] = sourceTuple[selected[column]];
    }
  }
}

/**
 * Read a single-piece .vtu in one shot.
 */
void readWholeVtu(const std::string& inputFile,
                  const std::vector<std::size_t>& selected,
                  GridData& gridData) {
  vtkNew<vtkXMLUnstructuredGridReader> reader;
  reader->SetFileName(inputFile.c_str());
  if (reader->CanReadFile(inputFile.c_str()) == 0) {
    throw std::runtime_error("ReadGridFromFile: cannot read file " + inputFile);
  }
  reader->Update();

  vtkUnstructuredGrid* grid = reader->GetOutput();
  if (grid == nullptr || grid->GetNumberOfCells() == 0) {
    throw std::runtime_error(
      "ReadGridFromFile: no unstructured-grid cells found in " + inputFile);
  }

  gridData.coordinates.clear();
  gridData.values.clear();
  gridData.setNumberOfVariables(selected.size());
  gridData.coordinates.reserve(static_cast<std::size_t>(grid->GetNumberOfCells()));
  gridData.values.reserve(static_cast<std::size_t>(grid->GetNumberOfCells()) * selected.size());

  appendGridCells(grid, selected, inputFile, gridData);
}

/**
 * Stream a multi-piece .pvtu one piece at a time. A fresh reader is used per
 * piece so that only one piece is resident at once: peak memory is bounded by
 * one piece plus the accumulating (compact) GridData, not by the whole grid.
 */
void readPvtuStreaming(const std::string& inputFile,
                       const std::vector<std::size_t>& selected,
                       GridData& gridData) {
  int numberOfPieces = 0;
  {
    vtkNew<vtkXMLPUnstructuredGridReader> summaryReader;
    summaryReader->SetFileName(inputFile.c_str());
    if (summaryReader->CanReadFile(inputFile.c_str()) == 0) {
      throw std::runtime_error("ReadGridFromFile: cannot read file " + inputFile);
    }
    summaryReader->UpdateInformation();
    numberOfPieces = summaryReader->GetNumberOfPieces();
  }
  if (numberOfPieces < 1) {
    throw std::runtime_error("ReadGridFromFile: " + inputFile + " declares no pieces");
  }

  gridData.coordinates.clear();
  gridData.values.clear();
  gridData.setNumberOfVariables(selected.size());

  for (int piece = 0; piece < numberOfPieces; ++piece) {
    // Fresh reader per piece: reusing one reader across UpdatePiece() calls
    // caches internal buffers and climbs back to the whole-grid footprint.
    vtkNew<vtkXMLPUnstructuredGridReader> reader;
    reader->SetFileName(inputFile.c_str());
    reader->UpdatePiece(piece, numberOfPieces, 0);
    vtkUnstructuredGrid* grid = reader->GetOutput();

    if (piece == 0 && grid != nullptr) {
      // Estimate the total from the first piece and reserve, so the growing
      // GridData is not reallocated (which would transiently double its memory).
      // The "+ numberOfPieces" covers the uneven remainder of an equal split.
      const std::size_t estimate =
        static_cast<std::size_t>(grid->GetNumberOfCells()) * numberOfPieces
        + static_cast<std::size_t>(numberOfPieces);
      gridData.coordinates.reserve(estimate);
      gridData.values.reserve(estimate * gridData.numberOfVariables);
    }

    appendGridCells(grid, selected, inputFile, gridData);
  }

  if (gridData.coordinates.empty()) {
    throw std::runtime_error("ReadGridFromFile: no cells found across the pieces of " + inputFile);
  }
}

}  // namespace

void ReadGridFromFile(const std::string& inputFile,
                      GridData& gridData,
                      const VariableMapping& variableMapping) {
  const std::vector<std::size_t> selected = selectedComponents(variableMapping);
  if (selected.empty()) {
    throw std::runtime_error(
      "ReadGridFromFile: the variable mapping selects no variable to read");
  }

  if (hasExtension(inputFile, ".pvtu")) {
    readPvtuStreaming(inputFile, selected, gridData);
  } else {
    readWholeVtu(inputFile, selected, gridData);
  }
}

}  // namespace exahype2::ccz4::ahf
