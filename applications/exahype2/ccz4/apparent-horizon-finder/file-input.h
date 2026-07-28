#ifndef APPARENT_HORIZON_FINDER_FILE_INPUT_H_
#define APPARENT_HORIZON_FINDER_FILE_INPUT_H_

#include "ahf-data.h"

#include <string>

namespace exahype2::ccz4::ahf {

/**
 * Load a grid dataset from VTK XML unstructured-grid data, keeping only the
 * variables selected by @p variableMapping.
 *
 * The input may be either a single-piece `.vtu` (read in one shot) or a
 * multi-piece `.pvtu` series (dispatched by the `.pvtu` extension). The `.pvtu`
 * path is streamed one piece at a time with a fresh reader per piece, so peak
 * memory is bounded by one piece plus the accumulating GridData rather than by
 * the whole grid. Use repartition_vtu.py to turn a single-piece snapshot into a
 * `.pvtu` series for low-memory reads.
 *
 * The CCZ4 finite-volume output stores its unknowns as a single multi-component
 * cell-data array (one finite volume per cell). Rather than copying all
 * components, this reader stores only the source components referenced by the
 * mapping. This keeps GridData small: a full CCZ4 mapping selects 14 of the 64
 * stored components, so GridData::values shrinks accordingly.
 *
 * Output layout. GridData::coordinates[i] is the centre of cell i. GridData is
 * compact: it has one variable column per *set* field of the mapping, in this
 * fixed order, skipping any field left at VariableMapping::NotStored:
 *   conformalMetric              xx, xy, xz, yy, yz, zz   (columns 0..5)
 *   conformalExtrinsicCurvature  xx, xy, xz, yy, yz, zz   (columns 6..11)
 *   traceK                                                (column 12)
 *   conformalFactor                                       (column 13)
 * GridData::numberOfVariables equals the number of set fields. The source
 * component index that feeds column k is therefore the k-th set field of the
 * mapping in the order above; the mapping's own indices refer to the file, not
 * to the compact GridData layout.
 *
 * @param inputFile       path to a .vtu file.
 * @param gridData        output; any previous content is overwritten.
 * @param variableMapping selects which source components to read; every set
 *                        index must be a valid component of the cell-data array.
 * @throws std::runtime_error if the file cannot be read, holds no cell data,
 *         the mapping selects no field, or a selected index is out of range.
 */
void ReadGridFromFile(const std::string& inputFile,
                      GridData& gridData,
                      const VariableMapping& variableMapping);

}  // namespace exahype2::ccz4::ahf

#endif  // APPARENT_HORIZON_FINDER_FILE_INPUT_H_
