#include "file-input.h"

namespace exahype2::ccz4::ahf {

void ReadGridFromFile(const std::string& inputFile) {
  (void)inputFile;
  // Loading grid data from a standard VTK/HDF5 file. 
  // This approach has the maximal compatibility, with the price of extra I/O to disk (thus speed). 
  // Parallelism can be used for performance (and to avoid exceeding the memory threshold). 
  //             Ref: https://docs.vtk.org/en/latest/design_documents/smptools.html.
  // This approach can be applied either on-the-fly or during post-processing.
  // – Input: VTK, HDF5...
  // – Output: Full grid dataset. A full grid dataset giving solution $\{Q_{i,j,k}\}$ on each cell/volume, 
  //           and the corresponding Cartesian coordinates $\{x_I,y_I,z_I\}$. 
  //           They may be stored in subdomains, managed by multiple cores when the domain is large.
}

}  // namespace exahype2::ccz4::ahf
