#include "io.h"

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

void ReadGridFromPeano(const PeanoGridData& peanoGridData) {
  (void)peanoGridData;
  // Loading grid data from the Peano data structure among timesteppings.
  // Must be done in parallel as each core only holds a subdomain of the grid. 
  // This approach requires binding codes to ExaGRyPE, reducing its generality.
  // – Input: Peano grid data
  // – Output: Full grid dataset. A full grid dataset giving solution $\{Q_{i,j,k}\}$ on each cell/volume, 
  //           and the corresponding Cartesian coordinates $\{x_I,y_I,z_I\}$. 
  //           They may be stored in subdomains, managed by multiple cores when the domain is large.
}

void WriteAHToFile(const std::vector<AHGuess>& initialAHGuesses,
                   const std::vector<ReferencePoint>& referencePoints,
                   const std::vector<double>& masses,
                   const std::string& outputFile) {
  (void)initialAHGuesses;
  (void)referencePoints;
  (void)masses;
  (void)outputFile;
  // Call SolvePDEforAH() to obtain the solution AH function ($ \{h^{(f)}(\rho _I, \sigma _I)\}$)
  //      for each valid reference point, and write it to an ASCII file. 
  // The ASCII file should contain the list reference point ${s_i}$, and the list of tuples of 
  //      $ \{h_i^{(f)}(\rho, \sigma)\}$ (i.e., the angular coordinates on the apparent horizon surface
  //      and their h value).
  // This function can be the external/public API for users. Users only need to call this function. 
  //      All other functions/operations are encapsulated and handled internally by this function.
  // – Input: Initial AH Guess functions ($ \{h_{i}^{(0)}(\rho _I, \sigma _I)\}$), Reference points 
  //          ($\{s_i\}$), The black hole mass associated with corresponding reference points ($\{M_i\}$).
  // – Output: ASCII file(s) contain all coordinate origins $\{s_i\}$ with a valid AH, and the 
  //           list of tuples of corresponding $ \{h^{(f)}_i(\rho, \sigma)\}$.
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
