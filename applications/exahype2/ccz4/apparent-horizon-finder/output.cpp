#include "output.h"

namespace exahype2::ccz4::ahf {

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
