#include "file-input.h"
#include "output.h"
#include "interpolation.h"
#include "assemble-matrix.h"
#include "iteration.h"

#include <cstdlib>
#include <iostream>
#include <string>

namespace {

void PrintUsage(const char* executable) {
  std::cerr
      << "Usage: " << executable << " --input <snapshot.vtu|snapshot.pvtu> --output <horizon.vtp>\n"
      << "\n"
      << "Reads a VTK snapshot produced from CCZ4 output and runs the offline\n"
      << "apparent-horizon finder pipeline. The current AHF implementation still\n"
      << "contains placeholder numerical kernels, so this frontend wires the\n"
      << "postprocessing entry point without changing the solver algorithms.\n";
}

}  // namespace

int main(int argc, char** argv) {
  std::string inputFile;
  std::string outputFile;

  for (int i = 1; i < argc; ++i) {
    const std::string argument(argv[i]);
    if (argument == "--input" && i + 1 < argc) {
      inputFile = argv[++i];
    } else if (argument == "--output" && i + 1 < argc) {
      outputFile = argv[++i];
    } else if (argument == "--help" || argument == "-h") {
      PrintUsage(argv[0]);
      return EXIT_SUCCESS;
    } else {
      std::cerr << "Unknown or incomplete argument: " << argument << "\n";
      PrintUsage(argv[0]);
      return EXIT_FAILURE;
    }
  }

  if (inputFile.empty() || outputFile.empty()) {
    PrintUsage(argv[0]);
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
