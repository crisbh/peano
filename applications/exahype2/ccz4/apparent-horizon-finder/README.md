# Apparent Horizon Finder Build Notes

This folder contains the CCZ4 apparent-horizon-finder (AHF) code. It can be compiled in two different ways:

- as part of a generated CCZ4/Peano simulation executable for on-the-fly use
- as a local CMake postprocessing executable for VTK/VTU files

The source split is:

- `peano-input.*`: on-the-fly input from Peano data structures
- `file-input.*`: offline input from files, currently intended for VTK/VTU data
- `output.*`: shared output helpers for both on-the-fly and offline modes
- `interpolation.*`, `assemble-matrix.*`, `iteration.*`: shared AHF numerical pipeline pieces
- `ahf-postprocess.cpp`: local offline command-line frontend

## Local Postprocessing Build

Build from inside this folder:

```bash
cd applications/exahype2/ccz4/apparent-horizon-finder
cmake -S . -B build
cmake --build build --target ahf-postprocess
```

Run the postprocessor:

```bash
./build/ahf-postprocess --input snapshot.vtu --output horizon.vtp
```

The local postprocessing build uses `file-input.cpp`, so CMake has to find the C++ VTK package. If VTK is not found, pass either `VTK_DIR` or `CMAKE_PREFIX_PATH`:

```bash
cmake -S . -B build -DVTK_DIR=/path/to/vtk/lib/cmake/vtk-*
```

or:

```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH=/path/to/vtk
```

## Splitting a Large Snapshot (optional)

A big `.vtu` snapshot is read in one go, which needs a lot of memory. To read it
more cheaply, first split it into smaller pieces with the helper script:

```bash
python3 repartition_vtu.py snapshot.vtu pieces/snapshot.pvtu -k 8
```

This produces `pieces/snapshot.pvtu` (a summary) plus its piece files. Give that `.pvtu` to the finder instead of the original `.vtu`: it then reads one piece at a time and uses far less memory. You only need to do this once per snapshot. `ahf-io-tests` test both `vtu` and `pvtu` (under a pieces directory) approach.
Requires Python with the `vtk` module. 

## Local Postprocessing Build With Tests

Tests are optional and disabled by default. Enable them with:

```bash
cd applications/exahype2/ccz4/apparent-horizon-finder
cmake -S . -B build -DENABLE_CCZ4_AHF_TESTS=ON
cmake --build build
ctest --test-dir build --output-on-failure
```

The test executables are:

```text
ahf-io-tests
ahf-interpolation-tests
ahf-assemble-matrix-tests
ahf-iteration-tests
ahf-harmonics-tests
ahf-spectral-solver-tests
ahf-spacetimes-tests
ahf-analytic-horizon-tests
ahf-brill-lindquist-tests
```

They can also be run directly:

```bash
./build/ahf-io-tests
./build/ahf-interpolation-tests
./build/ahf-assemble-matrix-tests
./build/ahf-iteration-tests
```

## On-The-Fly CCZ4 Build

For on-the-fly use, compile AHF through the normal `ccz4.py` generator:

```bash
cd applications/exahype2/ccz4

export PYTHONPATH=../../../python
export PYTHONPATH=$PYTHONPATH:../../../applications/exahype2/ccz4

# The final --with-ahf flag enables the apparent-horizon finder.
python3 ccz4.py \
  -impl fd4-rk1-dsl \
  -s gauge \
  -maxh 0.05 \
  -minh 0.05 \
  -ps 9 \
  -plt 0.1 \
  -et 0.5 \
  -exn gauge_ahf \
  --domain_r 0.5 \
  --ReSwi 0 \
  -cfl 0.1 \
  --KOSigma 8.0 \
  --with-ahf  # enable AHF on-the-fly
```

This generates and builds:

```bash
./peano4_gauge_ahf
```

The `--with-ahf` build adds only the on-the-fly/shared sources:

```text
apparent-horizon-finder/peano-input.cpp
apparent-horizon-finder/output.cpp
apparent-horizon-finder/interpolation.cpp
apparent-horizon-finder/assemble-matrix.cpp
apparent-horizon-finder/iteration.cpp
```

It deliberately does not compile `file-input.cpp`, so the on-the-fly simulation build does not require VTK.

## Cleaning

For a CMake build directory:

```bash
cmake --build build --target clean
```

To remove the local CMake build completely:

```bash
rm -rf build
```
