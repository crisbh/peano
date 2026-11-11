# Apparent Horizon Finder Build Notes

This folder contains ExaGRyPE's apparent-horizon-finder (AHF) code. At the moment it is a standalone CMake postprocessing tool for VTK/VTU snapshots, not part of the ExaGRyPE library; on-the-fly use inside a simulation is planned (see below).

The source split is:

- `peano-input.*`: on-the-fly input from Peano data structures
- `file-input.*`: offline input from files, currently intended for VTK/VTU data
- `output.*`: shared output helpers for both on-the-fly and offline modes
- `interpolation.*`, `assemble-matrix.*`, `iteration.*`: shared AHF numerical pipeline pieces
- `ahf-postprocess.cpp`: local offline command-line frontend

## Cartesian interpolation modes

`CartesianToSurfaceOptions::method` keeps the existing interpolation choices
and also provides an opt-in smooth-support choice:

- `WeightedPolynomial` selects exactly `numberOfNeighbours` nearest cell
  centres and applies the existing Gaussian moving-least-squares weights.
- `RadiusWeightedPolynomial` selects every cell centre at distance strictly
  less than `supportRadius` and applies the compact Wendland C2 weight
  `(1-q)^4(4q+1)` for `q=distance/supportRadius < 1`. The weight is zero at and
  outside the support boundary.

The radius-supported mode is useful when the interpolated target moves during
a Newton iteration: support membership can change while a cell's contribution
vanishes smoothly at the boundary. It requires a finite positive
`supportRadius` and enough geometrically independent cells inside that radius
for the requested `polynomialDegree`. The fixed-count method remains unchanged
and is still the default weighted-polynomial behaviour.

## Local Postprocessing Build

Build from inside this folder:

```bash
cd src/exagrype/ahf
cmake -S . -B build
cmake --build build --target ahf-postprocess
```

It can also be built as part of Peano's CMake build, by adding `-DENABLE_EXAGRYPE=ON -DENABLE_CCZ4_AHF_POSTPROCESSOR=ON` (and `-DENABLE_CCZ4_AHF_TESTS=ON` for the tests) to the usual ExaGRyPE configuration.

Run the postprocessor:

```bash
./build/ahf-postprocess --input snapshot.vtu --output horizon.vtp
```

The local postprocessing build uses `file-input.cpp`, so CMake has to find the C++ VTK package. Without VTK, the core and the tests that do not read `.vtu` files are still built, but `ahf-postprocess` and the `ahf-io`, `ahf-multi-horizon-dataset` and `ahf-sample-t0-user-workflow` tests are skipped. If VTK is installed but not found, pass either `VTK_DIR` or `CMAKE_PREFIX_PATH`:

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
cd src/exagrype/ahf
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

## On-The-Fly Use (planned)

On-the-fly use inside an ExaGRyPE simulation is not available yet. The `--with-ahf` flag of the old `applications/exahype2/ccz4/ccz4.py` script has been removed, as ExaGRyPE runs are now built through `exagrype.api` (see `benchmarks/exagrype`). The plan is to compile the core into the ExaGRyPE library and add an `exagrype.api` action set that switches the finder on with `-DIncludeApparentHorizonFinder`. The Peano-side input, `peano-input.*` (`ReadGridFromPeano`), is still a stub. The on-the-fly build will not compile `file-input.cpp`, so it will not require VTK.

## Cleaning

For a CMake build directory:

```bash
cmake --build build --target clean
```

To remove the local CMake build completely:

```bash
rm -rf build
```
