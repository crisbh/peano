cmake_minimum_required(VERSION 3.21)

project(CCZ4ApparentHorizonFinder LANGUAGES CXX)

option(ENABLE_CCZ4_AHF_TESTS "Build CCZ4 apparent-horizon-finder test executables" OFF)

if(NOT TARGET PeanoOptions)
  add_library(PeanoOptions INTERFACE)
  target_compile_features(PeanoOptions INTERFACE cxx_std_20)
endif()

set(AHF_CORE_SOURCES
  output.cpp
  interpolation.cpp
  surface-geometry.cpp
  harmonics.cpp
  expansion.cpp
  assemble-matrix.cpp
  spectral-solver.cpp
  iteration.cpp
)

set(AHF_CORE_HEADERS
  ahf-data.h
  apparent-horizon-finder.h
  peano-input.h
  output.h
  interpolation.h
  surface-geometry.h
  tensor-algebra.h
  harmonics.h
  expansion.h
  assemble-matrix.h
  dense-linear-algebra.h
  spectral-solver.h
  iteration.h
)

add_library(ApparentHorizonFinderCore STATIC
  ${AHF_CORE_SOURCES}
  ${AHF_CORE_HEADERS}
)

target_include_directories(ApparentHorizonFinderCore
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

target_link_libraries(ApparentHorizonFinderCore
  PUBLIC
    PeanoOptions
)

target_compile_definitions(ApparentHorizonFinderCore PUBLIC IncludeApparentHorizonFinder)

find_package(VTK REQUIRED COMPONENTS
  CommonCore
  CommonDataModel
  FiltersCore
  IOXML
)

add_library(ApparentHorizonFinderPostprocess STATIC
  file-input.cpp
  file-input.h
)

target_include_directories(ApparentHorizonFinderPostprocess
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

target_link_libraries(ApparentHorizonFinderPostprocess
  PUBLIC
    ApparentHorizonFinderCore
  PRIVATE
    VTK::CommonCore
    VTK::CommonDataModel
    VTK::FiltersCore
    VTK::IOXML
)

vtk_module_autoinit(
  TARGETS ApparentHorizonFinderPostprocess
  MODULES VTK::CommonCore VTK::CommonDataModel VTK::FiltersCore VTK::IOXML
)

add_executable(ahf-postprocess ahf-postprocess.cpp)
target_link_libraries(ahf-postprocess PRIVATE ApparentHorizonFinderPostprocess)

if(ENABLE_CCZ4_AHF_TESTS)
  enable_testing()

  add_executable(ahf-io-tests test/io/IoTests.cpp)
  target_link_libraries(ahf-io-tests PRIVATE ApparentHorizonFinderPostprocess)
  add_test(NAME ahf-io-tests COMMAND ahf-io-tests)

  add_executable(ahf-data-tests test/data/AhfDataTests.cpp)
  target_link_libraries(ahf-data-tests PRIVATE ApparentHorizonFinderCore)
  add_test(NAME ahf-data-tests COMMAND ahf-data-tests)

  add_executable(ahf-interpolation-tests test/interpolation/InterpolationTests.cpp)
  target_link_libraries(ahf-interpolation-tests PRIVATE ApparentHorizonFinderCore)
  add_test(NAME ahf-interpolation-tests COMMAND ahf-interpolation-tests)

  add_executable(ahf-surface-geometry-tests test/surface-geometry/SurfaceGeometryTests.cpp)
  target_link_libraries(ahf-surface-geometry-tests PRIVATE ApparentHorizonFinderCore)
  add_test(NAME ahf-surface-geometry-tests COMMAND ahf-surface-geometry-tests)

  add_executable(ahf-ssbh-integration-tests test/integration/SSBHInterpolationTests.cpp)
  target_link_libraries(ahf-ssbh-integration-tests PRIVATE ApparentHorizonFinderPostprocess)
  add_test(NAME ahf-ssbh-integration-tests COMMAND ahf-ssbh-integration-tests)
  set_tests_properties(ahf-ssbh-integration-tests PROPERTIES LABELS integration)

  add_executable(ahf-assemble-matrix-tests test/assemble-matrix/AssembleMatrixTests.cpp)
  target_link_libraries(ahf-assemble-matrix-tests PRIVATE ApparentHorizonFinderCore)
  add_test(NAME ahf-assemble-matrix-tests COMMAND ahf-assemble-matrix-tests)

  add_executable(ahf-harmonics-tests test/harmonics/HarmonicsTests.cpp)
  target_link_libraries(ahf-harmonics-tests PRIVATE ApparentHorizonFinderCore)
  add_test(NAME ahf-harmonics-tests COMMAND ahf-harmonics-tests)

  add_executable(ahf-spectral-solver-tests test/spectral-solver/SpectralSolverTests.cpp)
  target_link_libraries(ahf-spectral-solver-tests PRIVATE ApparentHorizonFinderCore)
  add_test(NAME ahf-spectral-solver-tests COMMAND ahf-spectral-solver-tests)

  add_executable(ahf-iteration-tests test/iteration/IterationTests.cpp)
  target_link_libraries(ahf-iteration-tests PRIVATE ApparentHorizonFinderCore)
  add_test(NAME ahf-iteration-tests COMMAND ahf-iteration-tests)
endif()
