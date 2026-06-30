cmake_minimum_required(VERSION 3.21)

project(CCZ4ApparentHorizonFinder LANGUAGES CXX)

if(NOT TARGET PeanoOptions)
  add_library(PeanoOptions INTERFACE)
  target_compile_features(PeanoOptions INTERFACE cxx_std_20)
endif()

set(AHF_CORE_SOURCES
  output.cpp
  interpolation.cpp
  assemble-matrix.cpp
  iteration.cpp
)

set(AHF_CORE_HEADERS
  apparent-horizon-finder.h
  peano-input.h
  output.h
  interpolation.h
  assemble-matrix.h
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
    VTK::IOXML
)

vtk_module_autoinit(
  TARGETS ApparentHorizonFinderPostprocess
  MODULES VTK::CommonCore VTK::CommonDataModel VTK::IOXML
)

add_executable(ahf-postprocess ahf-postprocess.cpp)
target_link_libraries(ahf-postprocess PRIVATE ApparentHorizonFinderPostprocess)
