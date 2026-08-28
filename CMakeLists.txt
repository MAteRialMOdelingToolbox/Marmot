#/* ---------------------------------------------------------------------
# *                                       _
# *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_
# * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
# * | | | | | | (_| | |  | | | | | | (_) | |_
# * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
# *
# * Unit of Strength of Materials and Structural Analysis
# * University of Innsbruck,
# * 2020 - today
# *
# * festigkeitslehre@uibk.ac.at
# *
# * This file is part of the MAteRialMOdellingToolbox (marmot).
# *
# * This library is free software; you can redistribute it and/or
# * modify it under the terms of the GNU Lesser General Public
# * License as published by the Free Software Foundation; either
# * version 2.1 of the License, or (at your option) any later version.
# *
# * The full text of the license can be found in the file LICENSE.md at
# * the top level directory of marmot.
# * ---------------------------------------------------------------------
# */

# cmake_minimum_required must be the very first command so that the correct
# CMake policies are active when the rest of the file is parsed.
cmake_minimum_required(VERSION 3.16)

message("+--------------------------------------------------------------------+")
message("|                                                                    |")
message("|                              MARMOT                                |")
message("|                                                                    |")
message("|                  -- MAteRialMOdellingToolbox --                    |")
message("|                                                                    |")
message("+--------------------------------------------------------------------+")
message(STATUS "This is CMake ${CMAKE_VERSION}")
message(STATUS "")

# Fortran is not used by Marmot; keep the language list to CXX only so that
# a Fortran compiler is not required.
project(Marmot VERSION 1.0.0 LANGUAGES CXX)

# avoid in-source compilation
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds are not allowed. Create and/or change to your build \
folder and run cmake from there, e.g.:\n  mkdir build && cd build && cmake ../\n\
You may need to remove CMakeCache.txt and CMakeFiles/ first.")
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# set install directories
include(GNUInstallDirs)

# Compiler-specific settings.
# Optimization/debug flags (-O3, -g, ...) come from CMake's per-config defaults.
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall -Wextra -Wno-unused-parameter -Wno-long-long -pedantic)
elseif(MSVC)
    if(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
        string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
    endif()
endif()

# If no build type was specified, Release is default
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
    message("--> setting to default CMAKE_BUILD_TYPE Release")
endif()

# Opt-in tuning for the host CPU. Fastor's SIMD kernels only vectorize beyond
# the SSE2 baseline when the target architecture allows it (AVX2/FMA), so this
# can matter a lot for tensor-heavy code paths. OFF by default because the
# resulting binaries are not portable to older CPUs.
option(MARMOT_MARCH_NATIVE "Compile with -march=native (non-portable binaries)" OFF)
if(MARMOT_MARCH_NATIVE)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        add_compile_options(-march=native)
    else()
        message(WARNING "MARMOT_MARCH_NATIVE is only supported for GCC/Clang and will be ignored.")
    endif()
endif()

# Disable (shared memory) parallelization within Marmot - such a parallelization is usually not desired
# add_compile_definitions(EIGEN_DONT_PARALLELIZE)                     # for Eigen
# add_compile_definitions(BLAZE_USE_SHARED_MEMORY_PARALLELIZATION=0)  # for Blaze

## find Eigen library
## The include directory is provided transitively via the Eigen3::Eigen imported target;
## no separate include_directories() call is needed.
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
message("--> found Eigen: ${EIGEN3_INCLUDE_DIR}")

## find autodiff library (https://autodiff.github.io/)
## autodiff ships a CMake package providing the autodiff::autodiff imported target.
find_package(autodiff REQUIRED)
message("--> found autodiff: ${autodiff_DIR}")

## find Fastor library (https://github.com/romeric/Fastor)
## Prefer Fastor's CMake package (Fastor::Fastor imported target); fall back to a
## plain header search for copies of Fastor installed without its package files.
find_package(Fastor QUIET)
if(TARGET Fastor::Fastor)
    message("--> found Fastor: ${Fastor_DIR}")
else()
    find_path(FASTOR_INCLUDE_DIR Fastor)
    if(NOT FASTOR_INCLUDE_DIR)
        message(FATAL_ERROR "Fastor not found (neither as a CMake package nor as plain headers)")
    endif()
    add_library(Fastor::Fastor INTERFACE IMPORTED)
    set_target_properties(Fastor::Fastor PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${FASTOR_INCLUDE_DIR}")
    message("--> found Fastor headers: ${FASTOR_INCLUDE_DIR}")
endif()

# Testing: include(CTest) provides the standard BUILD_TESTING option (default ON)
# and calls enable_testing() when it is enabled. Configure with -DBUILD_TESTING=OFF
# to skip building the test executables (e.g. for install-only builds).
include(CTest)

# Custom function to add a test executable.
# Include directories are propagated transitively from the Marmot library target
# (PUBLIC includes), so no explicit target_include_directories call is needed here.
function(add_marmot_test test_name test_source)
  add_executable(${test_name} ${test_source})
  target_link_libraries(${test_name} PRIVATE ${PROJECT_NAME})
  add_test(NAME ${test_name} COMMAND ${test_name} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endfunction()

# marmot_module_requires: check that every listed dependency of <module_name> is
# among the discovered INSTALLED_MODULES, so a module can skip itself gracefully
# instead of failing at link time. Emits a warning for each missing dependency
# and sets <result_var> to TRUE/FALSE in the calling scope.
#
# Usage (from a module.cmake):
#   marmot_module_requires(_deps_ok MyModule MarmotMechanicsCore MarmotCosseratCore)
#   if(_deps_ok)
#     # register include dir and sources
#   endif()
function(marmot_module_requires result_var module_name)
  set(_ok TRUE)
  foreach(_dep IN LISTS ARGN)
    if(NOT _dep IN_LIST INSTALLED_MODULES)
      message(WARNING "Module ${module_name} will NOT be built: required module '${_dep}' is not enabled.")
      set(_ok FALSE)
    endif()
  endforeach()
  set(${result_var} ${_ok} PARENT_SCOPE)
endfunction()

set(MODULES_DIR ${CMAKE_SOURCE_DIR}/modules)
set(SHARED_LIBRARIES "")

# Accumulated lists populated by module.cmake files
set(sources "")
set(INSTALLED_MODULE_INCLUDE_DIRS "")

# ── Module discovery helpers ──────────────────────────────────────────────────

# subdirlist: populate <result> with the names of immediate subdirectories of <curdir>.
function(subdirlist result curdir)
  file(GLOB _children CONFIGURE_DEPENDS RELATIVE "${curdir}" "${curdir}/*")
  set(_dirs "")
  foreach(_child ${_children})
    if(IS_DIRECTORY "${curdir}/${_child}")
      list(APPEND _dirs "${_child}")
    endif()
  endforeach()
  set(${result} ${_dirs} PARENT_SCOPE)
endfunction()

# register_module_category: scan <category_path> for subdirectories that contain a
# module.cmake file, apply the optional filter list stored in the variable named
# <filter_var>, and append matching entries to INSTALLED_MODULES and
# INSTALLED_MODULEPATHS in the calling scope.
#
# Usage:
#   register_module_category(<path> <filter-variable-name> "<label>")
#
# A module is included when either:
#   - the filter variable is undefined (defaults to "all"), or
#   - the filter variable equals "all", or
#   - the directory name appears in the filter list.
macro(register_module_category category_path filter_var category_label)
  message("+------------------------------------------------------------------------------+")
  message("Checking for ${category_label} modules in ${category_path}.")
  if(EXISTS "${category_path}")
    subdirlist(_rmc_dirs "${category_path}")
    foreach(_rmc_dir IN LISTS _rmc_dirs)
      set(_rmc_modulepath "${category_path}/${_rmc_dir}")
      if(EXISTS "${_rmc_modulepath}/module.cmake")
        if(NOT DEFINED ${filter_var})
          set(${filter_var} "all")
        endif()
        if("${${filter_var}}" STREQUAL "all" OR ${_rmc_dir} IN_LIST ${filter_var})
          message("--> found ${_rmc_dir}")
          list(APPEND INSTALLED_MODULES "${_rmc_dir}")
          list(APPEND INSTALLED_MODULEPATHS "${_rmc_modulepath}")
        endif()
      endif()
    endforeach()
    unset(_rmc_dirs)
    unset(_rmc_dir)
    unset(_rmc_modulepath)
  endif()
endmacro()

# ── Discover modules ──────────────────────────────────────────────────────────

set(INSTALLED_MODULES "")
set(INSTALLED_MODULEPATHS "")

register_module_category("${MODULES_DIR}/core"           CORE_MODULES          "Core")
register_module_category("${MODULES_DIR}/materials"      MATERIAL_MODULES      "Material")
register_module_category("${MODULES_DIR}/elements"       ELEMENT_MODULES       "Element")
register_module_category("${MODULES_DIR}/particles"      PARTICLE_MODULES      "Particle")
register_module_category("${MODULES_DIR}/materialpoints" MATERIALPOINT_MODULES "MaterialPoint")
register_module_category("${MODULES_DIR}/cells"          CELL_MODULES          "Cell")
register_module_category("${MODULES_DIR}/cellelements"   CELLELEMENT_MODULES   "CellElement")

message("+------------------------------------------------------------------------------+")
message("Check the list above: if a module is missing despite being in the correct")
message("directory, verify that it has a proper module.cmake file.")
message("+------------------------------------------------------------------------------+")

# Include all discovered module.cmake files.
# Each file appends to: sources, INSTALLED_MODULE_INCLUDE_DIRS.
foreach(MODULEPATH IN LISTS INSTALLED_MODULEPATHS)
    include("${MODULEPATH}/module.cmake")
endforeach()

message("+------------------------------------------------------------------------------+")
message(STATUS "CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")

# ── Build the shared library ──────────────────────────────────────────────────

add_library(${PROJECT_NAME} SHARED ${sources})

set_target_properties(${PROJECT_NAME} PROPERTIES
    VERSION   ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
)

# Include-directory strategy:
#   BUILD_INTERFACE  – source-tree paths used when building Marmot and its tests.
#   INSTALL_INTERFACE – the single installed include root used by find_package consumers.
# Making these PUBLIC means test executables that link against ${PROJECT_NAME} receive
# the include paths transitively without any explicit target_include_directories call.
# Dependency include paths come transitively from the linked imported targets.
target_include_directories(${PROJECT_NAME}
    PUBLIC
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
foreach(_inc_dir IN LISTS INSTALLED_MODULE_INCLUDE_DIRS)
    target_include_directories(${PROJECT_NAME}
        PUBLIC $<BUILD_INTERFACE:${_inc_dir}>)
endforeach()

# PUBLIC: the installed Marmot headers include Eigen, autodiff, and Fastor headers,
# so consumers need these interface targets too (resolved via MarmotConfig.cmake).
target_link_libraries(${PROJECT_NAME} PUBLIC Eigen3::Eigen autodiff::autodiff Fastor::Fastor)

if(SHARED_LIBRARIES)
    message("+------------------------------------------------------------------------------+")
    message("Marmot modules have requested to link against following shared libraries:")
    foreach(SHARED_LIBRARY IN LISTS SHARED_LIBRARIES)
        message("--> ${SHARED_LIBRARY}")
        target_link_libraries(${PROJECT_NAME} PRIVATE "${SHARED_LIBRARY}")
    endforeach()
    message("+------------------------------------------------------------------------------+")
endif()

# Include test definitions for all installed modules
if(BUILD_TESTING)
    foreach(MODULEPATH IN LISTS INSTALLED_MODULEPATHS)
        if(EXISTS "${MODULEPATH}/test.cmake")
            include("${MODULEPATH}/test.cmake")
        endif()
    endforeach()
endif()

# ── Install ───────────────────────────────────────────────────────────────────

install(TARGETS ${PROJECT_NAME}
    EXPORT  MarmotTargets
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
)

# Install every module's include tree (headers are laid out as
# <module>/include/Marmot/<Header>.h, so the Marmot/ prefix is preserved).
# This installs all module headers wholesale instead of relying on
# hand-maintained per-module header lists.
foreach(MODULEPATH IN LISTS INSTALLED_MODULEPATHS)
    if(EXISTS "${MODULEPATH}/include")
        install(DIRECTORY "${MODULEPATH}/include/"
            DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
            FILES_MATCHING PATTERN "*.h"
        )
    endif()
endforeach()

install(EXPORT MarmotTargets
    FILE        MarmotTargets.cmake
    NAMESPACE   Marmot::
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/Marmot"
)

include(CMakePackageConfigHelpers)

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/MarmotConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/MarmotConfig.cmake"
    INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/Marmot"
)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/MarmotConfigVersion.cmake"
    VERSION       ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/MarmotConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/MarmotConfigVersion.cmake"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/Marmot"
)
