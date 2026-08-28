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

# Compiler-specific settings
if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
    add_compile_options(-Wall -Wextra -Wno-unused-parameter -Wno-long-long -pedantic)
    # -O3 belongs only in Release builds, not Debug.
    add_compile_options($<$<CONFIG:Release>:-O3>)
    # add_compile_options(-mavx)
    # add_compile_options(-ffast-math)
    add_compile_options($<$<CONFIG:Debug>:-g>)
elseif(MSVC)
    if(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
        string(REGEX REPLACE "/W[0-4]" "/W0" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W0")
    endif()
    add_compile_definitions(_CRT_SECURE_NO_WARNINGS)
    add_compile_definitions(__PRETTY_FUNCTION__=__FUNCSIG__)
endif()

# If no build type was specified, Release is default
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
    message("--> setting to default CMAKE_BUILD_TYPE Release")
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
find_path(AUTODIFF_INCLUDE_DIR autodiff)
if(NOT AUTODIFF_INCLUDE_DIR)
    message(FATAL_ERROR "autodiff include directory not found")
endif()
message("--> found autodiff: ${AUTODIFF_INCLUDE_DIR}")

## find Fastor library (https://github.com/romeric/Fastor)
find_path(FASTOR_INCLUDE_DIR Fastor)
if(NOT FASTOR_INCLUDE_DIR)
    message(FATAL_ERROR "Fastor include directory not found")
endif()
message("--> found Fastor: ${FASTOR_INCLUDE_DIR}")

# enable testing
enable_testing()

# Custom function to add a test executable.
# Include directories are propagated transitively from the Marmot library target
# (PUBLIC includes), so no explicit target_include_directories call is needed here.
function(add_marmot_test test_name test_source)
  add_executable(${test_name} ${test_source})
  if(MSVC AND BUILD_TESTING)
    target_link_libraries(${test_name} PRIVATE -WHOLEARCHIVE:${PROJECT_NAME})
  else()
    target_link_libraries(${test_name} PRIVATE ${PROJECT_NAME})
  endif()
  add_test(NAME ${test_name} COMMAND ${test_name} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endfunction()

set(MODULES_DIR ${CMAKE_SOURCE_DIR}/modules)
set(SHARED_LIBRARIES "")

# Accumulated lists populated by module.cmake files
set(sources "")
set(publicheaders "")
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
# Each file appends to: sources, publicheaders, INSTALLED_MODULE_INCLUDE_DIRS.
foreach(MODULEPATH IN LISTS INSTALLED_MODULEPATHS)
    include("${MODULEPATH}/module.cmake")
endforeach()

message("+------------------------------------------------------------------------------+")
message(STATUS "CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")

# Build as STATIC on Windows when testing is enabled to allow tests to link internal symbols
if(MSVC AND BUILD_TESTING)
    add_library(${PROJECT_NAME} STATIC ${sources})
else()
    add_library(${PROJECT_NAME} SHARED ${sources})
endif()

# Generate Windows Export Header
include(GenerateExportHeader)
generate_export_header(${PROJECT_NAME}
    BASE_NAME MARMOT
    EXPORT_FILE_NAME "${CMAKE_BINARY_DIR}/include/Marmot/marmot_export.h"
)

# Append the generated export header to public headers so it gets installed
list(APPEND publicheaders "${CMAKE_BINARY_DIR}/include/Marmot/marmot_export.h")

set_target_properties(${PROJECT_NAME} PROPERTIES
    VERSION       ${PROJECT_VERSION}
    SOVERSION     ${PROJECT_VERSION_MAJOR}
    PUBLIC_HEADER "${publicheaders}"
)

# Include-directory strategy:
#   BUILD_INTERFACE  – source-tree paths used when building Marmot and its tests.
#   INSTALL_INTERFACE – the single installed include root used by find_package consumers.
# Making these PUBLIC means test executables that link against ${PROJECT_NAME} receive
# the include paths transitively without any explicit target_include_directories call.
target_include_directories(${PROJECT_NAME}
    PUBLIC
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
        $<BUILD_INTERFACE:${AUTODIFF_INCLUDE_DIR}>
        $<BUILD_INTERFACE:${FASTOR_INCLUDE_DIR}>
        $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
)
foreach(_inc_dir IN LISTS INSTALLED_MODULE_INCLUDE_DIRS)
    target_include_directories(${PROJECT_NAME}
        PUBLIC $<BUILD_INTERFACE:${_inc_dir}>)
endforeach()

target_link_libraries(${PROJECT_NAME} PUBLIC Eigen3::Eigen)

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
foreach(MODULEPATH IN LISTS INSTALLED_MODULEPATHS)
    if(EXISTS "${MODULEPATH}/test.cmake")
        include("${MODULEPATH}/test.cmake")
    endif()
endforeach()

# ── Install ───────────────────────────────────────────────────────────────────

install(TARGETS ${PROJECT_NAME}
    EXPORT        MarmotTargets
    LIBRARY       DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/Marmot"
)

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
