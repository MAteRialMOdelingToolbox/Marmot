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
# The release this source tree is (or leads up to), as in the release tag (v26.11 -> 26.11.0); set to the next release
# right after tagging. The conda workflow refuses a tag that does not match it.
project(Marmot VERSION 26.11.0 LANGUAGES CXX)

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
set(CMAKE_BUILD_RPATH "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# set install directories
include(GNUInstallDirs)

# Compiler-specific settings.
# Optimization/debug flags (-O3, -g, ...) come from CMake's per-configuration defaults.
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall -Wextra -Wno-unused-parameter -Wno-long-long -pedantic)
elseif(MSVC)
    if(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
        string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
    endif()
    # /bigobj: the template-heavy (Eigen/Fastor/autodiff) translation units exceed MSVC's default section limit.
    # /permissive-: standard-conforming two-phase lookup, as GCC/Clang do.
    # /Zc:__cplusplus: report the real language level in __cplusplus, which Eigen and Fastor check.
    add_compile_options(/bigobj /permissive- /Zc:__cplusplus)
    # Keep MSVC's secure-CRT deprecation warnings (C4996 for std::getenv, sscanf, ...) out of the /W4 output.
    add_compile_definitions(_CRT_SECURE_NO_WARNINGS)
endif()

# If no build type was specified, Release is default
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
    message("--> setting to default CMAKE_BUILD_TYPE Release")
endif()

# Opt-in tuning for the host CPU. Fastor's SIMD kernels only vectorize beyond the SSE2 baseline when the
# target architecture allows it (AVX2/FMA). OFF by default because the binaries are not portable.
option(MARMOT_MARCH_NATIVE "Compile with -march=native (non-portable binaries)" OFF)
if(MARMOT_MARCH_NATIVE)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        add_compile_options(-march=native)
    else()
        message(WARNING "MARMOT_MARCH_NATIVE is only supported for GCC/Clang and is ignored.")
    endif()
endif()

# Disable (shared memory) parallelization within Marmot - such a parallelization is usually not desired
# add_compile_definitions(EIGEN_DONT_PARALLELIZE)                     # for Eigen
# add_compile_definitions(BLAZE_USE_SHARED_MEMORY_PARALLELIZATION=0)  # for Blaze

## find Eigen, autodiff and Fastor as imported targets (Eigen3::Eigen, autodiff::autodiff, Fastor::Fastor).
## The same file is installed and used by MarmotConfig.cmake, so consumers resolve them alike.
include(cmake/MarmotDependencies.cmake)
if(NOT MARMOT_DEPENDENCIES_FOUND)
    message(FATAL_ERROR "${MARMOT_DEPENDENCIES_MESSAGE}")
endif()
message("--> found Eigen ${Eigen3_VERSION}")

# Testing: include(CTest) provides the standard BUILD_TESTING option (default ON) and calls enable_testing().
# Configure with -DBUILD_TESTING=OFF to skip the test executables, e.g. for install-only builds.
include(CTest)

# ── Options affecting the module objects ─────────────────────────────────────

# The library exports only what is marked MARMOT_API (see MarmotPortability.h), as a Windows DLL always does; on
# other platforms only with MARMOT_EXPORT_API_ONLY, to check the exported API without a Windows build.
option(MARMOT_EXPORT_API_ONLY "Export only the MARMOT_API-marked API from libMarmot, as on Windows" OFF)

# Opt-in: instruments the Marmot library with gcov. Test executables are instrumented
# too (see add_marmot_test) because templated classes get re-instantiated per
# translation unit, and an uninstrumented test binary's own copy of a template method
# otherwise wins symbol resolution over the library's instrumented one, hiding it from
# coverage entirely. Flags are kept PRIVATE: PUBLIC link options would be exported into
# MarmotTargets.cmake, silently forcing --coverage onto any downstream consumer that
# links Marmot::Marmot from a coverage-enabled install.
option(MARMOT_ENABLE_COVERAGE "Instrument libMarmot and its tests with gcov code coverage (GCC/Clang only)" OFF)

if(MARMOT_ENABLE_COVERAGE)
    if(NOT CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        message(FATAL_ERROR "MARMOT_ENABLE_COVERAGE requires GCC or Clang")
    endif()
    if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
        message(WARNING "MARMOT_ENABLE_COVERAGE is intended to be used with -DCMAKE_BUILD_TYPE=Debug for accurate line coverage")
    endif()
    message(STATUS "Instrumenting ${PROJECT_NAME} for code coverage (--coverage)")
endif()

include(cmake/MarmotModules.cmake)

set(MODULES_DIR ${CMAKE_SOURCE_DIR}/modules)

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
# INSTALLED_MODULEPATHS in the calling scope. Every module found, filtered or not, is
# appended to MARMOT_EXISTING_MODULES; the modules named in a filter to MARMOT_EXPLICIT_MODULES.
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
  if(NOT DEFINED ${filter_var})
    set(${filter_var} "all")
  endif()
  if(NOT "${${filter_var}}" STREQUAL "all")
    set(MARMOT_MODULES_FILTERED TRUE)
  endif()
  if(EXISTS "${category_path}")
    subdirlist(_rmc_dirs "${category_path}")
    foreach(_rmc_dir IN LISTS _rmc_dirs)
      set(_rmc_modulepath "${category_path}/${_rmc_dir}")
      if(EXISTS "${_rmc_modulepath}/module.cmake")
        list(APPEND MARMOT_EXISTING_MODULES "${_rmc_dir}")
        if("${${filter_var}}" STREQUAL "all" OR ${_rmc_dir} IN_LIST ${filter_var})
          message("--> found ${_rmc_dir}")
          list(APPEND INSTALLED_MODULES "${_rmc_dir}")
          if(NOT "${${filter_var}}" STREQUAL "all")
            list(APPEND MARMOT_EXPLICIT_MODULES "${_rmc_dir}")
          endif()
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
set(MARMOT_EXPLICIT_MODULES "")
set(MARMOT_EXISTING_MODULES "")
set(MARMOT_MODULES_FILTERED FALSE)

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

# Read all discovered module.cmake files, resolve the dependencies and create one target per module
# (see cmake/MarmotModules.cmake).
marmot_read_modules("${INSTALLED_MODULEPATHS}")
marmot_build_modules("${INSTALLED_MODULES}" "${MARMOT_EXISTING_MODULES}" "${MARMOT_EXPLICIT_MODULES}"
                     ${MARMOT_MODULES_FILTERED})

message("+------------------------------------------------------------------------------+")
message(STATUS "Building modules: ${MARMOT_BUILT_MODULES}")
message(STATUS "CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")

# ── Build the shared library ──────────────────────────────────────────────────

# All module objects, compiled once, from which both the shared library and the tests are linked.
set(MARMOT_MODULE_OBJECTS "")
foreach(_target IN LISTS MARMOT_MODULE_TARGETS)
    get_target_property(_type ${_target} TYPE)
    if(_type STREQUAL "OBJECT_LIBRARY")
        list(APPEND MARMOT_MODULE_OBJECTS $<TARGET_OBJECTS:${_target}>)
    endif()
endforeach()
if(NOT MARMOT_MODULE_OBJECTS)
    message(FATAL_ERROR "No module with sources is selected; check the *_MODULES filters.")
endif()

add_library(${PROJECT_NAME} SHARED ${MARMOT_MODULE_OBJECTS})

set_target_properties(${PROJECT_NAME} PROPERTIES
    VERSION   ${PROJECT_VERSION}
    # MAJOR.MINOR (libMarmot.so.26.11): the major number alone (26) is shared by the 26.05 and 26.11 releases, whose
    # ABIs differ, so the loader must not treat them as interchangeable.
    SOVERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
)

# Marmot's headers need these MSVC options in every translation unit that includes them (Eigen and Fastor check
# __cplusplus, which MSVC only reports correctly with /Zc:__cplusplus), so consumers of Marmot::Marmot get them too.
target_compile_options(${PROJECT_NAME} INTERFACE "$<$<CXX_COMPILER_ID:MSVC>:/bigobj;/permissive-;/Zc:__cplusplus>")

# Include-directory strategy:
#   BUILD_INTERFACE  – source-tree paths, for targets of this build linking ${PROJECT_NAME} (e.g. the bindings).
#   INSTALL_INTERFACE – the single installed include root used by find_package consumers.
# The shared library does not link the module targets: MARMOT_BUILDING_LIBRARY and the coverage flags stay
# private to the module objects.
target_include_directories(${PROJECT_NAME} PUBLIC $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
foreach(_inc_dir IN LISTS MARMOT_INCLUDE_DIRS)
    target_include_directories(${PROJECT_NAME} PUBLIC $<BUILD_INTERFACE:${_inc_dir}>)
endforeach()

# PUBLIC: Marmot's headers include Eigen, autodiff and Fastor, so consumers need them too
# (MarmotConfig.cmake recreates these targets on the consumer side).
target_link_libraries(${PROJECT_NAME} PUBLIC Eigen3::Eigen autodiff::autodiff Fastor::Fastor)
target_link_libraries(${PROJECT_NAME} PRIVATE ${MARMOT_MODULE_LINK_LIBRARIES})

if(MARMOT_ENABLE_COVERAGE)
    target_link_options(${PROJECT_NAME} PRIVATE --coverage)
endif()

# ── Tests ─────────────────────────────────────────────────────────────────────

# Custom function to add a test executable.
# Tests link Marmot's object files directly, not the shared library: they also exercise Marmot's internals, which
# the library does not export on Windows (or with MARMOT_EXPORT_API_ONLY), see MARMOT_API in MarmotPortability.h.
# The same on every platform, so that the tests are built alike everywhere. Whether a consumer can link the shared
# library is checked by the one test that does, tests/consumer.
# Many Marmot classes are templates (<int nDim, int nNodes>) that get implicitly re-instantiated inside each test's
# own translation unit; with MARMOT_ENABLE_COVERAGE the tests are instrumented too (through MarmotModuleSettings),
# as otherwise an uninstrumented test binary's own copy of such a template method wins symbol resolution over the
# library's instrumented copy, so gcov never sees it.
# A test sees the headers of its own module (and the modules that one requires) and of the module providing
# MarmotTesting.h; a test including headers of further modules, or creating their materials or elements by name
# through the factories, lists them, and is skipped if one is not built:
#   add_marmot_test(<name> <source> [REQUIRES <module>...] [TIMEOUT <seconds>])
# Every test has a timeout (default MARMOT_TEST_TIMEOUT, well below ctest's implicit 1500 s), so that a test which
# becomes too slow fails as such instead of passing on a fast runner and timing out on a slow one; it must hold also
# for the -O0 coverage build.
set(MARMOT_TESTING_MODULE MarmotMechanicsCore)
set(MARMOT_TEST_TIMEOUT 300 CACHE STRING "Default timeout of each Marmot test in seconds")
function(add_marmot_test test_name test_source)
  cmake_parse_arguments(PARSE_ARGV 2 _arg "" "TIMEOUT" "REQUIRES")
  if(_arg_UNPARSED_ARGUMENTS)
    message(FATAL_ERROR "add_marmot_test(${test_name}): unknown arguments: ${_arg_UNPARSED_ARGUMENTS}")
  endif()
  set(_modules ${MARMOT_CURRENT_MODULE} ${MARMOT_TESTING_MODULE} ${_arg_REQUIRES})
  foreach(_module IN LISTS _modules)
    if(NOT _module IN_LIST MARMOT_BUILT_MODULES)
      message(STATUS "Test ${test_name} skipped: it requires module ${_module}, which is not built")
      return()
    endif()
  endforeach()
  set(_usage "")
  foreach(_module IN LISTS _modules)
    marmot_module_target(_target ${_module})
    list(APPEND _usage ${_target})
  endforeach()
  list(REMOVE_DUPLICATES _usage)
  # Linking a module target directly links its objects; the objects of all other modules are added explicitly.
  set(_objects "")
  foreach(_target IN LISTS MARMOT_MODULE_TARGETS)
    get_target_property(_type ${_target} TYPE)
    if(_type STREQUAL "OBJECT_LIBRARY" AND NOT _target IN_LIST _usage)
      list(APPEND _objects $<TARGET_OBJECTS:${_target}>)
    endif()
  endforeach()
  add_executable(${test_name} ${test_source} ${_objects})
  # the objects of all modules are linked, so their libraries are needed too
  target_link_libraries(${test_name} PRIVATE ${_usage} ${MARMOT_MODULE_LINK_LIBRARIES})
  add_test(NAME ${test_name} COMMAND ${test_name} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  if(NOT _arg_TIMEOUT)
    set(_arg_TIMEOUT ${MARMOT_TEST_TIMEOUT})
  endif()
  set_tests_properties(${test_name} PROPERTIES TIMEOUT ${_arg_TIMEOUT})
  if(MARMOT_ENABLE_COVERAGE)
    target_link_options(${test_name} PRIVATE --coverage)
  endif()
endfunction()

if(BUILD_TESTING)
    # Include the test definitions of the built modules only.
    foreach(MARMOT_CURRENT_MODULE IN LISTS MARMOT_BUILT_MODULES)
        get_property(_dir GLOBAL PROPERTY MARMOT_MODULE_${MARMOT_CURRENT_MODULE}_DIR)
        if(EXISTS "${_dir}/test.cmake")
            include("${_dir}/test.cmake")
        endif()
    endforeach()

    # The consumer test, linking the shared library only (see tests/consumer/CMakeLists.txt). It creates a CPS4
    # element with a LINEARELASTIC material through the factories, so it needs both modules.
    if("DisplacementFiniteElement" IN_LIST MARMOT_BUILT_MODULES AND "LinearElastic" IN_LIST MARMOT_BUILT_MODULES)
        add_subdirectory(tests/consumer)
    else()
        message(STATUS "tests/consumer skipped: it requires the DisplacementFiniteElement and LinearElastic modules")
    endif()
endif()

# ── Python bindings ───────────────────────────────────────────────────────────

# Opt-in: building the bindings requires the Python development headers and
# fetches nanobind, neither of which should be a hard dependency of the C++ library.
option(MARMOT_BUILD_PYTHON_BINDINGS "Build the Marmot Python bindings" OFF)

if(MARMOT_BUILD_PYTHON_BINDINGS)
    # The bindings wrap Marmot API that is not marked MARMOT_API, which the shared library they link does not export
    # on Windows or with MARMOT_EXPORT_API_ONLY; they would fail to link (or to import), so fail early instead.
    if(WIN32 OR MARMOT_EXPORT_API_ONLY)
        message(FATAL_ERROR "The Marmot Python bindings are not supported on Windows or with MARMOT_EXPORT_API_ONLY: "
                            "they use API that the library does not export then. Build them on Linux or macOS "
                            "without MARMOT_EXPORT_API_ONLY.")
    endif()
    message(STATUS "Building the Marmot Python bindings")
    add_subdirectory(modules/python)
endif()

# ── Install ───────────────────────────────────────────────────────────────────

install(TARGETS ${PROJECT_NAME}
    EXPORT        MarmotTargets
    LIBRARY       DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    # Windows: the DLL is a RUNTIME artifact, its import library (.lib) an ARCHIVE one.
    RUNTIME       DESTINATION "${CMAKE_INSTALL_BINDIR}"
    ARCHIVE       DESTINATION "${CMAKE_INSTALL_LIBDIR}"
)

# Install the include tree of every built module wholesale. Headers are laid out as
# <module>/include/Marmot/<Header>.h, so they all land in <prefix>/include/Marmot/.
foreach(_inc_dir IN LISTS MARMOT_INCLUDE_DIRS)
    install(DIRECTORY "${_inc_dir}/"
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
        FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp"
    )
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
    # Releases are versioned YY.MM, and each has its own ABI: find_package(Marmot 26.11) accepts 26.11.x only.
    COMPATIBILITY SameMinorVersion
)

install(FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/MarmotDependencies.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/MarmotConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/MarmotConfigVersion.cmake"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/Marmot"
)
