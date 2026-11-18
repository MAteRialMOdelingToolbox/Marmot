#/* ---------------------------------------------------------------------
# *
# * Locates Marmot's header-only dependencies and provides them as imported targets:
# *
# *   Eigen3::Eigen        (Eigen's CMake package)
# *   autodiff::autodiff   (autodiff's CMake package, or a plain header search)
# *   Fastor::Fastor       (Fastor's CMake package, or a plain header search)
# *   Sacado::Sacado       (only if MARMOT_AD_VECTOR_BACKEND is Sacado: Trilinos' Sacado CMake package, or a
# *                         plain header search)
# *
# * Shared by Marmot's top-level CMakeLists.txt and the installed MarmotConfig.cmake, so that
# * Marmot and its consumers resolve the dependencies the same way.
# *
# * The header-search fallbacks cover copies installed without CMake package files (e.g. a
# * Fastor installed by copying its headers).
# *
# * Sets MARMOT_DEPENDENCIES_FOUND to TRUE or FALSE and, on failure,
# * MARMOT_DEPENDENCIES_MESSAGE; it never stops configuring itself, the includer decides.
# *
# * ---------------------------------------------------------------------
# */

set(MARMOT_DEPENDENCIES_FOUND TRUE)
set(MARMOT_DEPENDENCIES_MESSAGE "")

# Eigen: no version is passed to find_package(), so that an Eigen 3 found instead is reported with its
# version below rather than as "not found". Marmot requires Eigen 5.
find_package(Eigen3 QUIET NO_MODULE)
if(NOT TARGET Eigen3::Eigen)
    set(MARMOT_DEPENDENCIES_FOUND FALSE)
    string(APPEND MARMOT_DEPENDENCIES_MESSAGE "Eigen3 (CMake package) not found. ")
elseif(Eigen3_VERSION VERSION_LESS 5.0)
    set(MARMOT_DEPENDENCIES_FOUND FALSE)
    string(APPEND MARMOT_DEPENDENCIES_MESSAGE "Marmot requires Eigen >= 5.0, found ${Eigen3_VERSION}. ")
endif()

# marmot_find_header_only_dependency: provide <target> from <package>'s CMake package if it has
# one, otherwise from a search for the header directory <header>.
macro(marmot_find_header_only_dependency package target header)
    if(NOT TARGET ${target})
        find_package(${package} QUIET CONFIG)
    endif()
    if(NOT TARGET ${target})
        find_path(MARMOT_${package}_INCLUDE_DIR ${header})
        if(MARMOT_${package}_INCLUDE_DIR)
            add_library(${target} INTERFACE IMPORTED)
            set_target_properties(${target} PROPERTIES
                INTERFACE_INCLUDE_DIRECTORIES "${MARMOT_${package}_INCLUDE_DIR}")
        else()
            set(MARMOT_DEPENDENCIES_FOUND FALSE)
            string(APPEND MARMOT_DEPENDENCIES_MESSAGE
                "${package} not found (neither as a CMake package nor as plain headers). ")
        endif()
    endif()
endmacro()

marmot_find_header_only_dependency(autodiff autodiff::autodiff autodiff)
marmot_find_header_only_dependency(Fastor Fastor::Fastor Fastor)

# Sacado (optional vector-mode AD backend, see MarmotAutomaticDifferentiationTypes.h). Only its header-only part
# is used, so Sacado::Sacado carries just the include directories, not Trilinos' (static) sacado library. An
# installed Trilinos provides the package Sacado with the target Sacado::all_libs; a bare copy of the installed
# headers (including the generated Sacado_config.h) is found by a header search.
if(MARMOT_AD_VECTOR_BACKEND STREQUAL "Sacado")
    if(NOT TARGET Sacado::Sacado)
        find_package(Sacado QUIET CONFIG)
        if(TARGET Sacado::all_libs)
            get_target_property(_marmot_sacado_include_dirs Sacado::all_libs INTERFACE_INCLUDE_DIRECTORIES)
            add_library(Sacado::Sacado INTERFACE IMPORTED)
            set_target_properties(Sacado::Sacado PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${_marmot_sacado_include_dirs}")
        endif()
    endif()
    if(NOT TARGET Sacado::Sacado)
        find_path(MARMOT_Sacado_INCLUDE_DIR Sacado_No_Kokkos.hpp)
        if(MARMOT_Sacado_INCLUDE_DIR AND EXISTS "${MARMOT_Sacado_INCLUDE_DIR}/Sacado_config.h")
            add_library(Sacado::Sacado INTERFACE IMPORTED)
            set_target_properties(Sacado::Sacado PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${MARMOT_Sacado_INCLUDE_DIR}")
        else()
            set(MARMOT_DEPENDENCIES_FOUND FALSE)
            string(APPEND MARMOT_DEPENDENCIES_MESSAGE
                "Sacado not found (neither as a CMake package nor as installed headers with Sacado_config.h), "
                "but MARMOT_AD_VECTOR_BACKEND is Sacado. ")
        endif()
    endif()
elseif(NOT MARMOT_AD_VECTOR_BACKEND STREQUAL "autodiff")
    set(MARMOT_DEPENDENCIES_FOUND FALSE)
    string(APPEND MARMOT_DEPENDENCIES_MESSAGE
        "MARMOT_AD_VECTOR_BACKEND must be autodiff or Sacado, not '${MARMOT_AD_VECTOR_BACKEND}'. ")
endif()
