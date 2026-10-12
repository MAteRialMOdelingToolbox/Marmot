#/* ---------------------------------------------------------------------
# *
# * Locates Marmot's header-only dependencies and provides them as imported targets:
# *
# *   Eigen3::Eigen        (Eigen's CMake package)
# *   autodiff::autodiff   (autodiff's CMake package, or a plain header search)
# *   Fastor::Fastor       (Fastor's CMake package, or a plain header search)
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

# Eigen: no version is passed to find_package(): Eigen 3.4 and Eigen 5 both install SameMajorVersion
# config files, so neither a minimum of 3.3 nor a 3.3...5 range accepts both. Check the minimum here.
find_package(Eigen3 QUIET NO_MODULE)
if(NOT TARGET Eigen3::Eigen)
    set(MARMOT_DEPENDENCIES_FOUND FALSE)
    string(APPEND MARMOT_DEPENDENCIES_MESSAGE "Eigen3 (CMake package) not found. ")
elseif(Eigen3_VERSION VERSION_LESS 3.3)
    set(MARMOT_DEPENDENCIES_FOUND FALSE)
    string(APPEND MARMOT_DEPENDENCIES_MESSAGE "Marmot requires Eigen >= 3.3, found ${Eigen3_VERSION}. ")
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
