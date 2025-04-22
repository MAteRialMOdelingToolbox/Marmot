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
# * Matthias Neuner matthias.neuner@uibk.ac.at
# * Magdalena Schreter magdalena.schreter@uibk.ac.at
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

MESSAGE("+--------------------------------------------------------------------+")
MESSAGE("|                                                                    |")
MESSAGE("|                              MARMOT                                |")
MESSAGE("|                                                                    |")
MESSAGE("|                  -- MAteRialMOdellingToolbox --                    |")
MESSAGE("|                                                                    |")
MESSAGE("+--------------------------------------------------------------------+")

MESSAGE(STATUS "This is CMake ${CMAKE_VERSION}")
MESSAGE(STATUS "")

cmake_minimum_required(VERSION 3.12)

project(Marmot CXX Fortran)

# avoid in-source compilation
if (${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds are not allowed. Create and/or change to your build folder \
  and run then cmake from there, e.g. mkdir build && cd build && cmake ../ . \\n
  You first may need to remove CMakeCache.txt and CMakeFiles/*.")
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_FIND_LIBRARY_PREFIXES lib)

# enable testing
enable_testing()

# custom function to add tests
function(add_marmot_test test_name test_source)
  add_executable(${test_name} ${test_source})
  target_include_directories(${test_name} PRIVATE "${CMAKE_SOURCE_DIR}/include")
  target_include_directories(${test_name} PRIVATE "${INSTALLED_MODULEPATHS}/include")
  target_include_directories(${test_name} PRIVATE include)
  target_link_libraries(${test_name} PRIVATE ${PROJECT_NAME})
  add_test(NAME ${test_name} COMMAND ${test_name} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endfunction()

# set install directories
include(GNUInstallDirs)

# Compiler specific stettings
if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-parameter")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-long-long")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic")
    # only activate if you know what you are doing:
    # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 ")
    # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx")
    # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
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

# Disable (shared memory) parallelization within Marmot - such a parallelization is usually not desired
# add_compile_definitions(EIGEN_DONT_PARALLELIZE)                     # for Eigen
# add_compile_definitions(BLAZE_USE_SHARED_MEMORY_PARALLELIZATION=0)  # for Blaze

## find Eigen library
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
include_directories(${EIGEN3_INCLUDE_DIR})
message("--> found Eigen: ${EIGEN3_INCLUDE_DIR}")

## find autodiff library (https://autodiff.github.io/)
find_path(AUTODIFF_INCLUDE_DIR autodiff REQUIRED NO_MODULE)
include_directories(${AUTODIFF_INCLUDE_DIR})
message("--> found autodiff: ${AUTODIFF_INCLUDE_DIR}")

set(MODULES_DIR ${CMAKE_SOURCE_DIR}/modules)
set(SHARED_LIBRARIES "")
#
# Add sources
#
file(GLOB sources "src/*.cpp")
file(GLOB publicheaders "include/Marmot/*.h")



# INSTALL MODULES 

MACRO(SUBDIRLIST result curdir)
  FILE(GLOB children RELATIVE ${curdir} ${curdir}/*)
  SET(dirlist "")
  FOREACH(child ${children})
    IF(IS_DIRECTORY ${curdir}/${child})
      LIST(APPEND dirlist ${child})
    ENDIF()
  ENDFOREACH()
  SET(${result} ${dirlist})
ENDMACRO()

get_filename_component(COREPATH ${MODULES_DIR}/core ABSOLUTE) 
get_filename_component(MATERIALSPATH ${MODULES_DIR}/materials ABSOLUTE) 
get_filename_component(PARTICLESPATH ${MODULES_DIR}/particles ABSOLUTE) 
get_filename_component(MATERIALPOINTSPATH ${MODULES_DIR}/materialpoints ABSOLUTE) 
get_filename_component(ELEMENTSPATH ${MODULES_DIR}/elements ABSOLUTE) 
get_filename_component(CELLSPATH ${MODULES_DIR}/cells ABSOLUTE) 
get_filename_component(CELLELEMENTSPATH ${MODULES_DIR}/cellelements ABSOLUTE) 

SUBDIRLIST(COREDIRS ${COREPATH})
SUBDIRLIST(MATERIALSDIRS ${MATERIALSPATH})
SUBDIRLIST(PARTICLESDIRS ${PARTICLESPATH})
SUBDIRLIST(MATERIALPOINTSDIRS ${MATERIALPOINTSPATH})
SUBDIRLIST(ELEMENTSDIRS ${ELEMENTSPATH})
SUBDIRLIST(CELLSDIRS ${CELLSPATH})
SUBDIRLIST(CELLELEMENTSDIRS ${CELLELEMENTSPATH})

set(INSTALLED_MODULES "")
set(INSTALLED_MODULE_PATHS "")
message("+------------------------------------------------------------------------------+")

message("checking for modules in 
${COREPATH}. 
"       "If you want to install a Marmot Core module, make sure it is in this directory.")
foreach(DIRNAME ${COREDIRS})
    get_filename_component(MODULEPATH "${COREPATH}/${DIRNAME}" ABSOLUTE) 
    if(EXISTS ${MODULEPATH}/module.cmake)
        if (NOT DEFINED CORE_MODULES)
            set(CORE_MODULES "all")
        endif()
        if( CORE_MODULES STREQUAL "all" OR ${DIRNAME} IN_LIST CORE_MODULES) 
            message("--> found ${DIRNAME}")
            list(APPEND INSTALLED_MODULES ${DIRNAME})
            list(APPEND INSTALLED_MODULEPATHS ${MODULEPATH})
        endif()
    endif()
endforeach(DIRNAME)
message("+------------------------------------------------------------------------------+")

message("Checking for modules in 
${MATERIALSPATH}. 
"       "If you want to install a Marmot Material, make sure it is in this directory.")
foreach(DIRNAME ${MATERIALSDIRS})
    get_filename_component(MODULEPATH "${MATERIALSPATH}/${DIRNAME}" ABSOLUTE) 
    if(EXISTS ${MODULEPATH}/module.cmake)
        if (NOT DEFINED MATERIAL_MODULES)
            set(MATERIAL_MODULES "all")
        endif()
        if(MATERIAL_MODULES STREQUAL "all" OR ${DIRNAME} IN_LIST MATERIAL_MODULES)
            message("--> found ${DIRNAME}")
            list(APPEND INSTALLED_MODULES ${DIRNAME})
            list(APPEND INSTALLED_MODULEPATHS ${MODULEPATH})
        endif()
    endif()
endforeach(DIRNAME)
message("+------------------------------------------------------------------------------+")

message("Checking for modules in 
${ELEMENTSPATH}. 
"       "If you want to install a Marmot Finite Element, make sure it is in this directory.")
foreach(DIRNAME ${ELEMENTSDIRS})
    get_filename_component(MODULEPATH "${ELEMENTSPATH}/${DIRNAME}" ABSOLUTE) 
    if(EXISTS ${MODULEPATH}/module.cmake)
        if (NOT DEFINED ELEMENT_MODULES)
            set(ELEMENT_MODULES "all")
        endif()
        if(ELEMENT_MODULES STREQUAL "all" OR ${DIRNAME} IN_LIST ELEMENT_MODULES)
            message("--> found ${DIRNAME}")
            list(APPEND INSTALLED_MODULES ${DIRNAME})
            list(APPEND INSTALLED_MODULEPATHS ${MODULEPATH})
        endif()
    endif()
endforeach(DIRNAME)
message("+------------------------------------------------------------------------------+")

message("Checking for modules in 
${PARTICLESPATH}. 
"       "If you want to install a Marmot Particle module, make sure it is in this directory.")
foreach(DIRNAME ${PARTICLESDIRS})
    get_filename_component(MODULEPATH "${PARTICLESPATH}/${DIRNAME}" ABSOLUTE) 
    if(EXISTS ${MODULEPATH}/module.cmake)
        if (NOT DEFINED PARTICLE_MODULES)
            set(PARTICLE_MODULES "all")
        endif()
        if(PARTICLE_MODULES STREQUAL "all" OR ${DIRNAME} IN_LIST PARTICLE_MODULES)
            message("--> found ${DIRNAME}")
            list(APPEND INSTALLED_MODULES ${DIRNAME})
            list(APPEND INSTALLED_MODULEPATHS ${MODULEPATH})
        endif()
    endif()
endforeach(DIRNAME)
message("+------------------------------------------------------------------------------+")

message("Checking for modules in 
${MATERIALPOINTSPATH}. 
"       "If you want to install a Marmot Particle module, make sure it is in this directory.")
foreach(DIRNAME ${MATERIALPOINTSDIRS})
    get_filename_component(MODULEPATH "${MATERIALPOINTSPATH}/${DIRNAME}" ABSOLUTE) 
    if(EXISTS ${MODULEPATH}/module.cmake)
        if (NOT DEFINED MATERIALPOINT_MODULES)
            set(MATERIALPOINT_MODULES "all")
        endif()
        if(MATERIALPOINT_MODULES STREQUAL "all" OR ${DIRNAME} IN_LIST MATERIALPOINT_MODULES)
            message("--> found ${DIRNAME}")
            list(APPEND INSTALLED_MODULES ${DIRNAME})
            list(APPEND INSTALLED_MODULEPATHS ${MODULEPATH})
        endif()
    endif()
endforeach(DIRNAME)
message("+------------------------------------------------------------------------------+")

message("Checking for modules in 
${CELLSPATH}. 
"       "If you want to install a Marmot Cell module, make sure it is in this directory.")

foreach(DIRNAME ${CELLSDIRS})
    get_filename_component(MODULEPATH "${CELLSPATH}/${DIRNAME}" ABSOLUTE) 
    if(EXISTS ${MODULEPATH}/module.cmake)
        if (NOT DEFINED CELL_MODULES)
            set(CELL_MODULES "all")
        endif()
        if(CELL_MODULES STREQUAL "all" OR ${DIRNAME} IN_LIST CELL_MODULES)
            message("--> found ${DIRNAME}")
            list(APPEND INSTALLED_MODULES ${DIRNAME})
            list(APPEND INSTALLED_MODULEPATHS ${MODULEPATH})
        endif()
    endif()
endforeach(DIRNAME)
message("+------------------------------------------------------------------------------+")

message("Checking for modules in 
${CELLELEMENTSPATH}. 
"       "If you want to install a Marmot CellElement module, make sure it is in this directory.")

foreach(DIRNAME ${CELLELEMENTSDIRS})
    get_filename_component(MODULEPATH "${CELLELEMENTSPATH}/${DIRNAME}" ABSOLUTE) 
    if(EXISTS ${MODULEPATH}/module.cmake)
        if (NOT DEFINED CELLELEMENT_MODULES)
            set(CELLELEMENT_MODULES "all")
        endif()
        if(CELLELEMENT_MODULES STREQUAL "all" OR ${DIRNAME} IN_LIST CELLELEMENT_MODULES)
            message("--> found ${DIRNAME}")
            list(APPEND INSTALLED_MODULES ${DIRNAME})
            list(APPEND INSTALLED_MODULEPATHS ${MODULEPATH})
        endif()
    endif()
endforeach(DIRNAME)
message("+------------------------------------------------------------------------------+")

message("Now check if all of your modules are listed. If one is missing in this list 
"       "despite being present in the correct directory, also check if the module has a 
"       "proper module.cmake file.")

message("+------------------------------------------------------------------------------+")

foreach(MODULEPATH ${INSTALLED_MODULEPATHS})
    include(${MODULEPATH}/module.cmake)
endforeach(MODULEPATH)

message("+------------------------------------------------------------------------------+")
message(STATUS "CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")

# compile the shared library
add_library(${PROJECT_NAME} SHARED ${sources})
target_include_directories(${PROJECT_NAME} PRIVATE include)
target_link_libraries (${PROJECT_NAME} Eigen3::Eigen)

if (SHARED_LIBRARIES)
    message("+------------------------------------------------------------------------------+")
    message("Marmot modules have requested to link against following shared libraries:")
    foreach(SHARED_LIBRARY ${SHARED_LIBRARIES})
        message("--> ${SHARED_LIBRARY}")
        target_link_libraries (${PROJECT_NAME} ${SHARED_LIBRARY})
    endforeach()
    message("+------------------------------------------------------------------------------+")
endif()

foreach(MODULEPATH ${INSTALLED_MODULEPATHS})
    if(EXISTS ${MODULEPATH}/test.cmake)
        include(${MODULEPATH}/test.cmake)
    endif()
endforeach(MODULEPATH)

# and finally install the library to the system location
set_target_properties(${PROJECT_NAME} PROPERTIES
    PUBLIC_HEADER "${publicheaders}"
)


install(TARGETS ${PROJECT_NAME}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/Marmot
    )
