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
cmake_minimum_required(VERSION 3.3)  

project(Marmot CXX Fortran)

# avoid in-source compilation
if (${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds are not allowed. Create and/or change to your build folder \
  and run then cmake from there, e.g. mkdir build && cd build && cmake ../ . \\n
  You first may need to remove CMakeCache.txt and CMakeFiles/*.")
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_FIND_LIBRARY_PREFIXES lib)

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
# include_directories(${EIGEN3_INCLUDE_DIR})
message("--> found Eigen: ${EIGEN3_INCLUDE_DIR}")

## find autodiff library (https://autodiff.github.io/)
find_path(AUTODIFF_INCLUDE_DIR autodiff REQUIRED NO_MODULE)
include_directories(${AUTODIFF_INCLUDE_DIR})
message("--> found autodiff: ${AUTODIFF_INCLUDE_DIR}")

set(MODULES_DIR ${CMAKE_SOURCE_DIR}/modules)
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
get_filename_component(ELEMENTSPATH ${MODULES_DIR}/elements ABSOLUTE) 

SUBDIRLIST(COREDIRS ${COREPATH})
SUBDIRLIST(MATERIALSDIRS ${MATERIALSPATH})
SUBDIRLIST(ELEMENTSDIRS ${ELEMENTSPATH})

set(INSTALLED_MODULES "")

foreach(DIRNAME ${COREDIRS})
    get_filename_component(MODULEPATH "${COREPATH}/${DIRNAME}" ABSOLUTE) 
    if(EXISTS ${MODULEPATH}/module.cmake)
        if (NOT DEFINED CORE_MODULES)
            set(CORE_MODULES "all")
        endif()
        if( ${CORE_MODULES} STREQUAL "all" OR ${DIRNAME} IN_LIST CORE_MODULES) 
            message("--> found ${DIRNAME}")
            include(${MODULEPATH}/module.cmake)
            list(APPEND INSTALLED_MODULES ${DIRNAME})
        endif()
    endif()
endforeach(DIRNAME)

foreach(DIRNAME ${MATERIALSDIRS})
    get_filename_component(MODULEPATH "${MATERIALSPATH}/${DIRNAME}" ABSOLUTE) 
    if(EXISTS ${MODULEPATH}/module.cmake)
        if (NOT DEFINED MATERIAL_MODULES)
            set(MATERIAL_MODULES "all")
        endif()
        if(${MATERIAL_MODULES} STREQUAL "all" OR ${DIRNAME} IN_LIST MATERIAL_MODULES)
            message("--> found ${DIRNAME}")
            include(${MODULEPATH}/module.cmake)
            list(APPEND INSTALLED_MODULES ${DIRNAME})
        endif()
    endif()
endforeach(DIRNAME)

foreach(DIRNAME ${ELEMENTSDIRS})
    get_filename_component(MODULEPATH "${ELEMENTSPATH}/${DIRNAME}" ABSOLUTE) 
    if(EXISTS ${MODULEPATH}/module.cmake)
        if (NOT DEFINED ELEMENT_MODULES)
            set(ELEMENT_MODULES "all")
        endif()
        if(${ELEMENT_MODULES} STREQUAL "all" OR ${DIRNAME} IN_LIST ELEMENT_MODULES)
            message("--> found ${DIRNAME}")
            include(${MODULEPATH}/module.cmake)
            list(APPEND INSTALLED_MODULES ${DIRNAME})
        endif()
    endif()
endforeach(DIRNAME)

message(STATUS "CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")

# compile the shared library
add_library(${PROJECT_NAME} SHARED ${sources})
target_include_directories(${PROJECT_NAME} PRIVATE include)
target_link_libraries (${PROJECT_NAME} Eigen3::Eigen)

# and finally install the library to the system location
set_target_properties(${PROJECT_NAME} PROPERTIES
    PUBLIC_HEADER "${publicheaders}"
)
install(TARGETS ${PROJECT_NAME}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/Marmot
    )
