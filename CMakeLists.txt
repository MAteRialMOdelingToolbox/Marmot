# Copyright (C) 2019 Matthias Neuner <matthias.neuner@uibk.ac.at>

# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
cmake_minimum_required(VERSION 3.1)  

project(bftUserLibrary CXX Fortran)

include_directories(include)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_FIND_LIBRARY_PREFIXES lib)

## avoid in-source compilation
set(DEFAULT_OUT_OF_SOURCE_FOLDER "build")
if (${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds are not allowed. Create and/or change to ${DEFAULT_OUT_OF_SOURCE_FOLDER} \
 folder and run then cmake from there, e.g. cd build && cmake ../ . You first may need to remove CMakeCache.txt and CMakeFiles/*.")
endif ()

# Disable (shared memory) parallelization within bftUserLibrary - such a parallelization is usually not desired
# add_compile_definitions(EIGEN_DONT_PARALLELIZE)                     # for Eigen
# add_compile_definitions(BLAZE_USE_SHARED_MEMORY_PARALLELIZATION=0)  # for Blaze

set(bftMechanicsDir "modules/bftMechanics/")
set(optionalLibsRootDir "modules/")

set(optionalCoreLibraries
    bftFiniteElementCore
    bftCosseratCore
    bftMicromorphicCore
    )

set (materials
    Barodesy
    BarodesyGradientVoid
    CDP
    CDPFibreReinforced
    CDPM2
    BarodesyGradientDeformationModulus
    CosseratLinearElastic
    CosseratDruckerPrager
    DiffuseMicroCracking
    DruckerPrager 
    DruckerPragerMD
    GMCDP
    GMBiotElastic
    GMDruckerPrager
    GMNeoHooke
    HoekBrown
    JointedRock
    LinearElastic
    transverseisotropiclinearelastic
    LinearElasticNonLocal
    LinearElasticSolidificationCreep
    MCDP
    ModifiedCamClay
    ModLeon
    ModLeonSemiExplicit
    ModLeonAnalytical
    ModLeonAdaptive
    ModLeonSemiExplicitAdaptive
    ModLeonPlaneStress
    ModLeonNonLocal
    MohrCoulomb
    MisesTI
    ViscoPlasticShotcreteModel 
    PorousElastic
    RockDamagePlasticity
    RockDamagePlasticityNonLocal
    SandHypo
    SandHypoMicropolar
    SchaedlichSchweiger
    ShotLeon
    ShotLeonV2
    ShotLeonNonLocal
    ShotLeonV2NonLocal
    SolidificationCreep
    StVenantKirchhoffIsotropic
    UntereggerRockMassPlaxis
    )

set(elements
    uelCosserat
    uelDisplacement
    uelDisplacementTL
    uelDisplacementUL
    uelDisplacementEAS
    uelNonLocal
    uelNonLocalCosserat
    uelGradientEnhancedMicropolar
    uelGradientEnhancedMicropolarUL
    uelNonLocalUL
    uelNonLocalULFBar
    uelNonLocalEAS
    uelNonLocalMixed
    )

# Compiler specific stettings
if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-parameter")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-long-long")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-copy")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic")
    # only activate if you know what you are doing:
    # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
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

# Add sources
file(GLOB sources "src/*.cpp")
#message(${GLOB})
file(GLOB headers "include/*.h")

# mandatory bftMechanics
get_filename_component(bftMechanicsDir ${bftMechanicsDir} ABSOLUTE) 
if(EXISTS ${bftMechanicsDir})
    message("--> found bftMechanics")
    include(${bftMechanicsDir}/module.cmake)
else()
    message(FATAL_ERROR "--> could not find bftMechanics")
endif()

#
#optional core libraries
#
foreach(library ${optionalCoreLibraries})
    file(TO_CMAKE_PATH "${optionalLibsRootDir}/${library}" librarypath)
    get_filename_component(librarypath ${librarypath} ABSOLUTE) 
    if(EXISTS ${librarypath})
        message("--> found ${library}")
        include(${librarypath}/module.cmake)
    endif()
endforeach(library)

#
#optional materials 
#
foreach(library ${materials})
    file(TO_CMAKE_PATH "${optionalLibsRootDir}/materials/${library}" librarypath)
    get_filename_component(librarypath ${librarypath} ABSOLUTE) 
    if(EXISTS ${librarypath})
        message("--> found ${library}")
        include(${librarypath}/module.cmake)
    endif()
endforeach(library)

#
#optional elements 
#
foreach(library ${elements})
    file(TO_CMAKE_PATH "${optionalLibsRootDir}/elements/${library}" librarypath)
    get_filename_component(librarypath ${librarypath} ABSOLUTE) 
    if(EXISTS ${librarypath})
        message("--> found ${library}")
        include(${librarypath}/module.cmake)
    endif()
endforeach(library)

add_library(${PROJECT_NAME} SHARED ${sources})
install(TARGETS ${PROJECT_NAME} DESTINATION lib)
