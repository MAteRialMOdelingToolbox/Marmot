# Copyright (C) 2019 Matthias Neuner <matthias.neuner@uibk.ac.at>

# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

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

# Disable (shared memory) parallelization within Marmot - such a parallelization is usually not desired
# add_compile_definitions(EIGEN_DONT_PARALLELIZE)                     # for Eigen
# add_compile_definitions(BLAZE_USE_SHARED_MEMORY_PARALLELIZATION=0)  # for Blaze

## find Eigen library
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
# include_directories(${EIGEN3_INCLUDE_DIR})
message("--> found Eigen: ${EIGEN3_INCLUDE_DIR}")

set(MODULES_DIR "modules/")

if(NOT CORE_MODULES)
    set(CORE_MODULES
        MarmotMechanicsCore
        MarmotFiniteElementCore
        MarmotCosseratCore
        MarmotMicromorphicCore
        )
endif()

if(NOT MATERIAL_MODULES)
    set (MATERIAL_MODULES
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
        GMCDPFiniteStrain
        GMBiotElastic
        GMDruckerPrager
        GMNeoHooke
        GradientEnhancedDruckerPrager
        HoekBrown
        JointedRock
        JointedHoekBrown
        LinearElastic
        transverseisotropiclinearelastic
        OrthotropicLinearElastic
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
        ViscoPlasticShotcreteModel 
        PorousElastic
        RockDamagePlasticity
        RockDamagePlasticityNonLocal
        ORDP
        SandHypo
        SandHypoMicropolar
        SCDPV3NonLocal
        SchaedlichSchweiger
        ShotLeon
        ShotLeonV2
        ShotLeonNonLocal
        ShotLeonV2NonLocal
        ShotcreteSolidificationCreep
        ShotcreteHydration
        SolidificationCreep
        StVenantKirchhoffIsotropic
        UntereggerRockMassPlaxis
        )
endif()

if(NOT ELEMENT_MODULES)
    set(ELEMENT_MODULES
        CosseratFiniteElement
        DisplacementFiniteElement
        DisplacementTLFiniteElement
        DisplacementULFiniteElement
        DisplacementEASFiniteElement
        GradientEnhancedDisplacementFiniteElement
        GradientEnhancedDisplacementULFiniteElement
        GradientEnhancedDisplacementULFBarFiniteElement
        GradientEnhancedDisplacementEASFiniteElement
        GradientEnhancedDisplacementMixedFiniteElement
        GradientEnhancedCosseratFiniteElement
        GradientEnhancedMicropolarFiniteElement
        GradientEnhancedMicropolarULFiniteElement
        )
endif()

#
# Add sources
#
file(GLOB sources "src/*.cpp")
file(GLOB publicheaders "include/Marmot/*.h")

# optional core modules 
foreach(library ${CORE_MODULES})
    file(TO_CMAKE_PATH "${MODULES_DIR}/${library}" librarypath)
    get_filename_component(librarypath ${librarypath} ABSOLUTE) 
    if(EXISTS ${librarypath}/module.cmake)
        message("--> found ${library}")
        include(${librarypath}/module.cmake)
    endif()
endforeach(library)

#optional material modules
foreach(library ${MATERIAL_MODULES})
    file(TO_CMAKE_PATH "${MODULES_DIR}/materials/${library}" librarypath)
    get_filename_component(librarypath ${librarypath} ABSOLUTE) 
    if(EXISTS ${librarypath}/module.cmake)
        message("--> found ${library}")
        include(${librarypath}/module.cmake)
    endif()
endforeach(library)


#optional elements modules
foreach(library ${ELEMENT_MODULES})
    file(TO_CMAKE_PATH "${MODULES_DIR}/elements/${library}" librarypath)
    get_filename_component(librarypath ${librarypath} ABSOLUTE) 
    if(EXISTS ${librarypath}/module.cmake)
        message("--> found ${library}")
        include(${librarypath}/module.cmake)
    endif()
endforeach(library)

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
