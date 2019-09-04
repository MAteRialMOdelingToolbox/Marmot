# matthias neuner, 2016, 
# 
# depends on EIGEN library 3.2.7(+) 

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


if(CMAKE_VERSION VERSION_LESS "3.9.0") 
    set(outdatedcmake TRUE)
    message("--> Please consider to switch to CMake 3.9.0")
endif()

set(bftMechanicsDir "modules/bftMechanics/")
set(bftFiniteElementCoreDir "modules/bftFiniteElementCore/")

set(optionalLibsRootDir "modules/")

set (materials
    Barodesy
    BarodesyGradientVoid
    LinearElastic
    LinearElasticNonLocal
    LinearElasticSolidificationCreep
    HoekBrown
    ModLeon
    ModLeonSemiExplicit
    ModLeonAnalytical
    ModLeonAdaptive
    ModLeonSemiExplicitAdaptive
    ModLeonPlaneStress
    ModLeonNonLocal
    MohrCoulomb
    Meschke
    SchaedlichSchweiger
    ShotLeon
    ShotLeonV2
    ShotLeonNonLocal
    ShotLeonV2NonLocal
    SolidificationCreep
    StVenantKirchhoffIsotropic
    UntereggerRockMass
    UntereggerRockMassNonLocal
    UntereggerRockMassPlaxis)

set(elements
    uelDisplacement
    uelDisplacementTL
    uelDisplacementUL
    uelDisplacementEAS
    uelNonLocal
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
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic")
    # only activate if you know what you are doing:
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
elseif(MSVC)
    if(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
        string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
    endif()
endif()

# If not build type was specified, Release is default
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
    message("--> setting to default CMAKE_BUILD_TYPE Release")
endif()

# Add sources
file(GLOB sources "src/*.cpp" "src/*.f90")
#message(${GLOB})
file(GLOB headers "include/*.h")

#optional materials 
foreach(library ${materials})
    file(TO_CMAKE_PATH "${optionalLibsRootDir}/materials/${library}" librarypath)
    get_filename_component(librarypath ${librarypath} ABSOLUTE) 
    if(EXISTS ${librarypath})
        message("--> found ${library}")
        include_directories(${librarypath}/include)
        file(GLOB libsources "${librarypath}/src/*.cpp" "${librarypath}/src/*.f90")
        set(sources ${sources} ${libsources})
        string(TOUPPER ${library} LIBRARY)
        add_definitions(-D${LIBRARY})
    endif()
endforeach(library)

#optional elements 
foreach(library ${elements})
    file(TO_CMAKE_PATH "${optionalLibsRootDir}/elements/${library}" librarypath)
    get_filename_component(librarypath ${librarypath} ABSOLUTE) 
    if(EXISTS ${librarypath})
        message("--> found ${library}")
        set(bftElements_are_present ON)
        include_directories(${librarypath}/include)
        file(GLOB libsources "${librarypath}/src/*.cpp")
        set(sources ${sources} ${libsources})
        string(TOUPPER ${library} LIBRARY)
        add_definitions(-D${LIBRARY})
    endif()
endforeach(library)

#  bftFiniteElementCore, mandatory if bftElements are used 
get_filename_component(bftFiniteElementCoreDir ${bftFiniteElementCoreDir} ABSOLUTE) 
if(EXISTS ${bftFiniteElementCoreDir})
    message("--> found bftFiniteElementCore")
    include_directories(${bftFiniteElementCoreDir}/include)
    file(GLOB libsources "${bftFiniteElementCoreDir}/src/*.cpp")
    set(sources ${sources} ${libsources})
else()
    if (DEFINED bftElements_are_present)
        message(FATAL_ERROR "found bftElements, but no bftFiniteElementCore")
    endif()
endif()

# mandatory bftMechanics
get_filename_component(bftMechanicsDir ${bftMechanicsDir} ABSOLUTE) 
if(EXISTS ${bftMechanicsDir})
    message("--> found bftMechanics")
    include_directories(${bftMechanicsDir}/include)
    file(GLOB libsources "${bftMechanicsDir}/src/*.cpp")
    set(sources ${sources} ${libsources})
else()
    message(FATAL_ERROR "--> could not find bftMechanics")
endif()

add_library(${PROJECT_NAME} SHARED ${sources})
install(TARGETS ${PROJECT_NAME} DESTINATION lib)
