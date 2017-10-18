# matthias neuner, 2016, standard static lib cmake file for BFT
# 
# depends on EIGEN library 3.2.7(+) 

cmake_minimum_required(VERSION 3.1)  
project(bftUserLibrary)

include_directories(include)

file(GLOB sources "src/*.cpp")
file(GLOB headers "include/*.h")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall ")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_library(${PROJECT_NAME} SHARED ${sources})

set(bftMechanicsDir "../bftMechanics/")

set(optionalLibsRootDir "../")
set (optionalLibraries 
    LinearElastic
    LinearElasticSolidificationCreep
    HoekBrown
    ModLeon
    ModLeonSemiExplicit
    ModLeonAdaptive
    ModLeonSemiExplicitAdaptive
    ModLeonPlaneStress
    ModLeonNonLocal
    Meschke
    SchaedlichSchweiger
    ShotLeon
    ShotLeonV2
    ShotLeonNonLocal
    ShotLeonV2NonLocal
    SolidificationCreep
    UntereggerRockMass
    UntereggerRockMassNonLocal

    uelDisplacement
    uelDisplacementEAS
    uelNonLocal

    )

#optional libraries
foreach(library ${optionalLibraries})
    file(TO_CMAKE_PATH "${optionalLibsRootDir}/${library}" librarypath)
    get_filename_component(librarypath ${librarypath} ABSOLUTE) # make the path absolute, which is is sometimes necessary for some reason
    if(EXISTS ${librarypath})
        message("--> found ${library}")
        include_directories(${librarypath}/include)
        string(TOUPPER ${library} LIBRARY)
        add_definitions(-D${LIBRARY})
        set(FOUND_LIB "FOUND_LIB-NOTFOUND")
        find_library(FOUND_LIB ${library} ${librarypath}/lib)
        if(NOT FOUND_LIB)
            message(FATAL_ERROR "found directory of ${library}, but no binary")
        endif()
        target_link_libraries(${PROJECT_NAME} ${FOUND_LIB})
    endif()
endforeach(library)

# mandatory bftMechanics
get_filename_component(bftMechanicsDir ${bftMechanicsDir} ABSOLUTE) # make the path absolute, which is is sometimes necessary for some reason
if(EXISTS ${bftMechanicsDir})
    message("--> found bftMechanics")
    include_directories(${bftMechanicsDir}/include)
    find_library(bftMechanicsLib bftMechanics ${bftMechanicsDir}/lib)
    target_link_libraries(${PROJECT_NAME} ${bftMechanicsLib})
else()
    message(FATAL_ERROR "could not find bftMechanics")
endif()

if(CMAKE_CONFIGURATION_TYPES)
   set(CMAKE_CONFIGURATION_TYPES Release)
endif()

install(TARGETS ${PROJECT_NAME} DESTINATION lib)
