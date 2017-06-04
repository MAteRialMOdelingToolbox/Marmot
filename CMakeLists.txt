# matthias neuner, 2016, standard static lib cmake file for BFT
# 
# depends on EIGEN library 3.2.7(+) 

cmake_minimum_required(VERSION 3.1)  

project(bftUserLibrary)

include_directories(include)
#include_directories(/path/to/additional/headers) <-- extend if necessary...

file(GLOB sources "src/*.cpp")
file(GLOB headers "include/*.h")

set(CMAKE_BUILD_TYPE Release)
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
    HoekBrown
    ModLeon
    ModLeonNonLocal
    Meschke
    SchaedlichSchweiger
    ShotLeon
    SolidificationCreep
    UntereggerRockMass
    UntereggerRockMassNonLocal

    uelCPS4
    uelCPS4NonLocal
    uelCPS8R
    uelCPS8RNonLocal

    uelC3D8
    uelC3D8NonLocal
    
    uelCPE4
    uelCPE4NonLocal
    uelCPE8R
    uelCPE8RNonLocal
    )

#optional libraries
foreach(library ${optionalLibraries})
    file(TO_CMAKE_PATH "${optionalLibsRootDir}/${library}" librarypath)
    if(EXISTS ${librarypath})
        message("--> found ${library}")
        include_directories(${librarypath}/include)
        add_definitions(-D${library})
        set(FOUND_LIB "FOUND_LIB-NOTFOUND")
        find_library(FOUND_LIB ${library} ${librarypath}/lib)
        if(NOT FOUND_LIB)
            message(FATAL_ERROR "found directory of ${library}, but no binary")
        endif()
        target_link_libraries(${PROJECT_NAME} ${FOUND_LIB})
    endif()
endforeach(library)

# mandatory bftMechanics
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
