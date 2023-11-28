set(MODULE_NAME
    "DisplacementFinitStrainULElement")

set(MODULES_DEPENDENCIES
    MarmotMechanicsCore
    MarmotFiniteElementCore 
    MarmotFiniteStrainMechanicsCore 
    )

set(DEPENDECIES_FULLFILLED TRUE)
foreach( DEPENDENCY ${MODULES_DEPENDENCIES} )
    if (NOT DEPENDENCY IN_LIST INSTALLED_MODULES)
        message("----> " "module ${MODULE_NAME} dependency not fulfilled: ${DEPENDENCY}")
        set(DEPENDECIES_FULLFILLED FALSE)
    endif()
endforeach()

if ( DEPENDECIES_FULLFILLED )
    set(CMAKE_CXX_STANDARD 20)
    include_directories(${CMAKE_CURRENT_LIST_DIR}/include)
    file(GLOB sources_material "${CMAKE_CURRENT_LIST_DIR}/src/*.cpp")
    list(APPEND sources ${sources_material})
endif()
