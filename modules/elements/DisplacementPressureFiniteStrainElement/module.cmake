set(MODULE_NAME "DisplacementPressureFiniteStrainElement")

set(MODULE_DEPENDENCIES
    MarmotMechanicsCore
    MarmotFiniteElementCore
    MarmotFiniteStrainMechanicsCore
    )

set(DEPENDENCIES_FULFILLED TRUE)
foreach(DEPENDENCY IN LISTS MODULE_DEPENDENCIES)
    if(NOT DEPENDENCY IN_LIST INSTALLED_MODULES)
        message("----> module ${MODULE_NAME} dependency not fulfilled: ${DEPENDENCY}")
        set(DEPENDENCIES_FULFILLED FALSE)
    endif()
endforeach()

if(DEPENDENCIES_FULFILLED)
    list(APPEND INSTALLED_MODULE_INCLUDE_DIRS "${CMAKE_CURRENT_LIST_DIR}/include")
    file(GLOB module_sources CONFIGURE_DEPENDS "${CMAKE_CURRENT_LIST_DIR}/src/*.cpp")
    list(APPEND sources ${module_sources})
endif()
