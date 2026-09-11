list(APPEND INSTALLED_MODULE_INCLUDE_DIRS "${CMAKE_CURRENT_LIST_DIR}/include")
file(GLOB module_sources CONFIGURE_DEPENDS "${CMAKE_CURRENT_LIST_DIR}/src/*.cpp")
list(APPEND sources ${module_sources})
list(APPEND publicheaders
    "${CMAKE_CURRENT_LIST_DIR}/include/Marmot/MarmotMaterialFiniteStrain.h"
    "${CMAKE_CURRENT_LIST_DIR}/include/Marmot/MarmotMaterialFiniteStrainFactory.h"
    "${CMAKE_CURRENT_LIST_DIR}/include/Marmot/MarmotMaterialPointSolverFiniteStrain.h"
    "${CMAKE_CURRENT_LIST_DIR}/include/Marmot/MarmotMaterialHughesWinget.h"
    )
