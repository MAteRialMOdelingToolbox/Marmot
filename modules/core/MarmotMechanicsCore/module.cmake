list(APPEND INSTALLED_MODULE_INCLUDE_DIRS "${CMAKE_CURRENT_LIST_DIR}/include")
file(GLOB module_sources CONFIGURE_DEPENDS "${CMAKE_CURRENT_LIST_DIR}/src/*.cpp")
list(APPEND sources ${module_sources})
list(APPEND publicheaders
    "${CMAKE_CURRENT_LIST_DIR}/include/Marmot/MarmotMaterialHypoElastic.h"
    "${CMAKE_CURRENT_LIST_DIR}/include/Marmot/MarmotMaterialHypoElasticFactory.h"
    "${CMAKE_CURRENT_LIST_DIR}/include/Marmot/MarmotVoigt.h"
    "${CMAKE_CURRENT_LIST_DIR}/include/Marmot/MarmotMaterialPointSolverHypoElastic.h"
    )
