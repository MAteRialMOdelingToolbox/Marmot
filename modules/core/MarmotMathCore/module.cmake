list(APPEND INSTALLED_MODULE_INCLUDE_DIRS "${CMAKE_CURRENT_LIST_DIR}/include")
file(GLOB module_sources CONFIGURE_DEPENDS "${CMAKE_CURRENT_LIST_DIR}/src/*.cpp")
list(APPEND sources ${module_sources})

list(APPEND publicheaders
    "${CMAKE_CURRENT_LIST_DIR}/include/Marmot/MarmotMath.h"
    "${CMAKE_CURRENT_LIST_DIR}/include/Marmot/MarmotConstants.h"
    "${CMAKE_CURRENT_LIST_DIR}/include/Marmot/MarmotFastorTensorBasics.h"
    "${CMAKE_CURRENT_LIST_DIR}/include/Marmot/MarmotNumericalDifferentiation.h"
    "${CMAKE_CURRENT_LIST_DIR}/include/Marmot/MarmotTensor.h"
    "${CMAKE_CURRENT_LIST_DIR}/include/Marmot/MarmotTypedefs.h"
    )
