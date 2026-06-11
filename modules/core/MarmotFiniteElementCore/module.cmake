list(APPEND INSTALLED_MODULE_INCLUDE_DIRS "${CMAKE_CURRENT_LIST_DIR}/include")
file(GLOB module_sources CONFIGURE_DEPENDS "${CMAKE_CURRENT_LIST_DIR}/src/*.cpp")
list(APPEND sources ${module_sources})
list(APPEND publicheaders
    "${CMAKE_CURRENT_LIST_DIR}/include/Marmot/MarmotElement.h"
    "${CMAKE_CURRENT_LIST_DIR}/include/Marmot/MarmotElementFactory.h"
    "${CMAKE_CURRENT_LIST_DIR}/include/Marmot/MarmotElementProperty.h"
    )
