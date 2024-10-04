include_directories(${CMAKE_CURRENT_LIST_DIR}/include)
file(GLOB sources_material "${CMAKE_CURRENT_LIST_DIR}/src/*.cpp")
list(APPEND sources ${sources_material})

list(APPEND publicheaders 
    "${CMAKE_CURRENT_LIST_DIR}/include/Marmot/MarmotMath.h" 
    "${CMAKE_CURRENT_LIST_DIR}/include/Marmot/MarmotConstants.h" 
    "${CMAKE_CURRENT_LIST_DIR}/include/Marmot/MarmotTensor.h" 
    "${CMAKE_CURRENT_LIST_DIR}/include/Marmot/MarmotTypedefs.h" 
    )
