include_directories(${CMAKE_CURRENT_LIST_DIR}/include)
file(GLOB sources_material "${CMAKE_CURRENT_LIST_DIR}/src/*.cpp")
list(APPEND sources ${sources_material})
list(APPEND publicheaders 
    "${CMAKE_CURRENT_LIST_DIR}/include/Marmot/MarmotMaterialFiniteStrain.h" 
    )
