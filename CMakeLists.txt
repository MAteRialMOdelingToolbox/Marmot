include_directories(${CMAKE_CURRENT_LIST_DIR}/include)
file(GLOB sources_element "${CMAKE_CURRENT_LIST_DIR}/src/*.cpp")
list(APPEND sources ${sources_element})
