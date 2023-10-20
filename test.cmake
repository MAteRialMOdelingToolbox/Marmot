add_executable("MarmotMathCore_executable" "${CMAKE_CURRENT_LIST_DIR}/test/test.cpp")
target_link_libraries("MarmotMathCore_executable" Marmot)
add_test(NAME "MarmotMathCore" COMMAND "MarmotMathCore_executable")
