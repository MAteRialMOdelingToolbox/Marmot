# add current directory to the source for the tests
SET(CURR_TEST_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/test")

# Tests for MicroplaneM7 material
add_marmot_test("TestMicroplaneM7" "${CURR_TEST_SOURCE_DIR}/testMicroplaneM7.cpp")
