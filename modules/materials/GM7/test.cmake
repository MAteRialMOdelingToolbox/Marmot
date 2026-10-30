# add current directory to the source for the tests
SET(CURR_TEST_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/test")

# Tests for GM7 material
add_marmot_test("TestGM7" "${CURR_TEST_SOURCE_DIR}/testGM7.cpp" REQUIRES MicroplaneM7)
