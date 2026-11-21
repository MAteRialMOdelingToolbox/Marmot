# add current directory to the source for the tests
SET(CURR_TEST_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/test")

# Tests for GM7 material
add_marmot_test("TestGM7" "${CURR_TEST_SOURCE_DIR}/testGM7.cpp" REQUIRES MicroplaneM7)

# Tests for GM7 in finite strain through the corotational wrapper (GM7/COROTATIONAL)
add_marmot_test("TestGM7Corotational" "${CURR_TEST_SOURCE_DIR}/testGM7Corotational.cpp" REQUIRES MicroplaneM7)
