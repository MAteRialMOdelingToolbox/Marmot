# add current directory to the source for the tests
SET(CURR_TEST_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/test")

# Tests for CompressibleFiniteStrainFractionalSchapery material
add_marmot_test("TestCompressibleFiniteStrainFractionalSchapery" "${CURR_TEST_SOURCE_DIR}/test.cpp")
