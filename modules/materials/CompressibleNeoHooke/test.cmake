# add current directory to the source for the tests
SET(CURR_TEST_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/test")

# Tests for VonMises material
add_marmot_test("TestCompressibleNeoHooke" "${CURR_TEST_SOURCE_DIR}/test.cpp")
