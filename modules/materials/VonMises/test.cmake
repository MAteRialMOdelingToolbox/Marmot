# add current directory to the source for the tests
SET(CURR_TEST_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/test")

# Tests for VonMises material
add_marmot_test("TestVonMises" "${CURR_TEST_SOURCE_DIR}/test.cpp")

# a single large strain increment, VonMises and ADVonMises: the return mapping converges on the relative tolerance
add_marmot_test("TestVonMisesLargeIncrement" "${CURR_TEST_SOURCE_DIR}/testLargeIncrement.cpp" REQUIRES ADVonMises)
