# add current directory to the source for the tests
SET(CURR_TEST_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/test")

# Tests for the ModelCode2010BondSlip law
add_marmot_test("TestModelCode2010BondSlip" "${CURR_TEST_SOURCE_DIR}/test.cpp")
