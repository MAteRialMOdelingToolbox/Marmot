# add current directory to the source for the tests
SET(CURR_TEST_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/test")

# Tests for the ModelCode2010CoulombBondSlip law
add_marmot_test("TestModelCode2010CoulombBondSlip" "${CURR_TEST_SOURCE_DIR}/test.cpp" REQUIRES ModelCode2010BondSlip)
