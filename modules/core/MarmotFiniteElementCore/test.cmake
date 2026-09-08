# add current directory to the source for the tests
SET(CURR_TEST_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/test")

# Tests for MarmotBulkViscosity
add_marmot_test("TestMarmotBulkViscosity" "${CURR_TEST_SOURCE_DIR}/TestMarmotBulkViscosity.cpp")
