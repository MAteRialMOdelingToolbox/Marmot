# add current directory to the source for the tests
SET(CURR_TEST_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/test")

# Tests for MarmotJournal
add_marmot_test("TestMarmotJournal" "${CURR_TEST_SOURCE_DIR}/TestMarmotJournal.cpp")

# Tests for MarmotTesting
add_marmot_test("TestMarmotTesting" "${CURR_TEST_SOURCE_DIR}/TestMarmotTesting.cpp")
