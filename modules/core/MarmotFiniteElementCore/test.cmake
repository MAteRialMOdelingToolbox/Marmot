# add current directory to the source for the tests
SET(CURR_TEST_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/test")

# Tests for MarmotFiniteElementBoundary
add_marmot_test("TestMarmotFiniteElementBoundary" "${CURR_TEST_SOURCE_DIR}/TestMarmotFiniteElementBoundary.cpp")
