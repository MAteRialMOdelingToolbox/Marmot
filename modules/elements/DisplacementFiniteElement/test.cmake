# add current directory to the source for the tests
SET(CURR_TEST_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/test")

# Tests for DisplacementFiniteElement
add_marmot_test("TestDisplacementFiniteElement" "${CURR_TEST_SOURCE_DIR}/TestDisplacementFiniteElement.cpp")
