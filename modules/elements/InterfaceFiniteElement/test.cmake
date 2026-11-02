# add current directory to the source for the tests
SET(CURR_TEST_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/test")

# Tests for InterfaceFiniteElement
add_marmot_test("TestInterfaceFiniteElement" "${CURR_TEST_SOURCE_DIR}/TestInterfaceFiniteElement.cpp")
