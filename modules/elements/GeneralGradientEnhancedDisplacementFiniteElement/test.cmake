# add current directory to the source for the tests
SET(CURR_TEST_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/test")

# Tests for GeneralGradientEnhancedDisplacementFiniteElement
add_marmot_test("TestGeneralGradientEnhancedDisplacementFiniteElement" "${CURR_TEST_SOURCE_DIR}/TestGeneralGradientEnhancedDisplacementFiniteElement.cpp")
