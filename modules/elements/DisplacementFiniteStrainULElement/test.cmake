# add current directory to the source for the tests
SET(CURR_TEST_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/test")

# Tests for DisplacementFiniteStrainULElement
add_marmot_test("TestDisplacementFiniteStrainULElement" "${CURR_TEST_SOURCE_DIR}/TestDisplacementFiniteStrainULElement.cpp" REQUIRES CompressibleHyperelasticity)
add_marmot_test("TestDisplacementFiniteStrainULElementSmallStrainLimit" "${CURR_TEST_SOURCE_DIR}/TestDisplacementFiniteStrainULElementSmallStrainLimit.cpp" REQUIRES CompressibleHyperelasticity DisplacementFiniteElement LinearElastic)
