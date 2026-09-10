# add current directory to the source for the tests
SET(CURR_TEST_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/test")

# Tests for MarmotFiniteElementBoundary
add_marmot_test("TestMarmotFiniteElementBoundary" "${CURR_TEST_SOURCE_DIR}/TestMarmotFiniteElementBoundary.cpp")

# Tests for MarmotFiniteElementSpatialWrapper
add_marmot_test("TestMarmotFiniteElementSpatialWrapper" "${CURR_TEST_SOURCE_DIR}/TestMarmotFiniteElementSpatialWrapper.cpp")

# Tests for MarmotDofLayoutTools
add_marmot_test("TestMarmotDofLayoutTools" "${CURR_TEST_SOURCE_DIR}/TestMarmotDofLayoutTools.cpp")

# Tests for MarmotEnhancedAssumedStrain
add_marmot_test("TestMarmotEnhancedAssumedStrain" "${CURR_TEST_SOURCE_DIR}/TestMarmotEnhancedAssumedStrain.cpp")
