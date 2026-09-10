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

# Tests for MarmotFiniteElement2D
add_marmot_test("TestMarmotFiniteElement2D" "${CURR_TEST_SOURCE_DIR}/TestMarmotFiniteElement2D.cpp")

# Tests for MarmotGeometryElement
add_marmot_test("TestMarmotGeometryElement" "${CURR_TEST_SOURCE_DIR}/TestMarmotGeometryElement.cpp")

# Tests for MarmotFiniteElementBasic
add_marmot_test("TestMarmotFiniteElementBasic" "${CURR_TEST_SOURCE_DIR}/TestMarmotFiniteElementBasic.cpp")
