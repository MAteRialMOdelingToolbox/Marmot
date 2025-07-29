# add current directory to the source for the tests
SET(CURR_TEST_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/test")

# Tests for MarmotDeformationMeasures
add_marmot_test("TestMarmotDeformationMeasures" "${CURR_TEST_SOURCE_DIR}/TestMarmotDeformationMeasures.cpp")

# Tests for MarmotEnergyDensityFunctions
add_marmot_test("TestMarmotEnergyDensityFunctions" "${CURR_TEST_SOURCE_DIR}/TestMarmotEnergyDensityFunctions.cpp")

# Tests for MarmotFiniteStrainPlasticity
add_marmot_test("TestMarmotFiniteStrainPlasticity" "${CURR_TEST_SOURCE_DIR}/TestMarmotFiniteStrainPlasticity.cpp")

# Tests for MarmotStressMeasures
add_marmot_test("TestMarmotStressMeasures" "${CURR_TEST_SOURCE_DIR}/TestMarmotStressMeasures.cpp")
