# add current directory to the source for the tests
SET(CURR_TEST_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/test")

# Tests for MarmotAutomaticDifferentiation
add_marmot_test("TestMarmotAutomaticDifferentiation_executable" "${CURR_TEST_SOURCE_DIR}/TestMarmotAutomaticDifferentiation.cpp")

# Tests for MarmotMath
add_marmot_test("TestMarmotMath_executable" "${CURR_TEST_SOURCE_DIR}/TestMarmotMath.cpp")

# Tests for MarmotNumericalDifferentiation
add_marmot_test("TestMarmotNumericalDifferentiation_executable" "${CURR_TEST_SOURCE_DIR}/TestMarmotNumericalDifferentiation.cpp")

# Tests for MarmotNumericalIntegration
add_marmot_test("TestMarmotNumericalIntegration_executable" "${CURR_TEST_SOURCE_DIR}/TestMarmotNumericalIntegration.cpp")

# Tests for MarmotTensor
add_marmot_test("TestMarmotTensor_executable" "${CURR_TEST_SOURCE_DIR}/TestMarmotTensor.cpp")

# Tests for NewtonConvergenceChecker
add_marmot_test("TestNewtonConvergenceChecker_executable" "${CURR_TEST_SOURCE_DIR}/TestNewtonConvergenceChecker.cpp")
