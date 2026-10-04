# add current directory to the source for the tests
SET(CURR_TEST_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/test")

# Tests for MarmotAutomaticDifferentiation
add_marmot_test("TestMarmotAutomaticDifferentiation" "${CURR_TEST_SOURCE_DIR}/TestMarmotAutomaticDifferentiation.cpp" REQUIRES MarmotMechanicsCore)

# Tests for MarmotAutomaticDifferentiation
add_marmot_test("TestMarmotAutomaticDifferentiationForFastor" "${CURR_TEST_SOURCE_DIR}/TestMarmotAutomaticDifferentiationForFastor.cpp" REQUIRES MarmotFiniteStrainMechanicsCore)

# Tests for MarmotMath
add_marmot_test("TestMarmotMath" "${CURR_TEST_SOURCE_DIR}/TestMarmotMath.cpp" REQUIRES MarmotMechanicsCore)

# Tests for MarmotNumericalDifferentiation
add_marmot_test("TestMarmotNumericalDifferentiation" "${CURR_TEST_SOURCE_DIR}/TestMarmotNumericalDifferentiation.cpp" REQUIRES MarmotMechanicsCore)

# Tests for MarmotNumericalDifferentiation
add_marmot_test("TestMarmotNumericalDifferentiationForFastor" "${CURR_TEST_SOURCE_DIR}/TestMarmotNumericalDifferentiationForFastor.cpp" REQUIRES MarmotMechanicsCore)

# Tests for MarmotNumericalIntegration
add_marmot_test("TestMarmotNumericalIntegration" "${CURR_TEST_SOURCE_DIR}/TestMarmotNumericalIntegration.cpp" REQUIRES MarmotMechanicsCore)

# Tests for MarmotTensor
add_marmot_test("TestMarmotTensor" "${CURR_TEST_SOURCE_DIR}/TestMarmotTensor.cpp" REQUIRES MarmotMechanicsCore)

# Tests for NewtonConvergenceChecker
add_marmot_test("TestNewtonConvergenceChecker" "${CURR_TEST_SOURCE_DIR}/TestNewtonConvergenceChecker.cpp" REQUIRES MarmotMechanicsCore)

# Tests for MarmotTensorExponential
add_marmot_test("TestMarmotTensorExponential" "${CURR_TEST_SOURCE_DIR}/TestMarmotTensorExponential.cpp" REQUIRES MarmotMechanicsCore)

# Tests for MarmotEigenSystems
add_marmot_test("TestMarmotEigenSystems" "${CURR_TEST_SOURCE_DIR}/TestMarmotEigenSystems.cpp" REQUIRES MarmotMechanicsCore)

# Tests for MarmotFastorTensorBasics
add_marmot_test("TestMarmotFastorTensorBasics" "${CURR_TEST_SOURCE_DIR}/TestMarmotFastorTensorBasics.cpp" REQUIRES MarmotMechanicsCore)
