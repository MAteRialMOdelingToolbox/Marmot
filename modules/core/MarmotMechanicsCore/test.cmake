# add current directory to the source for the tests
SET(CURR_TEST_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/test")

# Tests for HaighWestergaard
add_marmot_test("TestHaighWestergaard" "${CURR_TEST_SOURCE_DIR}/TestHaighWestergaard.cpp")

# Tests for MarmotElasticity
add_marmot_test("TestMarmotElasticity" "${CURR_TEST_SOURCE_DIR}/TestMarmotElasticity.cpp")

# Tests for MarmotKelvinChain
add_marmot_test("TestMarmotKelvinChain" "${CURR_TEST_SOURCE_DIR}/TestMarmotKelvinChain.cpp")

# Tests for MarmotKinematics
add_marmot_test("TestMarmotKinematics" "${CURR_TEST_SOURCE_DIR}/TestMarmotKinematics.cpp")

# Tests for MarmotLowerDimensionalStress
add_marmot_test("TestMarmotLowerDimensionalStress" "${CURR_TEST_SOURCE_DIR}/TestMarmotLowerDimensionalStress.cpp")

# Tests for MarmotLocalization
add_marmot_test("TestMarmotLocalization" "${CURR_TEST_SOURCE_DIR}/TestMarmotLocalization.cpp")

# Tests for MarmotPronySeries
add_marmot_test("TestMarmotPronySeries" "${CURR_TEST_SOURCE_DIR}/TestMarmotPronySeries.cpp")

# Tests for MarmotStateVarVectorManager
add_marmot_test("TestMarmotStateVarVectorManager" "${CURR_TEST_SOURCE_DIR}/TestMarmotStateVarVectorManager.cpp")

# Tests for MarmotViscoelasticity
add_marmot_test("TestMarmotViscoelasticity" "${CURR_TEST_SOURCE_DIR}/TestMarmotViscoelasticity.cpp")

# Tests for MarmotVoigt
add_marmot_test("TestMarmotVoigt" "${CURR_TEST_SOURCE_DIR}/TestMarmotVoigt.cpp")

# Tests for MenetreyWillam
add_marmot_test("TestMenetreyWillam" "${CURR_TEST_SOURCE_DIR}/TestMenetreyWillam.cpp")

# Tests for YieldSurfaceCombinatioNmanager
add_marmot_test("TestYieldSurfaceCombinationManager" "${CURR_TEST_SOURCE_DIR}/TestYieldSurfaceCombinationManager.cpp")

# Tests for NewmarkBetaIntegrator
add_marmot_test("TestNewmarkBetaIntegrator" "${CURR_TEST_SOURCE_DIR}/TestNewmarkBetaIntegrator.cpp")
