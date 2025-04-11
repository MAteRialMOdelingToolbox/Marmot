# add current directory to the source for the tests
SET(CURR_TEST_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/test")

# Tests for HaighWestergaard
add_marmot_test("TestHaighWestergaard" "${CURR_TEST_SOURCE_DIR}/TestHaighWestergaard.cpp")

# Tests for HughesWinget
add_marmot_test("TestHughesWinget" "${CURR_TEST_SOURCE_DIR}/TestHughesWinget.cpp")

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

# Tests for MarmotMaterialHypoElastic
add_marmot_test("TestMarmotMaterialHypoElastic" "${CURR_TEST_SOURCE_DIR}/TestMarmotMaterialHypoElastic.cpp")

# Tests for MarmotMaterialHypoElasticAD
add_marmot_test("TestMarmotMaterialHypoElasticAD" "${CURR_TEST_SOURCE_DIR}/TestMarmotMaterialHypoElasticAD.cpp")

# Tests for MarmotMaterialHypoElasticNonLocal
add_marmot_test("TestMarmotMaterialHypoElasticNonLocal" "${CURR_TEST_SOURCE_DIR}/TestMarmotMaterialHypoElasticNonLocal.cpp")

# Tests for MarmotMaterialHyperElastic
add_marmot_test("TestMarmotMaterialHyperElastic" "${CURR_TEST_SOURCE_DIR}/TestMarmotMaterialHyperElastic.cpp")

# Tests for MarmotMaterialMechanical
add_marmot_test("TestMarmotMaterialMechanical" "${CURR_TEST_SOURCE_DIR}/TestMarmotMaterialMechanical.cpp")

# Tests for MarmotPronySeries
add_marmot_test("TestMarmotPronySeries" "${CURR_TEST_SOURCE_DIR}/TestMarmotPronySeries.cpp")

# Tests for MarmotStateVarVectorManager
add_marmot_test("TestMarmotStateVarVectorManager" "${CURR_TEST_SOURCE_DIR}/TestMarmotStateVarVectorManager.cpp")

# Tests for MarmotUtility
add_marmot_test("TestMarmotUtility" "${CURR_TEST_SOURCE_DIR}/TestMarmotUtility.cpp")

# Tests for MarmotViscoelasticity
add_marmot_test("TestMarmotViscoelasticity" "${CURR_TEST_SOURCE_DIR}/TestMarmotViscoelasticity.cpp")

# Tests for MarmotVoigt
add_marmot_test("TestMarmotVoigt" "${CURR_TEST_SOURCE_DIR}/TestMarmotVoigt.cpp")

# Tests for MenetreyWillam
add_marmot_test("TestMenetreyWillam" "${CURR_TEST_SOURCE_DIR}/TestMenetreyWillam.cpp")

# Tests for YieldSurfaceCombinatioNmanager
add_marmot_test("TestYieldSurfaceCombinationManager" "${CURR_TEST_SOURCE_DIR}/TestYieldSurfaceCombinationManager.cpp")
