# add current directory to the source for the tests
SET(CURR_TEST_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/test")

# Tests for MarmotPhaseFieldEnergyDegradation
add_marmot_test("TestMarmotPhaseFieldEnergyDegradation" "${CURR_TEST_SOURCE_DIR}/TestMarmotPhaseFieldEnergyDegradation.cpp")

# Tests for MarmotDecreasingInteractions
add_marmot_test("TestMarmotDecreasingInteractions" "${CURR_TEST_SOURCE_DIR}/TestMarmotDecreasingInteractions.cpp")
