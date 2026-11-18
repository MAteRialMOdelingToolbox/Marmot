# add current directory to the source for the tests
SET(CURR_TEST_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/test")

# Tests for the TransportFiniteElement (heat conduction elements DC*T)
add_marmot_test("TestTransportFiniteElement" "${CURR_TEST_SOURCE_DIR}/TestTransportFiniteElement.cpp"
                REQUIRES HeatConduction)

# Tests for the moisture diffusion elements DC*H
add_marmot_test("TestMoistureTransportElements" "${CURR_TEST_SOURCE_DIR}/TestMoistureTransportElements.cpp"
                REQUIRES MoistureDiffusion)
