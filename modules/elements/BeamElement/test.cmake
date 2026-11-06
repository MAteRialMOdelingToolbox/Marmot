# add current directory to the source for the tests
SET(CURR_TEST_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/test")

# Tests for the BeamElement
add_marmot_test("TestBeamElement" "${CURR_TEST_SOURCE_DIR}/test.cpp" REQUIRES LinearElastic VonMises)

# Tests for the co-rotational beam B2D2CR
add_marmot_test("TestCorotationalBeamElement" "${CURR_TEST_SOURCE_DIR}/testCorotational.cpp" REQUIRES LinearElastic VonMises)

# Tests for the co-rotational beam B3D2CR
add_marmot_test("TestCorotationalBeamElement3D" "${CURR_TEST_SOURCE_DIR}/testCorotational3D.cpp" REQUIRES LinearElastic VonMises)

# Tests of the generic element interface and of invalid input for all beam types
add_marmot_test("TestBeamElementInterface" "${CURR_TEST_SOURCE_DIR}/testInterface.cpp" REQUIRES LinearElastic)
