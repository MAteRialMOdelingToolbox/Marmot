# add current directory to the source for the tests
SET(CURR_TEST_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/test")

# Tests for the EmbeddedBondElement
add_marmot_test("TestEmbeddedBondElement" "${CURR_TEST_SOURCE_DIR}/test.cpp" REQUIRES LinearElasticBondSlip ModelCode2010BondSlip)

# Tests of the generic element interface of the bond elements (all host shapes)
add_marmot_test("TestEmbeddedBondElementInterface" "${CURR_TEST_SOURCE_DIR}/testInterface.cpp" REQUIRES LinearElasticBondSlip)
