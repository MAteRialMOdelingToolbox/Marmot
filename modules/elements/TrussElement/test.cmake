# add current directory to the source for the tests
SET(CURR_TEST_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/test")

# Tests for the TrussElement
add_marmot_test("TestTrussElement" "${CURR_TEST_SOURCE_DIR}/test.cpp" REQUIRES LinearElastic VonMises CompressibleNeoHooke FiniteStrainJ2Plasticity)

# Tests of the generic element interface of all truss types
add_marmot_test("TestTrussElementInterface" "${CURR_TEST_SOURCE_DIR}/testInterface.cpp" REQUIRES LinearElastic CompressibleNeoHooke)
