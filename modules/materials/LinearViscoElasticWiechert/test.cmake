SET(CURR_TEST_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/test")

add_marmot_test("TestLinearViscoElasticWiechert" "${CURR_TEST_SOURCE_DIR}/test.cpp")
