# add current directory to the source for the tests
SET(CURR_TEST_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/test")

add_marmot_test("TestGradientEnhancedFiniteStrainParticle" "${CURR_TEST_SOURCE_DIR}/test.cpp")
# glibc fills freed and fresh heap memory with a pattern: a member that is read before it is set shows up as garbage
# instead of as zero (the particles have been bitten by this: second moments, VCI basis and volume, density)
set_tests_properties("TestGradientEnhancedFiniteStrainParticle" PROPERTIES ENVIRONMENT "MALLOC_PERTURB_=165")
