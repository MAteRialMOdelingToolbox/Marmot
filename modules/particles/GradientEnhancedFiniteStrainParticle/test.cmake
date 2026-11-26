# add current directory to the source for the tests
SET(CURR_TEST_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/test")

add_marmot_test("TestGradientEnhancedFiniteStrainParticle" "${CURR_TEST_SOURCE_DIR}/test.cpp" REQUIRES GradientEnhancedCompressibleNeoHookeDamage GM7 MicroplaneM7)
# glibc fills freed and fresh heap memory with a pattern: a member that is read before it is set shows up as garbage
# instead of as zero (the particles have been bitten by this: second moments, VCI basis and volume, density)
if(TEST "TestGradientEnhancedFiniteStrainParticle") # not defined if a required module is not built
  set_tests_properties("TestGradientEnhancedFiniteStrainParticle" PROPERTIES ENVIRONMENT "MALLOC_PERTURB_=165")
endif()

# the axisymmetric particles, of both formulations, against patch tests, Lame's thick-walled cylinder and CX8RUL
add_marmot_test("TestAxisymmetricParticles" "${CURR_TEST_SOURCE_DIR}/testAxisymmetric.cpp" REQUIRES DisplacementParticle CompressibleNeoHooke GradientEnhancedCompressibleNeoHookeDamage DisplacementFiniteStrainULElement)
if(TEST "TestAxisymmetricParticles")
  set_tests_properties("TestAxisymmetricParticles" PROPERTIES ENVIRONMENT "MALLOC_PERTURB_=165")
endif()
