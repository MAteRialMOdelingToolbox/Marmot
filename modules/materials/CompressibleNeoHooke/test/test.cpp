#include "Marmot/CompressibleNeoHooke.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTesting.h"

using namespace Marmot::Testing;
using namespace Marmot::Materials;
using namespace Marmot::FastorStandardTensors;
using namespace Marmot::FastorIndices;

// 1. Test: Undeformed configuration
void testUndeformedResponse()
{
  // idx 0 - Bulk modulus K, idx 1 - Shear modulus G
  std::array< double, 2 > materialProperties_ = { 3500, 1500 };
  const double            nMaterialProperties = 2;
  const int               elLabel             = 1;

  // Create material instance
  CompressibleNeoHooke mat = CompressibleNeoHooke( &materialProperties_[0], nMaterialProperties, elLabel );

  // Create deformation, time increment, response and tangent objects required for stress computation
  CompressibleNeoHooke::Deformation< 3 > def;
  CompressibleNeoHooke::TimeIncrement    timeInc = { 0, 0.1 };

  CompressibleNeoHooke::ConstitutiveResponse< 3 > response;
  CompressibleNeoHooke::AlgorithmicModuli< 3 >    tangent;

  // Prescribe a deformation gradient tensor F for the considered load case
  def.F = Marmot::FastorStandardTensors::Spatial3D::I;

  // Compute stress response
  mat.computeStress( response, tangent, def, timeInc );

  // Target Kirchhoff stress values for the given deformation gradient and material parameters
  Tensor33d stressTarget( 0.0 );

  // Compare computed stress to target stress values
  throwExceptionOnFailure( checkIfEqual( response.tau, stressTarget, 1e-10 ),
                           "1. Test: Undeformed configuration failed for CompressibleNeoHooke material in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// Test 2: Finite strain simple shear load case
void testFiniteStrainShearResponse()
{

  std::array< double, 2 > materialProperties_ = { 3500, 1500 };
  const double            nMaterialProperties = 2;
  const int               elLabel             = 1;

  CompressibleNeoHooke mat = CompressibleNeoHooke( &materialProperties_[0], nMaterialProperties, elLabel );

  CompressibleNeoHooke::Deformation< 3 > def;
  CompressibleNeoHooke::TimeIncrement    timeInc = { 0, 0.1 };

  CompressibleNeoHooke::ConstitutiveResponse< 3 > response;
  CompressibleNeoHooke::AlgorithmicModuli< 3 >    tangent;

  def.F = Marmot::FastorStandardTensors::Spatial3D::I;
  def.F( 1, 0 ) += 0.2;

  mat.computeStress( response, tangent, def, timeInc );

  Tensor33d stressTarget( 0.0 );
  stressTarget( 0, 0 ) = -20;
  stressTarget( 0, 1 ) = 300;
  stressTarget( 1, 0 ) = 300;
  stressTarget( 1, 1 ) = 40;
  stressTarget( 2, 2 ) = -20;

  throwExceptionOnFailure( checkIfEqual( response.tau, stressTarget, 1e-10 ),
                           "Test 2: Finite strain simple shear load case failed for CompressibleNeoHooke material "
                           "in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// Test 3: Small strain simple shear load case
void testSmallStrainShearResponse()
{

  std::array< double, 2 > materialProperties_ = { 3500, 1500 };
  const double            nMaterialProperties = 2;
  const int               elLabel             = 1;

  CompressibleNeoHooke mat = CompressibleNeoHooke( &materialProperties_[0], nMaterialProperties, elLabel );

  CompressibleNeoHooke::Deformation< 3 > def;
  CompressibleNeoHooke::TimeIncrement    timeInc = { 0, 0.1 };

  CompressibleNeoHooke::ConstitutiveResponse< 3 > response;
  CompressibleNeoHooke::AlgorithmicModuli< 3 >    tangent;

  def.F = Marmot::FastorStandardTensors::Spatial3D::I;
  def.F( 1, 0 ) += 1e-06;

  mat.computeStress( response, tangent, def, timeInc );

  Tensor33d stressTarget( 0.0 );
  stressTarget( 0, 0 ) = -4.99994712299667e-10;
  stressTarget( 0, 1 ) = 0.0015;
  stressTarget( 1, 0 ) = 0.0015;
  stressTarget( 1, 1 ) = 9.99793777126387e-10;
  stressTarget( 2, 2 ) = -4.99994712299667e-10;

  throwExceptionOnFailure( checkIfEqual( response.tau, stressTarget, 1e-10 ),
                           "Test 3: Small strain simple shear failed for CompressibleNeoHooke material in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// Test 4: Hydrostatic load case
void testHydrostaticResponse()
{

  std::array< double, 2 > materialProperties_ = { 3500, 1500 };
  const double            nMaterialProperties = 2;
  const int               elLabel             = 1;

  CompressibleNeoHooke mat = CompressibleNeoHooke( &materialProperties_[0], nMaterialProperties, elLabel );

  CompressibleNeoHooke::Deformation< 3 > def;
  CompressibleNeoHooke::TimeIncrement    timeInc = { 0, 0.1 };

  CompressibleNeoHooke::ConstitutiveResponse< 3 > response;
  CompressibleNeoHooke::AlgorithmicModuli< 3 >    tangent;

  def.F = Marmot::FastorStandardTensors::Spatial3D::I;
  def.F( 0, 0 ) += 0.02;
  def.F( 1, 1 ) += 0.02;
  def.F( 2, 2 ) += 0.02;

  mat.computeStress( response, tangent, def, timeInc );

  Tensor33d stressTarget( 0.0 );
  stressTarget( 0, 0 ) = 208.417157443082;
  stressTarget( 1, 1 ) = 208.417157443082;
  stressTarget( 2, 2 ) = 208.417157443082;

  throwExceptionOnFailure( checkIfEqual( response.tau, stressTarget, 1e-10 ),
                           "Test 4: Hydrostatic load case computation failed for CompressibleNeoHooke material "
                           "in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// Test 5: General deformation load case
void testGeneralDeformationResponse()
{

  std::array< double, 2 > materialProperties_ = { 3500, 1500 };
  const double            nMaterialProperties = 2;
  const int               elLabel             = 1;

  CompressibleNeoHooke mat = CompressibleNeoHooke( &materialProperties_[0], nMaterialProperties, elLabel );

  CompressibleNeoHooke::Deformation< 3 > def;
  CompressibleNeoHooke::TimeIncrement    timeInc = { 0, 0.1 };

  CompressibleNeoHooke::ConstitutiveResponse< 3 > response;
  CompressibleNeoHooke::AlgorithmicModuli< 3 >    tangent;

  def.F( 0, 0 ) = 1.01;
  def.F( 0, 1 ) = 0.06;
  def.F( 0, 2 ) = -0.03;
  def.F( 1, 0 ) = 0.06;
  def.F( 1, 1 ) = 1.02;
  def.F( 1, 2 ) = 0.04;
  def.F( 2, 0 ) = -0.03;
  def.F( 2, 1 ) = 0.04;
  def.F( 2, 2 ) = 0.95;

  mat.computeStress( response, tangent, def, timeInc );

  Tensor33d stressTarget;
  stressTarget( 0, 0 ) = -47.0953127005558;
  stressTarget( 0, 1 ) = 184.282786939341;
  stressTarget( 0, 2 ) = -86.1819998621794;
  stressTarget( 1, 0 ) = 184.282786939341;
  stressTarget( 1, 1 ) = -15.0062701986801;
  stressTarget( 1, 2 ) = 117.659822506876;
  stressTarget( 2, 0 ) = -86.1819998621794;
  stressTarget( 2, 1 ) = 117.659822506876;
  stressTarget( 2, 2 ) = -229.85004999695;

  throwExceptionOnFailure( checkIfEqual( response.tau, stressTarget, 1e-10 ),
                           "Test 5a: General deformation response computation failed for CompressibleNeoHooke "
                           "material in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  for ( int i = 0; i < 3; i++ )
    for ( int j = 0; j < 3; j++ )

      throwExceptionOnFailure( checkIfEqual( response.tau( i, j ), response.tau( j, i ), 1e-10 ),
                               "Test 5b: Kirchhoff stress tensor symmetry check for the arbitrary deformation load "
                               "case failed for CompressibleNeoHooke "
                               "material in " +
                                 std::string( __PRETTY_FUNCTION__ ) );
}

// Test 6: Computation of the algorithmic tangent
void testAlgorithmicTangent()
{

  std::array< double, 2 > materialProperties_ = { 3500, 1500 };
  const double            nMaterialProperties = 2;
  const int               elLabel             = 1;

  CompressibleNeoHooke mat = CompressibleNeoHooke( &materialProperties_[0], nMaterialProperties, elLabel );

  CompressibleNeoHooke::Deformation< 3 > def;
  CompressibleNeoHooke::TimeIncrement    timeInc = { 0, 0.1 };

  CompressibleNeoHooke::ConstitutiveResponse< 3 > response;
  CompressibleNeoHooke::AlgorithmicModuli< 3 >    tangent;

  def.F = Marmot::FastorStandardTensors::Spatial3D::I;
  def.F( 0, 0 ) += 0.01;
  def.F( 1, 1 ) += 0.02;
  def.F( 2, 2 ) += 0.03;

  mat.computeStress( response, tangent, def, timeInc );

  Tensor3333d tangentTarget( 0.0 );

  tangentTarget( 0, 0, 0, 0 ) = 5450.8251046444;
  tangentTarget( 0, 0, 1, 1 ) = 2494.28143117259;
  tangentTarget( 0, 0, 2, 2 ) = 2450.93382241823;
  tangentTarget( 0, 1, 0, 1 ) = 1470.68247507597;
  tangentTarget( 0, 1, 1, 0 ) = 1456.26401943797;
  tangentTarget( 0, 2, 0, 2 ) = 1485.10093071397;
  tangentTarget( 0, 2, 2, 0 ) = 1456.26401943797;
  tangentTarget( 1, 0, 0, 1 ) = 1470.68247507597;
  tangentTarget( 1, 0, 1, 0 ) = 1456.26401943797;
  tangentTarget( 1, 1, 0, 0 ) = 2518.97728692678;
  tangentTarget( 1, 1, 1, 1 ) = 5416.51601207936;
  tangentTarget( 1, 1, 2, 2 ) = 2431.98918491329;
  tangentTarget( 1, 2, 1, 2 ) = 1485.10093071397;
  tangentTarget( 1, 2, 2, 1 ) = 1470.68247507597;
  tangentTarget( 2, 0, 0, 2 ) = 1485.10093071397;
  tangentTarget( 2, 0, 2, 0 ) = 1456.26401943797;
  tangentTarget( 2, 1, 1, 2 ) = 1485.10093071397;
  tangentTarget( 2, 1, 2, 1 ) = 1470.68247507597;
  tangentTarget( 2, 2, 0, 0 ) = 2499.46716543642;
  tangentTarget( 2, 2, 1, 1 ) = 2455.83221613793;
  tangentTarget( 2, 2, 2, 2 ) = 5383.05976216137;

  throwExceptionOnFailure( checkIfEqual( tangent.dTau_dF, tangentTarget, 1e-10 ),
                           "Test 6: Algorithmic tangent computation failed for CompressibleNeoHooke material "
                           "in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// Test 7.
void testRotation()
{
  // Test 8a: Pure rotation about z-axis
  {
    std::array< double, 2 > materialProperties_ = { 3500, 1500 };
    const double            nMaterialProperties = 2;
    const int               elLabel             = 1;

    CompressibleNeoHooke mat = CompressibleNeoHooke( &materialProperties_[0], nMaterialProperties, elLabel );

    CompressibleNeoHooke::Deformation< 3 > def;
    CompressibleNeoHooke::TimeIncrement    timeInc = { 0, 0.1 };

    CompressibleNeoHooke::ConstitutiveResponse< 3 > response;
    CompressibleNeoHooke::AlgorithmicModuli< 3 >    tangent;

    for ( int phi_deg = 0; phi_deg <= 180; phi_deg++ ) {
      double phi = Marmot::Math::degToRad( phi_deg );

      def.F( 0, 0 ) = cos( phi );
      def.F( 0, 1 ) = -sin( phi );
      def.F( 0, 2 ) = 0;
      def.F( 1, 0 ) = sin( phi );
      def.F( 1, 1 ) = cos( phi );
      def.F( 1, 2 ) = 0;
      def.F( 2, 0 ) = 0;
      def.F( 2, 1 ) = 0;
      def.F( 2, 2 ) = 1;

      mat.computeStress( response, tangent, def, timeInc );

      Tensor33d stressTarget( 0.0 );

      throwExceptionOnFailure( checkIfEqual( response.tau, stressTarget, 1e-10 ),
                               "Test 7a: Pure rotation around the z-axis computation failed for CompressibleNeoHooke "
                               "material in " +
                                 std::string( __PRETTY_FUNCTION__ ) );
    }
  }

  // Test 7b: Objectivity test for arbitrary deformation and rotations about the z-axis
  {
    std::array< double, 2 > materialProperties_ = { 3500, 1500 };
    const double            nMaterialProperties = 2;
    const int               elLabel             = 1;

    CompressibleNeoHooke mat = CompressibleNeoHooke( &materialProperties_[0], nMaterialProperties, elLabel );

    CompressibleNeoHooke::Deformation< 3 > def;
    CompressibleNeoHooke::TimeIncrement    timeInc = { 0, 0.1 };

    CompressibleNeoHooke::ConstitutiveResponse< 3 > response;
    CompressibleNeoHooke::AlgorithmicModuli< 3 >    tangent;

    def.F( 0, 0 ) = 1.01;
    def.F( 0, 1 ) = 0.06;
    def.F( 0, 2 ) = -0.03;
    def.F( 1, 0 ) = 0.06;
    def.F( 1, 1 ) = 1.02;
    def.F( 1, 2 ) = 0.04;
    def.F( 2, 0 ) = -0.03;
    def.F( 2, 1 ) = 0.04;
    def.F( 2, 2 ) = 0.95;

    mat.computeStress( response, tangent, def, timeInc );

    Tensor33d stressUnrotated = response.tau;
    Tensor33d F_unrotated     = def.F;

    for ( int phi_deg = 0; phi_deg <= 180; phi_deg++ ) {
      double phi = Marmot::Math::degToRad( phi_deg );

      Tensor33d Q( 0.0 );
      Q( 0, 0 ) = cos( phi );
      Q( 0, 1 ) = -sin( phi );
      Q( 1, 0 ) = sin( phi );
      Q( 1, 1 ) = cos( phi );
      Q( 2, 2 ) = 1;

      // Fr = Q * F -> Fr_ij = Q_ik F_kj

      Tensor33d F_rotated = einsum< ik, kj, to_ij >( Q, F_unrotated );

      def.F = F_rotated;

      mat.computeStress( response, tangent, def, timeInc );

      Tensor33d stressNew = response.tau;

      // Tau (Q*F) = Q * Tau(F) * Q^T -> Tau(Q*F)_ij = Q_iI Tau(F)_IJ Q_jJ
      Tensor33d stressRotated = einsum< iI, IJ, jJ, to_ij >( Q, stressUnrotated, Q );

      throwExceptionOnFailure( checkIfEqual( stressNew, stressRotated, 1e-10 ),
                               "Test 7b: Objectivity test for arbitrary deformation and rotation around z-axis failed "
                               "for "
                               "CompressibleNeoHooke "
                               "material in " +
                                 std::string( __PRETTY_FUNCTION__ ) );
    }
  }

  // Test 7c: Isotropy test for arbitrary deformation
  {
    std::array< double, 2 > materialProperties_ = { 3500, 1500 };
    const double            nMaterialProperties = 2;
    const int               elLabel             = 1;

    CompressibleNeoHooke mat = CompressibleNeoHooke( &materialProperties_[0], nMaterialProperties, elLabel );

    CompressibleNeoHooke::Deformation< 3 > def;
    CompressibleNeoHooke::TimeIncrement    timeInc = { 0, 0.1 };

    CompressibleNeoHooke::ConstitutiveResponse< 3 > response;
    CompressibleNeoHooke::AlgorithmicModuli< 3 >    tangent;

    def.F( 0, 0 ) = 1.01;
    def.F( 0, 1 ) = 0.06;
    def.F( 0, 2 ) = -0.03;
    def.F( 1, 0 ) = 0.06;
    def.F( 1, 1 ) = 1.02;
    def.F( 1, 2 ) = 0.04;
    def.F( 2, 0 ) = -0.03;
    def.F( 2, 1 ) = 0.04;
    def.F( 2, 2 ) = 0.95;

    mat.computeStress( response, tangent, def, timeInc );

    Tensor33d stressUnrotated = response.tau;
    Tensor33d F_unrotated     = def.F;

    for ( int phi_deg = 0; phi_deg <= 180; phi_deg++ ) {
      double phi = Marmot::Math::degToRad( phi_deg );

      Tensor33d Q( 0.0 );
      Q( 0, 0 ) = cos( phi );
      Q( 0, 1 ) = -sin( phi );
      Q( 1, 0 ) = sin( phi );
      Q( 1, 1 ) = cos( phi );
      Q( 2, 2 ) = 1;

      // Fr = F * Q -> Fr_ij = F_ik Q_kj

      Tensor33d F_rotated = einsum< ik, kj, to_ij >( F_unrotated, Q );

      def.F = F_rotated;

      mat.computeStress( response, tangent, def, timeInc );

      Tensor33d stressNew = response.tau;

      throwExceptionOnFailure( checkIfEqual( stressNew, stressUnrotated, 1e-10 ),
                               "Test 7c: Isotropy test for arbitrary deformation failed "
                               "for "
                               "CompressibleNeoHooke "
                               "material in " +
                                 std::string( __PRETTY_FUNCTION__ ) );
    }
  }
}

int main()
{

  auto tests = std::vector< std::function< void() > >{ testUndeformedResponse,
                                                       testFiniteStrainShearResponse,
                                                       testSmallStrainShearResponse,
                                                       testHydrostaticResponse,
                                                       testGeneralDeformationResponse,
                                                       testAlgorithmicTangent,
                                                       testRotation };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
