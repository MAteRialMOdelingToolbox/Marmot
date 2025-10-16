#include "Marmot/CompressibleNeoHooke.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTesting.h"

using namespace Marmot::Testing;
using namespace Marmot::Materials;
using namespace Marmot::FastorStandardTensors;
using namespace Marmot::FastorIndices;

void testSetup( const std::string& testName,
                const Tensor33d&   inputF,
                const Tensor33d&   targetStress,
                bool               checkTangent     = false,
                const Tensor3333d& targetTangent    = Tensor3333d( 0.0 ),
                bool               ObjectivityCheck = false,
                bool               IsotropyCheck    = false )
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
  def.F = inputF;

  // Compute stress response
  mat.computeStress( response, tangent, def, timeInc );

  if ( ObjectivityCheck == false && IsotropyCheck == false ) {

    // Compare computed stress to target stress values
    throwExceptionOnFailure( checkIfEqual( response.tau, targetStress, 1e-10 ),
                             testName + " - Kirchhoff stress tensor (tau) computation failed" +
                               " for CompressibleNeoHooke material in " + std::string( __PRETTY_FUNCTION__ ) );

    for ( int i = 0; i < 3; i++ )
      for ( int j = 0; j < 3; j++ )

        throwExceptionOnFailure( checkIfEqual( response.tau( i, j ), response.tau( j, i ), 1e-10 ),
                                 testName + " - Kirchhoff stress tensor symmetry check failed" +
                                   " for CompressibleNeoHooke material in " + std::string( __PRETTY_FUNCTION__ ) );

    if ( checkTangent ) {
      // Compare algorithmic tangent to target tangent values
      throwExceptionOnFailure( checkIfEqual( tangent.dTau_dF, targetTangent, 1e-10 ),
                               testName + " - Algorithmic tangent tensor computation failed" +
                                 " for CompressibleNeoHooke material in " + std::string( __PRETTY_FUNCTION__ ) );
    }
  }

  // Check objectivity and isotropy if requested
  if ( ObjectivityCheck ) {
    // Use already computed stress response and current F
    Tensor33d stressUnrotated = response.tau;
    Tensor33d F_unrotated     = def.F;

    for ( int phi_deg = 0; phi_deg <= 180; phi_deg += 30 ) {

      double phi = Marmot::Math::degToRad( phi_deg );

      Tensor33d Q( 0.0 );
      Q( 0, 0 ) = cos( phi );
      Q( 0, 1 ) = -sin( phi );
      Q( 1, 0 ) = sin( phi );
      Q( 1, 1 ) = cos( phi );
      Q( 2, 2 ) = 1;

      // Fr = Q * F -> Fr_ij = Q_ik F_kj
      Tensor33d F_rotated = einsum< ik, kj, to_ij >( Q, F_unrotated );
      def.F               = F_rotated;

      mat.computeStress( response, tangent, def, timeInc );

      Tensor33d stressNew = response.tau;

      // Tau (Q*F) = Q * Tau(F) * Q^T -> Tau(Q*F)_ij = Q_iI Tau(F)_IJ Q_jJ
      Tensor33d stressRotated = einsum< iI, IJ, jJ, to_ij >( Q, stressUnrotated, Q );

      throwExceptionOnFailure( checkIfEqual( stressNew, stressRotated, 1e-10 ),
                               testName + " - Objectivity test failed (phi_deg=" + std::to_string( phi_deg ) +
                                 ") for CompressibleNeoHooke material in " + std::string( __PRETTY_FUNCTION__ ) );
    }
  }

  if ( IsotropyCheck ) {
    // Use already computed stress response and current deformed state
    Tensor33d stressUnrotated = response.tau;
    Tensor33d F_unrotated     = def.F;

    for ( int phi_deg = 0; phi_deg <= 180; phi_deg += 30 ) {
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
                               testName + " - Isotropy test failed (phi_deg=" + std::to_string( phi_deg ) +
                                 ") for CompressibleNeoHooke material in " + std::string( __PRETTY_FUNCTION__ ) );
    }
  }
}

// Test I-1: Undeformed configuration
void testUndeformedResponse()
{
  Tensor33d inputF = Marmot::FastorStandardTensors::Spatial3D::I;
  Tensor33d stressTarget( 0.0 );
  testSetup( "I-1: Undeformed configuration", inputF, stressTarget );
}

void testDeformationResponse()
{
  // Test I-2a: Finite strain simple shear load case
  {
    Tensor33d inputF = Marmot::FastorStandardTensors::Spatial3D::I;
    inputF( 1, 0 ) += 0.2;

    Tensor33d stressTarget( 0.0 );
    stressTarget( 0, 0 ) = -20;
    stressTarget( 0, 1 ) = 300;
    stressTarget( 1, 0 ) = 300;
    stressTarget( 1, 1 ) = 40;
    stressTarget( 2, 2 ) = -20;

    testSetup( "I-2a: Finite strain simple shear", inputF, stressTarget );
  }

  // Test I-2b: Small strain simple shear load case
  {
    Tensor33d inputF = Marmot::FastorStandardTensors::Spatial3D::I;
    inputF( 1, 0 ) += 1e-06;

    Tensor33d stressTarget( 0.0 );
    stressTarget( 0, 0 ) = -4.99994712299667e-10;
    stressTarget( 0, 1 ) = 0.0015;
    stressTarget( 1, 0 ) = 0.0015;
    stressTarget( 1, 1 ) = 9.99793777126387e-10;
    stressTarget( 2, 2 ) = -4.99994712299667e-10;

    testSetup( "I-2b: Small strain simple shear", inputF, stressTarget );
  }

  // Test I-2c: Hydrostatic load case
  {
    Tensor33d inputF = Marmot::FastorStandardTensors::Spatial3D::I;
    inputF( 0, 0 ) += 0.02;
    inputF( 1, 1 ) += 0.02;
    inputF( 2, 2 ) += 0.02;

    Tensor33d stressTarget( 0.0 );
    stressTarget( 0, 0 ) = 208.417157443082;
    stressTarget( 1, 1 ) = 208.417157443082;
    stressTarget( 2, 2 ) = 208.417157443082;
    testSetup( "I-2c: Hydrostatic", inputF, stressTarget );
  }

  // Test I-2d: Arbitrary deformation load case
  {
    Tensor33d inputF( 0.0 );
    inputF( 0, 0 ) = 1.01;
    inputF( 0, 1 ) = 0.06;
    inputF( 0, 2 ) = -0.03;
    inputF( 1, 0 ) = 0.06;
    inputF( 1, 1 ) = 1.02;
    inputF( 1, 2 ) = 0.04;
    inputF( 2, 0 ) = -0.03;
    inputF( 2, 1 ) = 0.04;
    inputF( 2, 2 ) = 0.95;

    Tensor33d stressTarget( 0.0 );
    stressTarget( 0, 0 ) = -47.0953127005558;
    stressTarget( 0, 1 ) = 184.282786939341;
    stressTarget( 0, 2 ) = -86.1819998621794;
    stressTarget( 1, 0 ) = 184.282786939341;
    stressTarget( 1, 1 ) = -15.0062701986801;
    stressTarget( 1, 2 ) = 117.659822506876;
    stressTarget( 2, 0 ) = -86.1819998621794;
    stressTarget( 2, 1 ) = 117.659822506876;
    stressTarget( 2, 2 ) = -229.85004999695;
    testSetup( "I-2d: Arbitrary deformation", inputF, stressTarget );
  }
}
// Test I-3: Computation of the algorithmic tangent
void testAlgorithmicTangent()
{
  Tensor33d inputF = Marmot::FastorStandardTensors::Spatial3D::I;
  inputF( 0, 0 ) += 0.01;
  inputF( 1, 1 ) += 0.02;
  inputF( 2, 2 ) += 0.03;

  Tensor33d stressTarget( 0.0 );
  stressTarget( 0, 0 ) = 178.712770583994;
  stressTarget( 1, 1 ) = 207.982235529133;
  stressTarget( 2, 2 ) = 237.540069587033;

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

  bool checkTangent = true;

  testSetup( "I-3: Algorithmic tangent", inputF, stressTarget, checkTangent, tangentTarget );
}

// Test I-4: Rotation tests
void testRotation()
{
  // Test I-4a: Pure rotation about z-axis
  {
    for ( int phi_deg = 0; phi_deg <= 180; phi_deg++ ) {
      double    phi = Marmot::Math::degToRad( phi_deg );
      Tensor33d inputF( 0.0 );
      inputF( 0, 0 ) = cos( phi );
      inputF( 0, 1 ) = -sin( phi );
      inputF( 0, 2 ) = 0;
      inputF( 1, 0 ) = sin( phi );
      inputF( 1, 1 ) = cos( phi );
      inputF( 1, 2 ) = 0;
      inputF( 2, 0 ) = 0;
      inputF( 2, 1 ) = 0;
      inputF( 2, 2 ) = 1;

      Tensor33d stressTarget( 0.0 );
      testSetup( "I-4a: Pure rotation (phi_deg=" + std::to_string( phi_deg ) + ")", inputF, stressTarget );
    }
  }

  // Test I-4b & c: Objectivity and Isotropy tests for arbitrary deformation and rotations about the z-axis
  {
    Tensor33d inputF( 0.0 );
    inputF( 0, 0 ) = 1.01;
    inputF( 0, 1 ) = 0.06;
    inputF( 0, 2 ) = -0.03;
    inputF( 1, 0 ) = 0.06;
    inputF( 1, 1 ) = 1.02;
    inputF( 1, 2 ) = 0.04;
    inputF( 2, 0 ) = -0.03;
    inputF( 2, 1 ) = 0.04;
    inputF( 2, 2 ) = 0.95;

    testSetup( "I-4b: Objectivity test", inputF, Tensor33d( 0.0 ), false, Tensor3333d( 0.0 ), true, false );
    testSetup( "I-4c: Isotropy test", inputF, Tensor33d( 0.0 ), false, Tensor3333d( 0.0 ), false, true );
  }
}

int main()
{

  auto tests = std::vector< std::function< void() > >{ testUndeformedResponse,
                                                       testDeformationResponse,
                                                       testAlgorithmicTangent,
                                                       testRotation };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
