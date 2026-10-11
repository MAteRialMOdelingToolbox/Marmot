#include "Marmot/CompressibleNeoHooke.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include "Marmot/MarmotMaterialPointSolverFiniteStrain.h"
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
  const int               nMaterialProperties = 2;
  const int               elLabel             = 1;

  // Create material instance
  const CompressibleNeoHooke mat = CompressibleNeoHooke( &materialProperties_[0], nMaterialProperties, elLabel );

  // Create deformation, time increment, response and tangent objects required for stress computation
  CompressibleNeoHooke::Deformation< 3 > def     = { Tensor33d( 0.0 ) };
  CompressibleNeoHooke::TimeIncrement    timeInc = { 0, 0.1 };

  CompressibleNeoHooke::ConstitutiveResponse< 3 > response( Tensor33d( 0.0 ), 0.0, 0.0, nullptr );
  CompressibleNeoHooke::AlgorithmicModuli< 3 >    tangent = { Tensor3333d( 0.0 ) };

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

void testWithMPSolver()
{
  using namespace Marmot::Solvers;
  auto        materialProperties = std::vector< double >{ 3500, 1500 };
  auto        solveropts         = MarmotMaterialPointSolverFiniteStrain::SolverOptions();
  std::string matName            = "COMPRESSIBLENEOHOOKE";
  auto        solver             = MarmotMaterialPointSolverFiniteStrain( matName,
                                                       materialProperties.data(),
                                                       materialProperties.size(),
                                                       solveropts );

  // create a step with controlled shear strain increment
  MarmotMaterialPointSolverFiniteStrain::Step step;
  step.gradUIncrementTarget         = Tensor33d( 0.0 );
  step.gradUIncrementTarget( 0, 1 ) = 0.0001;
  step.stressIncrementTarget        = Tensor33d( 0.0 );

  // set only shear gradU component to be controlled
  // both gradU12 and gradU21 are controlled
  step.isGradUComponentControlled         = Tensor33t< bool >( false );
  step.isGradUComponentControlled( 0, 1 ) = true;
  step.isGradUComponentControlled( 1, 0 ) = true;

  // set only shear stress component to be controlled
  // both tau12 and tau21 are not controlled
  step.isStressComponentControlled         = Tensor33t< bool >( true );
  step.isStressComponentControlled( 0, 1 ) = false;
  step.isStressComponentControlled( 1, 0 ) = false;
  step.timeStart                           = 0.0;
  step.timeEnd                             = 1;
  step.dTStart                             = .1;
  step.dTMax                               = 1;
  step.dTMin                               = 0.1;

  // add step to solver
  solver.addStep( step );

  // solve the material point problem
  solver.solve();

  // get the final stress state
  auto history     = solver.getHistory();
  auto finalStress = history.back().stress;

  Tensor33d stressTarget( 0.0 );
  stressTarget( 0, 1 ) = .15;
  stressTarget( 1, 0 ) = .15;

  throwExceptionOnFailure( checkIfEqual( finalStress, stressTarget, 1e-8 ),
                           "I-5: Material Point Solver simple shear test failed for CompressibleNeoHooke material in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}
// ─────────────────────────────────────────────────────────────────────────────
// The following tests exercise MarmotMaterialFiniteStrain's default (base-class)
// implementations, none of which CompressibleNeoHooke overrides: computeStressExplicit,
// computeStress(..., eigenDeformation), computePlaneStrain(Explicit) (with and without
// eigenDeformation), computePlaneStress(Explicit), getStateView, getMaximumWaveSpeed,
// setCharacteristicElementLength and findEigenDeformationForEigenStress.
// ─────────────────────────────────────────────────────────────────────────────

namespace {
  std::array< double, 2 > materialProperties_ = { 3500, 1500 };
  CompressibleNeoHooke    makeMaterial()
  {
    return CompressibleNeoHooke( &materialProperties_[0], 2, 1 );
  }
} // namespace

void testComputeStressExplicitMatchesComputeStress()
{
  CompressibleNeoHooke matA = makeMaterial();
  CompressibleNeoHooke matB = makeMaterial();

  Tensor33d inputF = Marmot::FastorStandardTensors::Spatial3D::I;
  inputF( 0, 1 ) += 0.05;

  CompressibleNeoHooke::Deformation< 3 >          def     = { inputF };
  CompressibleNeoHooke::TimeIncrement             timeInc = { 0, 0.1 };
  CompressibleNeoHooke::ConstitutiveResponse< 3 > resA( Tensor33d( 0.0 ), 0.0, 0.0, nullptr );
  CompressibleNeoHooke::AlgorithmicModuli< 3 >    tanA = { Tensor3333d( 0.0 ) };
  matA.computeStress( resA, tanA, def, timeInc );

  CompressibleNeoHooke::ConstitutiveResponse< 3 > resB( Tensor33d( 0.0 ), 0.0, 0.0, nullptr );
  matB.computeStressExplicit( resB, def, timeInc );

  throwExceptionOnFailure( checkIfEqual( resA.tau, resB.tau, 1e-14 ) &&
                             checkIfEqual( resA.elasticEnergyDensity, resB.elasticEnergyDensity, 1e-14 ),
                           "computeStressExplicit() must match computeStress() ignoring the tangent in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

void testComputeStressWithEigenDeformation()
{
  CompressibleNeoHooke mat = makeMaterial();

  Tensor33d inputF = Marmot::FastorStandardTensors::Spatial3D::I;
  inputF( 0, 0 )   = 1.01;
  inputF( 1, 1 )   = 1.02;
  inputF( 2, 2 )   = 0.95;

  const double F0_XX = 1.03, F0_YY = 0.98, F0_ZZ = 1.01;

  // Reference: compute directly at the eigen-deformation-scaled F (only the diagonal is scaled).
  Tensor33d fRef = inputF;
  fRef( 0, 0 ) *= F0_XX;
  fRef( 1, 1 ) *= F0_YY;
  fRef( 2, 2 ) *= F0_ZZ;

  CompressibleNeoHooke::Deformation< 3 >          defRef  = { fRef };
  CompressibleNeoHooke::TimeIncrement             timeInc = { 0, 0.1 };
  CompressibleNeoHooke::ConstitutiveResponse< 3 > resRef( Tensor33d( 0.0 ), 0.0, 0.0, nullptr );
  CompressibleNeoHooke::AlgorithmicModuli< 3 >    tanRef = { Tensor3333d( 0.0 ) };
  mat.computeStress( resRef, tanRef, defRef, timeInc );

  CompressibleNeoHooke::Deformation< 3 >          def = { inputF };
  CompressibleNeoHooke::ConstitutiveResponse< 3 > res( Tensor33d( 0.0 ), 0.0, 0.0, nullptr );
  CompressibleNeoHooke::AlgorithmicModuli< 3 >    tan = { Tensor3333d( 0.0 ) };
  // CompressibleNeoHooke's own 4-arg computeStress() override hides the base class's 5-arg
  // (eigenDeformation) overload for a concretely-typed instance -- call through an explicit
  // base-class reference, as production code does (always via a MarmotMaterialFiniteStrain*).
  MarmotMaterialFiniteStrain& matBase = mat;
  matBase.computeStress( res, tan, def, timeInc, { F0_XX, F0_YY, F0_ZZ } );

  throwExceptionOnFailure( checkIfEqual( res.tau, resRef.tau, 1e-10 ),
                           "computeStress() with eigenDeformation must match computing directly at the "
                           "eigen-deformation-scaled F in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  // The tangent's diagonal-vs-diagonal blocks are additionally chain-ruled by F0 (dF'_kk/dF_kk = F0_kk);
  // the remaining entries are unaffected, since only the diagonal of F is rescaled.
  const double F0[3]     = { F0_XX, F0_YY, F0_ZZ };
  bool         tangentOk = true;
  for ( int i = 0; i < 3; i++ )
    for ( int j = 0; j < 3; j++ )
      for ( int k = 0; k < 3; k++ )
        tangentOk = tangentOk && checkIfEqual( tan.dTau_dF( i, j, k, k ), tanRef.dTau_dF( i, j, k, k ) * F0[k], 1e-8 );

  throwExceptionOnFailure( tangentOk,
                           "computeStress() with eigenDeformation did not chain-rule the tangent's diagonal "
                           "blocks by F0 in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

void testComputePlaneStrainMatchesComputeStress()
{
  CompressibleNeoHooke mat = makeMaterial();

  Tensor33d inputF = Marmot::FastorStandardTensors::Spatial3D::I;
  inputF( 0, 0 )   = 1.02;
  inputF( 1, 1 )   = 0.99;
  inputF( 2, 2 )   = 1.0; // plane strain: no out-of-plane stretch

  CompressibleNeoHooke::Deformation< 3 > def     = { inputF };
  CompressibleNeoHooke::TimeIncrement    timeInc = { 0, 0.1 };

  CompressibleNeoHooke::ConstitutiveResponse< 3 > resStress( Tensor33d( 0.0 ), 0.0, 0.0, nullptr );
  CompressibleNeoHooke::AlgorithmicModuli< 3 >    tanStress = { Tensor3333d( 0.0 ) };
  mat.computeStress( resStress, tanStress, def, timeInc );

  CompressibleNeoHooke::ConstitutiveResponse< 3 > resStrain( Tensor33d( 0.0 ), 0.0, 0.0, nullptr );
  CompressibleNeoHooke::AlgorithmicModuli< 3 >    tanStrain = { Tensor3333d( 0.0 ) };
  mat.computePlaneStrain( resStrain, tanStrain, def, timeInc );

  throwExceptionOnFailure( checkIfEqual( resStress.tau, resStrain.tau, 1e-14 ) &&
                             checkIfEqual( tanStress.dTau_dF, tanStrain.dTau_dF, 1e-14 ),
                           "computePlaneStrain() must match computeStress() for MarmotMaterialFiniteStrain's "
                           "default (pass-through) implementation in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  // ... and the eigen-deformation overload must match computeStress(..., eigenDeformation) likewise.
  CompressibleNeoHooke::ConstitutiveResponse< 3 > resStressEigen( Tensor33d( 0.0 ), 0.0, 0.0, nullptr );
  CompressibleNeoHooke::AlgorithmicModuli< 3 >    tanStressEigen  = { Tensor3333d( 0.0 ) };
  MarmotMaterialFiniteStrain&                     matBaseForEigen = mat;
  matBaseForEigen.computeStress( resStressEigen, tanStressEigen, def, timeInc, { 1.01, 0.99, 1.02 } );

  CompressibleNeoHooke::ConstitutiveResponse< 3 > resStrainEigen( Tensor33d( 0.0 ), 0.0, 0.0, nullptr );
  CompressibleNeoHooke::AlgorithmicModuli< 3 >    tanStrainEigen = { Tensor3333d( 0.0 ) };
  mat.computePlaneStrain( resStrainEigen, tanStrainEigen, def, timeInc, { 1.01, 0.99, 1.02 } );

  throwExceptionOnFailure( checkIfEqual( resStressEigen.tau, resStrainEigen.tau, 1e-14 ),
                           "computePlaneStrain() with eigenDeformation must match computeStress() with "
                           "eigenDeformation in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  // computePlaneStrainExplicit (both overloads) must match computePlaneStrain ignoring the tangent.
  CompressibleNeoHooke::ConstitutiveResponse< 3 > resExplicit( Tensor33d( 0.0 ), 0.0, 0.0, nullptr );
  mat.computePlaneStrainExplicit( resExplicit, def, timeInc );
  throwExceptionOnFailure( checkIfEqual( resExplicit.tau, resStrain.tau, 1e-14 ),
                           "computePlaneStrainExplicit() must match computePlaneStrain() in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  CompressibleNeoHooke::ConstitutiveResponse< 3 > resExplicitEigen( Tensor33d( 0.0 ), 0.0, 0.0, nullptr );
  mat.computePlaneStrainExplicit( resExplicitEigen, def, timeInc, { 1.01, 0.99, 1.02 } );
  throwExceptionOnFailure( checkIfEqual( resExplicitEigen.tau, resStrainEigen.tau, 1e-14 ),
                           "computePlaneStrainExplicit() with eigenDeformation must match computePlaneStrain() "
                           "with eigenDeformation in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

void testComputePlaneStressIsNotYetImplemented()
{
  CompressibleNeoHooke mat = makeMaterial();

  Fastor::Tensor< double, 2, 2 > inputF2D( 0.0 );
  inputF2D( 0, 0 )                                        = 1.0;
  inputF2D( 1, 1 )                                        = 1.0;
  CompressibleNeoHooke::Deformation< 2 >          def2D   = { inputF2D };
  CompressibleNeoHooke::TimeIncrement             timeInc = { 0, 0.1 };
  CompressibleNeoHooke::ConstitutiveResponse< 2 > res2D( Fastor::Tensor< double, 2, 2 >( 0.0 ), 0.0, 0.0, nullptr );
  CompressibleNeoHooke::AlgorithmicModuli< 2 >    tan2D = { Fastor::Tensor< double, 2, 2, 2, 2 >( 0.0 ) };

  bool threw = false;
  try {
    mat.computePlaneStress( res2D, tan2D, def2D, timeInc );
  }
  catch ( const std::invalid_argument& ) {
    threw = true;
  }
  throwExceptionOnFailure( threw,
                           "computePlaneStress() (base class default) must throw std::invalid_argument in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  // computePlaneStressExplicit() delegates to computePlaneStress() and must propagate the throw.
  bool threwExplicit = false;
  try {
    mat.computePlaneStressExplicit( res2D, def2D, timeInc );
  }
  catch ( const std::invalid_argument& ) {
    threwExplicit = true;
  }
  throwExceptionOnFailure( threwExplicit,
                           "computePlaneStressExplicit() must propagate computePlaneStress()'s throw in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

void testGetStateViewThrowsForMaterialWithNoStateVars()
{
  // CompressibleNeoHooke requires no state variables, so any name lookup on its (empty)
  // state layout is necessarily unknown -- this documents that behaviour rather than
  // asserting a specific state variable exists.
  CompressibleNeoHooke mat = makeMaterial();

  bool threw = false;
  try {
    mat.getStateView( "anything", nullptr );
  }
  catch ( const std::runtime_error& ) {
    threw = true;
  }
  throwExceptionOnFailure( threw,
                           "getStateView() must throw for a material with no registered state variables in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

void testGetMaximumWaveSpeedMatchesAnalyticalTangentAtIdentity()
{
  // getDensity() requires a 3rd material property (density), unlike the K,G-only materials used
  // elsewhere in this file.
  std::array< double, 3 > propsWithDensity = { 3500, 1500, 2400 };
  CompressibleNeoHooke    mat( propsWithDensity.data(), 3, 1 );

  Tensor33d inputF = Marmot::FastorStandardTensors::Spatial3D::I;

  CompressibleNeoHooke::Deformation< 3 >          def     = { inputF };
  CompressibleNeoHooke::TimeIncrement             timeInc = { 0, 0.1 };
  CompressibleNeoHooke::ConstitutiveResponse< 3 > res( Tensor33d( 0.0 ), 0.0, 0.0, nullptr );
  CompressibleNeoHooke::AlgorithmicModuli< 3 >    tan = { Tensor3333d( 0.0 ) };
  mat.computeStress( res, tan, def, timeInc );

  const double maxDiag = std::max(
    { tan.dTau_dF( 0, 0, 0, 0 ), tan.dTau_dF( 1, 1, 1, 1 ), tan.dTau_dF( 2, 2, 2, 2 ) } );
  const double density  = mat.getDensity( nullptr );
  const double expected = std::sqrt( maxDiag / density );

  const double waveSpeed = mat.getMaximumWaveSpeed( nullptr, inputF );
  throwExceptionOnFailure( checkIfEqual( waveSpeed, expected, 1e-5 ),
                           "getMaximumWaveSpeed() does not match sqrt(max(C_ii)/rho) from the analytical tangent "
                           "at F=I in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

void testSetCharacteristicElementLengthIsANoOp()
{
  CompressibleNeoHooke mat = makeMaterial();
  // The base class default does nothing; this call exists purely to exercise that no-op body.
  mat.setCharacteristicElementLength( 0.5 );
}

void testFindEigenDeformationForEigenStressReproducesTargetStress()
{
  CompressibleNeoHooke mat = makeMaterial();

  const std::tuple< double, double, double > initialGuess = { 1.0, 1.0, 1.0 };
  const std::tuple< double, double, double > targetStress = { 100.0, 150.0, 200.0 };

  const auto [F0_XX, F0_YY, F0_ZZ] = mat.findEigenDeformationForEigenStress( initialGuess, targetStress, nullptr );

  Tensor33d fFound( 0.0 );
  fFound( 0, 0 ) = F0_XX;
  fFound( 1, 1 ) = F0_YY;
  fFound( 2, 2 ) = F0_ZZ;

  CompressibleNeoHooke::Deformation< 3 >          def     = { fFound };
  CompressibleNeoHooke::TimeIncrement             timeInc = { 0, 0.1 };
  CompressibleNeoHooke::ConstitutiveResponse< 3 > res( Tensor33d( 0.0 ), 0.0, 0.0, nullptr );
  CompressibleNeoHooke::AlgorithmicModuli< 3 >    tan = { Tensor3333d( 0.0 ) };
  mat.computeStress( res, tan, def, timeInc );

  throwExceptionOnFailure( checkIfEqual( res.tau( 0, 0 ), std::get< 0 >( targetStress ), 1e-6 ) &&
                             checkIfEqual( res.tau( 1, 1 ), std::get< 1 >( targetStress ), 1e-6 ) &&
                             checkIfEqual( res.tau( 2, 2 ), std::get< 2 >( targetStress ), 1e-6 ),
                           "findEigenDeformationForEigenStress() did not converge to a deformation reproducing "
                           "the target diagonal stress in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

int main()
{

  auto tests = std::vector< std::function< void() > >{
    testUndeformedResponse,
    testDeformationResponse,
    testAlgorithmicTangent,
    testRotation,
    testWithMPSolver,
    testComputeStressExplicitMatchesComputeStress,
    testComputeStressWithEigenDeformation,
    testComputePlaneStrainMatchesComputeStress,
    testComputePlaneStressIsNotYetImplemented,
    testGetStateViewThrowsForMaterialWithNoStateVars,
    testGetMaximumWaveSpeedMatchesAnalyticalTangentAtIdentity,
    testSetCharacteristicElementLengthIsANoOp,
    testFindEigenDeformationForEigenStressReproducesTargetStress,
  };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
