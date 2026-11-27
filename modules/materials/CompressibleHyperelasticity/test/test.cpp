#include "Marmot/CompressibleHyperelasticity.h"
#include "Marmot/MarmotDeformationMeasures.h"
#include "Marmot/MarmotEnergyDensityFunctions.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include "Marmot/MarmotMaterialPointSolverFiniteStrain.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTesting.h"
#include <array>
#include <limits>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

using namespace Marmot::Testing;
using namespace Marmot::Materials;
using namespace Marmot::Math::Tensors::FastorTensors;
using namespace Marmot::Math::Tensors::FastorIndices;

void testSetup( const std::string& testName,
                const Tensor33d&   inputF,
                const Tensor33d&   targetStress,
                bool               checkTangent     = false,
                const Tensor3333d& targetTangent    = Tensor3333d( 0.0 ),
                bool               ObjectivityCheck = false,
                bool               IsotropyCheck    = false )
{

  // isochoric NeoHooke (0) with mu = 1500, volumetric DollSchweizerhof (6) with kappa = 3500, alpha = 1, beta = 3:
  // the compressible Neo-Hooke model of Pence & Gou (variant B) with K = 3500, G = 1500
  std::array< double, 6 > materialProperties_ = { 0, 6, 1500, 3500, 1, 3 };
  const int               nMaterialProperties = 6;
  const int               elLabel             = 1;

  // Create material instance
  const CompressibleHyperelasticity mat = CompressibleHyperelasticity( &materialProperties_[0],
                                                                       nMaterialProperties,
                                                                       elLabel );

  // Create deformation, time increment, response and tangent objects required for stress computation
  CompressibleHyperelasticity::Deformation< 3 > def     = { Tensor33d( 0.0 ) };
  CompressibleHyperelasticity::TimeIncrement    timeInc = { 0, 0.1 };

  CompressibleHyperelasticity::ConstitutiveResponse< 3 > response( Tensor33d( 0.0 ), 0.0, 0.0, nullptr );
  CompressibleHyperelasticity::AlgorithmicModuli< 3 >    tangent = { Tensor3333d( 0.0 ) };

  // Prescribe a deformation gradient tensor F for the considered load case
  def.F = inputF;

  // Compute stress response
  mat.computeStress( response, tangent, def, timeInc );

  if ( ObjectivityCheck == false && IsotropyCheck == false ) {

    // Compare computed stress to target stress values
    throwExceptionOnFailure( checkIfEqual( response.tau, targetStress, 1e-10 ),
                             testName + " - Kirchhoff stress tensor (tau) computation failed" +
                               " for CompressibleHyperelasticity material in " + std::string( __PRETTY_FUNCTION__ ) );

    for ( int i = 0; i < 3; i++ )
      for ( int j = 0; j < 3; j++ )

        throwExceptionOnFailure( checkIfEqual( response.tau( i, j ), response.tau( j, i ), 1e-10 ),
                                 testName + " - Kirchhoff stress tensor symmetry check failed" +
                                   " for CompressibleHyperelasticity material in " +
                                   std::string( __PRETTY_FUNCTION__ ) );

    if ( checkTangent ) {
      // Compare algorithmic tangent to target tangent values
      throwExceptionOnFailure( checkIfEqual( tangent.dTau_dF, targetTangent, 1e-10 ),
                               testName + " - Algorithmic tangent tensor computation failed" +
                                 " for CompressibleHyperelasticity material in " + std::string( __PRETTY_FUNCTION__ ) );
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
                                 ") for CompressibleHyperelasticity material in " +
                                 std::string( __PRETTY_FUNCTION__ ) );
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
                                 ") for CompressibleHyperelasticity material in " +
                                 std::string( __PRETTY_FUNCTION__ ) );
    }
  }
}

// Test I-1: Undeformed configuration
void testUndeformedResponse()
{
  Tensor33d inputF = Marmot::Math::Tensors::FastorTensors::Spatial3D::I;
  Tensor33d stressTarget( 0.0 );
  testSetup( "I-1: Undeformed configuration", inputF, stressTarget );
}

void testDeformationResponse()
{
  // Test I-2a: Finite strain simple shear load case
  {
    Tensor33d inputF = Marmot::Math::Tensors::FastorTensors::Spatial3D::I;
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
    Tensor33d inputF = Marmot::Math::Tensors::FastorTensors::Spatial3D::I;
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
    Tensor33d inputF = Marmot::Math::Tensors::FastorTensors::Spatial3D::I;
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
  Tensor33d inputF = Marmot::Math::Tensors::FastorTensors::Spatial3D::I;
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
  using namespace Marmot::Materials;
  auto        materialProperties = std::vector< double >{ 0, 6, 1500, 3500, 1, 3 };
  auto        solveropts         = Marmot::Materials::MaterialPointSolverFiniteStrain::SolverOptions();
  std::string matName            = "COMPRESSIBLEHYPERELASTICITY";
  auto        solver             = Marmot::Materials::MaterialPointSolverFiniteStrain( matName,
                                                                    materialProperties.data(),
                                                                    static_cast< int >( materialProperties.size() ),
                                                                    solveropts );

  // create a step with controlled shear strain increment
  Marmot::Materials::MaterialPointSolverFiniteStrain::Step step;
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
                           "I-5: Material Point Solver simple shear test failed for CompressibleHyperelasticity "
                           "material in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}
// ─────────────────────────────────────────────────────────────────────────────
// The following tests exercise Marmot::Materials::MaterialFiniteStrain's default (base-class)
// implementations, none of which CompressibleHyperelasticity overrides: computeStressExplicit,
// computeStress(..., eigenDeformation), computePlaneStrain(Explicit) (with and without
// eigenDeformation), computePlaneStress(Explicit), getStateView, getMaximumWaveSpeed,
// setCharacteristicElementLength and findEigenDeformationForEigenStress.
// ─────────────────────────────────────────────────────────────────────────────

namespace {
  std::array< double, 6 >     materialProperties_ = { 0, 6, 1500, 3500, 1, 3 };
  CompressibleHyperelasticity makeMaterial()
  {
    return CompressibleHyperelasticity( &materialProperties_[0], 6, 1 );
  }
} // namespace

void testComputeStressExplicitMatchesComputeStress()
{
  CompressibleHyperelasticity matA = makeMaterial();
  CompressibleHyperelasticity matB = makeMaterial();

  Tensor33d inputF = Marmot::Math::Tensors::FastorTensors::Spatial3D::I;
  inputF( 0, 1 ) += 0.05;

  CompressibleHyperelasticity::Deformation< 3 >          def     = { inputF };
  CompressibleHyperelasticity::TimeIncrement             timeInc = { 0, 0.1 };
  CompressibleHyperelasticity::ConstitutiveResponse< 3 > resA( Tensor33d( 0.0 ), 0.0, 0.0, nullptr );
  CompressibleHyperelasticity::AlgorithmicModuli< 3 >    tanA = { Tensor3333d( 0.0 ) };
  matA.computeStress( resA, tanA, def, timeInc );

  CompressibleHyperelasticity::ConstitutiveResponse< 3 > resB( Tensor33d( 0.0 ), 0.0, 0.0, nullptr );
  matB.computeStressExplicit( resB, def, timeInc );

  throwExceptionOnFailure( checkIfEqual( resA.tau, resB.tau, 1e-14 ) &&
                             checkIfEqual( resA.elasticEnergyDensity, resB.elasticEnergyDensity, 1e-14 ),
                           "computeStressExplicit() must match computeStress() ignoring the tangent in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

void testComputeStressWithEigenDeformation()
{
  CompressibleHyperelasticity mat = makeMaterial();

  Tensor33d inputF = Marmot::Math::Tensors::FastorTensors::Spatial3D::I;
  inputF( 0, 0 )   = 1.01;
  inputF( 1, 1 )   = 1.02;
  inputF( 2, 2 )   = 0.95;

  const double F0_XX = 1.03, F0_YY = 0.98, F0_ZZ = 1.01;

  // Reference: compute directly at the eigen-deformation-scaled F (only the diagonal is scaled).
  Tensor33d fRef = inputF;
  fRef( 0, 0 ) *= F0_XX;
  fRef( 1, 1 ) *= F0_YY;
  fRef( 2, 2 ) *= F0_ZZ;

  CompressibleHyperelasticity::Deformation< 3 >          defRef  = { fRef };
  CompressibleHyperelasticity::TimeIncrement             timeInc = { 0, 0.1 };
  CompressibleHyperelasticity::ConstitutiveResponse< 3 > resRef( Tensor33d( 0.0 ), 0.0, 0.0, nullptr );
  CompressibleHyperelasticity::AlgorithmicModuli< 3 >    tanRef = { Tensor3333d( 0.0 ) };
  mat.computeStress( resRef, tanRef, defRef, timeInc );

  CompressibleHyperelasticity::Deformation< 3 >          def = { inputF };
  CompressibleHyperelasticity::ConstitutiveResponse< 3 > res( Tensor33d( 0.0 ), 0.0, 0.0, nullptr );
  CompressibleHyperelasticity::AlgorithmicModuli< 3 >    tan = { Tensor3333d( 0.0 ) };
  // CompressibleHyperelasticity's own 4-arg computeStress() override hides the base class's 5-arg
  // (eigenDeformation) overload for a concretely-typed instance -- call through an explicit
  // base-class reference, as production code does (always via a Marmot::Materials::MaterialFiniteStrain*).
  Marmot::Materials::MaterialFiniteStrain& matBase = mat;
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
  CompressibleHyperelasticity mat = makeMaterial();

  Tensor33d inputF = Marmot::Math::Tensors::FastorTensors::Spatial3D::I;
  inputF( 0, 0 )   = 1.02;
  inputF( 1, 1 )   = 0.99;
  inputF( 2, 2 )   = 1.0; // plane strain: no out-of-plane stretch

  CompressibleHyperelasticity::Deformation< 3 > def     = { inputF };
  CompressibleHyperelasticity::TimeIncrement    timeInc = { 0, 0.1 };

  CompressibleHyperelasticity::ConstitutiveResponse< 3 > resStress( Tensor33d( 0.0 ), 0.0, 0.0, nullptr );
  CompressibleHyperelasticity::AlgorithmicModuli< 3 >    tanStress = { Tensor3333d( 0.0 ) };
  mat.computeStress( resStress, tanStress, def, timeInc );

  CompressibleHyperelasticity::ConstitutiveResponse< 3 > resStrain( Tensor33d( 0.0 ), 0.0, 0.0, nullptr );
  CompressibleHyperelasticity::AlgorithmicModuli< 3 >    tanStrain = { Tensor3333d( 0.0 ) };
  mat.computePlaneStrain( resStrain, tanStrain, def, timeInc );

  throwExceptionOnFailure( checkIfEqual( resStress.tau, resStrain.tau, 1e-14 ) &&
                             checkIfEqual( tanStress.dTau_dF, tanStrain.dTau_dF, 1e-14 ),
                           "computePlaneStrain() must match computeStress() for "
                           "Marmot::Materials::MaterialFiniteStrain's "
                           "default (pass-through) implementation in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  // ... and the eigen-deformation overload must match computeStress(..., eigenDeformation) likewise.
  CompressibleHyperelasticity::ConstitutiveResponse< 3 > resStressEigen( Tensor33d( 0.0 ), 0.0, 0.0, nullptr );
  CompressibleHyperelasticity::AlgorithmicModuli< 3 >    tanStressEigen  = { Tensor3333d( 0.0 ) };
  Marmot::Materials::MaterialFiniteStrain&               matBaseForEigen = mat;
  matBaseForEigen.computeStress( resStressEigen, tanStressEigen, def, timeInc, { 1.01, 0.99, 1.02 } );

  CompressibleHyperelasticity::ConstitutiveResponse< 3 > resStrainEigen( Tensor33d( 0.0 ), 0.0, 0.0, nullptr );
  CompressibleHyperelasticity::AlgorithmicModuli< 3 >    tanStrainEigen = { Tensor3333d( 0.0 ) };
  mat.computePlaneStrain( resStrainEigen, tanStrainEigen, def, timeInc, { 1.01, 0.99, 1.02 } );

  throwExceptionOnFailure( checkIfEqual( resStressEigen.tau, resStrainEigen.tau, 1e-14 ),
                           "computePlaneStrain() with eigenDeformation must match computeStress() with "
                           "eigenDeformation in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  // computePlaneStrainExplicit (both overloads) must match computePlaneStrain ignoring the tangent.
  CompressibleHyperelasticity::ConstitutiveResponse< 3 > resExplicit( Tensor33d( 0.0 ), 0.0, 0.0, nullptr );
  mat.computePlaneStrainExplicit( resExplicit, def, timeInc );
  throwExceptionOnFailure( checkIfEqual( resExplicit.tau, resStrain.tau, 1e-14 ),
                           "computePlaneStrainExplicit() must match computePlaneStrain() in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  CompressibleHyperelasticity::ConstitutiveResponse< 3 > resExplicitEigen( Tensor33d( 0.0 ), 0.0, 0.0, nullptr );
  mat.computePlaneStrainExplicit( resExplicitEigen, def, timeInc, { 1.01, 0.99, 1.02 } );
  throwExceptionOnFailure( checkIfEqual( resExplicitEigen.tau, resStrainEigen.tau, 1e-14 ),
                           "computePlaneStrainExplicit() with eigenDeformation must match computePlaneStrain() "
                           "with eigenDeformation in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

void testComputePlaneStressIsNotYetImplemented()
{
  CompressibleHyperelasticity mat = makeMaterial();

  Fastor::Tensor< double, 2, 2 > inputF2D( 0.0 );
  inputF2D( 0, 0 )                                               = 1.0;
  inputF2D( 1, 1 )                                               = 1.0;
  CompressibleHyperelasticity::Deformation< 2 >          def2D   = { inputF2D };
  CompressibleHyperelasticity::TimeIncrement             timeInc = { 0, 0.1 };
  CompressibleHyperelasticity::ConstitutiveResponse< 2 > res2D( Fastor::Tensor< double, 2, 2 >( 0.0 ),
                                                                0.0,
                                                                0.0,
                                                                nullptr );
  CompressibleHyperelasticity::AlgorithmicModuli< 2 >    tan2D = { Fastor::Tensor< double, 2, 2, 2, 2 >( 0.0 ) };

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
  // CompressibleHyperelasticity requires no state variables, so any name lookup on its (empty)
  // state layout is necessarily unknown -- this documents that behaviour rather than
  // asserting a specific state variable exists.
  CompressibleHyperelasticity mat = makeMaterial();

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
  // getDensity() requires the density after the parameters of the potential, unlike the materials used
  // elsewhere in this file.
  std::array< double, 7 >     propsWithDensity = { 0, 6, 1500, 3500, 1, 3, 2400 };
  CompressibleHyperelasticity mat( propsWithDensity.data(), 7, 1 );

  Tensor33d inputF = Marmot::Math::Tensors::FastorTensors::Spatial3D::I;

  CompressibleHyperelasticity::Deformation< 3 >          def     = { inputF };
  CompressibleHyperelasticity::TimeIncrement             timeInc = { 0, 0.1 };
  CompressibleHyperelasticity::ConstitutiveResponse< 3 > res( Tensor33d( 0.0 ), 0.0, 0.0, nullptr );
  CompressibleHyperelasticity::AlgorithmicModuli< 3 >    tan = { Tensor3333d( 0.0 ) };
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
  CompressibleHyperelasticity mat = makeMaterial();
  // The base class default does nothing; this call exists purely to exercise that no-op body.
  mat.setCharacteristicElementLength( 0.5 );
}

void testFindEigenDeformationForEigenStressReproducesTargetStress()
{
  CompressibleHyperelasticity mat = makeMaterial();

  const std::tuple< double, double, double > initialGuess = { 1.0, 1.0, 1.0 };
  const std::tuple< double, double, double > targetStress = { 100.0, 150.0, 200.0 };

  const auto [F0_XX, F0_YY, F0_ZZ] = mat.findEigenDeformationForEigenStress( initialGuess, targetStress, nullptr );

  Tensor33d fFound( 0.0 );
  fFound( 0, 0 ) = F0_XX;
  fFound( 1, 1 ) = F0_YY;
  fFound( 2, 2 ) = F0_ZZ;

  CompressibleHyperelasticity::Deformation< 3 >          def     = { fFound };
  CompressibleHyperelasticity::TimeIncrement             timeInc = { 0, 0.1 };
  CompressibleHyperelasticity::ConstitutiveResponse< 3 > res( Tensor33d( 0.0 ), 0.0, 0.0, nullptr );
  CompressibleHyperelasticity::AlgorithmicModuli< 3 >    tan = { Tensor3333d( 0.0 ) };
  mat.computeStress( res, tan, def, timeInc );

  throwExceptionOnFailure( checkIfEqual( res.tau( 0, 0 ), std::get< 0 >( targetStress ), 1e-6 ) &&
                             checkIfEqual( res.tau( 1, 1 ), std::get< 1 >( targetStress ), 1e-6 ) &&
                             checkIfEqual( res.tau( 2, 2 ), std::get< 2 >( targetStress ), 1e-6 ),
                           "findEigenDeformationForEigenStress() did not converge to a deformation reproducing "
                           "the target diagonal stress in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// Material properties of every pair of isochoric and volumetric base, with the parameters of both parts.
static std::vector< std::vector< double > > propertiesOfAllBasePairs()
{
  const std::vector< std::vector< double > > isochoric  = { { 0, 1500. },
                                                            { 1, 700., -50., 20. },
                                                            { 2, 500., 200. },
                                                            { 3, 1500., 3. } };
  const std::vector< std::vector< double > > volumetric = { { 1, 3500. },
                                                            { 2, 3500. },
                                                            { 3, 3500. },
                                                            { 4, 3500., -2. },
                                                            { 5, 3500. },
                                                            { 6, 3500., 1.5, 2.5 },
                                                            { 7, 3500. },
                                                            { 8, 3500. } };
  std::vector< std::vector< double > >       properties;
  for ( const auto& iso : isochoric )
    for ( const auto& vol : volumetric ) {
      std::vector< double > p = { iso[0], vol[0] };
      p.insert( p.end(), iso.begin() + 1, iso.end() );
      p.insert( p.end(), vol.begin() + 1, vol.end() );
      properties.push_back( p );
    }
  return properties;
}

static Tensor33d deformedF()
{
  Tensor33d F = Spatial3D::I;
  F( 0, 0 ) += 0.08;
  F( 1, 1 ) -= 0.02;
  F( 2, 2 ) -= 0.015;
  F( 0, 1 ) = 0.03;
  F( 1, 2 ) = -0.02;
  F( 2, 0 ) = 0.01;
  return F;
}

// For every pair of bases: the stress is F ( 2 dPsi/dC ) F^T with the separately evaluated FirstOrder potentials of
// both parts, the energy is the sum of both parts, and the tangent matches central differences of the stress.
void testEveryBasePairMatchesSharedPotentials()
{
  namespace EDF = Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions;

  const Tensor33d F = deformedF();
  const Tensor33d C = Marmot::ContinuumMechanics::Kinematics::DeformationMeasures::rightCauchyGreen( F );

  for ( const auto& props : propertiesOfAllBasePairs() ) {
    const std::string pair = "isochoric base " + std::to_string( static_cast< int >( props[0] ) ) +
                             ", volumetric base " + std::to_string( static_cast< int >( props[1] ) );
    const double* iso = props.data() + 2;
    const double* vol = props.data() + props.size() - ( props[1] == 4 ? 2 : props[1] == 6 ? 3 : 1 );

    std::tuple< double, Tensor33d > isoPart, volPart;
    switch ( static_cast< int >( props[0] ) ) {
    case 0: isoPart = EDF::Isochoric::FirstOrder::NeoHookePotential( C, iso[0] ); break;
    case 1: isoPart = EDF::Isochoric::FirstOrder::YeohPotential( C, iso[0], iso[1], iso[2] ); break;
    case 2: isoPart = EDF::Isochoric::FirstOrder::MooneyRivlinPotential( C, iso[0], iso[1] ); break;
    default: isoPart = EDF::Isochoric::FirstOrder::ArrudaBoyce8ChainPotential( C, iso[0], iso[1] );
    }
    switch ( static_cast< int >( props[1] ) ) {
    case 1: volPart = EDF::Volumetric::FirstOrder::SussmanBathePotential( C, vol[0] ); break;
    case 2: volPart = EDF::Volumetric::FirstOrder::PenaltyPotential( C, vol[0] ); break;
    case 3: volPart = EDF::Volumetric::FirstOrder::SimoTaylorPotential( C, vol[0] ); break;
    case 4: volPart = EDF::Volumetric::FirstOrder::OgdenPotential( C, vol[0], vol[1] ); break;
    case 5: volPart = EDF::Volumetric::FirstOrder::LiuHofstetterMangPotential( C, vol[0] ); break;
    case 6: volPart = EDF::Volumetric::FirstOrder::DollSchweizerhofPotential( C, vol[0], vol[1], vol[2] ); break;
    case 7: volPart = EDF::Volumetric::FirstOrder::DollSchweizerhofExponentialPotential( C, vol[0] ); break;
    default: volPart = EDF::Volumetric::FirstOrder::DollSchweizerhofLogarithmicPotential( C, vol[0] );
    }
    const Tensor33d PK2       = 2. * ( std::get< 1 >( isoPart ) + std::get< 1 >( volPart ) );
    const Tensor33d tauTarget = einsum< iI, IJ, jJ, to_ij >( F, PK2, F );

    const CompressibleHyperelasticity mat( props.data(), static_cast< int >( props.size() ), 1 );
    const auto                        stressAt = [&mat]( const Tensor33d& F_ ) {
      CompressibleHyperelasticity::Deformation< 3 >          def     = { F_ };
      CompressibleHyperelasticity::TimeIncrement             timeInc = { 0, 0.1 };
      CompressibleHyperelasticity::ConstitutiveResponse< 3 > response( Tensor33d( 0.0 ), 0.0, 0.0, nullptr );
      CompressibleHyperelasticity::AlgorithmicModuli< 3 >    tangent = { Tensor3333d( 0.0 ) };
      mat.computeStress( response, tangent, def, timeInc );
      return std::tuple{ response.tau, tangent.dTau_dF, response.elasticEnergyDensity };
    };
    const auto [tau, dTau_dF, psi] = stressAt( F );

    throwExceptionOnFailure( checkIfEqual( tau, tauTarget, 1e-9 ),
                             "stress does not match the shared potentials for " + pair + " in " +
                               std::string( __PRETTY_FUNCTION__ ) );
    throwExceptionOnFailure( checkIfEqual( psi, std::get< 0 >( isoPart ) + std::get< 0 >( volPart ), 1e-10 ),
                             "energy density does not match the shared potentials for " + pair + " in " +
                               std::string( __PRETTY_FUNCTION__ ) );

    const double h = 1e-6;
    for ( int k = 0; k < 3; k++ )
      for ( int l = 0; l < 3; l++ ) {
        Tensor33d Fp = F, Fm = F;
        Fp( k, l ) += h;
        Fm( k, l ) -= h;
        const Tensor33d dTau = ( std::get< 0 >( stressAt( Fp ) ) - std::get< 0 >( stressAt( Fm ) ) ) / ( 2. * h );
        for ( int i = 0; i < 3; i++ )
          for ( int j = 0; j < 3; j++ )
            throwExceptionOnFailure( std::abs( dTau( i, j ) - dTau_dF( i, j, k, l ) ) < 1e-4,
                                     "tangent does not match central differences for " + pair + " in " +
                                       std::string( __PRETTY_FUNCTION__ ) );
      }
  }
}

// NeoHooke with DollSchweizerhof( K, 1, 3 ) is the compressible Neo-Hooke model of Pence & Gou (variant B), the
// former material CompressibleNeoHooke with the properties K, G.
void testNeoHookeDollSchweizerhofIsPenceGouB()
{
  const double                      K = 3500., G = 1500.;
  const std::array< double, 6 >     props = { 0, 6, G, K, 1, 3 };
  const CompressibleHyperelasticity mat( props.data(), 6, 1 );

  const Tensor33d F = deformedF();
  const Tensor33d C = Marmot::ContinuumMechanics::Kinematics::DeformationMeasures::rightCauchyGreen( F );

  CompressibleHyperelasticity::Deformation< 3 >          def     = { F };
  CompressibleHyperelasticity::TimeIncrement             timeInc = { 0, 0.1 };
  CompressibleHyperelasticity::ConstitutiveResponse< 3 > response( Tensor33d( 0.0 ), 0.0, 0.0, nullptr );
  CompressibleHyperelasticity::AlgorithmicModuli< 3 >    tangent = { Tensor3333d( 0.0 ) };
  mat.computeStress( response, tangent, def, timeInc );

  const auto [psi, dPsi_dC] = Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions::FirstOrder::
    PenceGouPotentialB( C, K, G );
  const Tensor33d PK2       = 2. * dPsi_dC;
  const Tensor33d tauTarget = einsum< iI, IJ, jJ, to_ij >( F, PK2, F );

  throwExceptionOnFailure( checkIfEqual( response.tau, tauTarget, 1e-10 ) &&
                             checkIfEqual( response.elasticEnergyDensity, psi, 1e-10 ),
                           "NeoHooke with DollSchweizerhof( K, 1, 3 ) is not PenceGouPotentialB in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// Unknown bases, too few properties and an invalid locking stretch are rejected; the density follows the
// parameters of the potential.
void testMaterialPropertiesAreValidated()
{
  const auto throwsInvalidArgument = []( const std::vector< double >& props ) {
    try {
      CompressibleHyperelasticity mat( props.data(), static_cast< int >( props.size() ), 1 );
    }
    catch ( const std::invalid_argument& ) {
      return true;
    }
    return false;
  };

  throwExceptionOnFailure( throwsInvalidArgument( { 4, 2, 1500., 3500. } ),
                           "unknown isochoric base accepted in " + std::string( __PRETTY_FUNCTION__ ) );
  throwExceptionOnFailure( throwsInvalidArgument( { 0, 0, 1500., 3500. } ) &&
                             throwsInvalidArgument( { 0, 9, 1500., 3500. } ),
                           "unknown volumetric base accepted in " + std::string( __PRETTY_FUNCTION__ ) );
  throwExceptionOnFailure( throwsInvalidArgument( { 0 } ) &&
                             throwsInvalidArgument( { 1, 6, 700., -50., 20., 3500., 1. } ),
                           "too few material properties accepted in " + std::string( __PRETTY_FUNCTION__ ) );
  throwExceptionOnFailure( throwsInvalidArgument( { 3, 2, 1500., 1.0, 3500. } ),
                           "locking stretch <= 1 accepted in " + std::string( __PRETTY_FUNCTION__ ) );
  throwExceptionOnFailure( throwsInvalidArgument( { 0.5, 2, 1500., 3500. } ) &&
                             throwsInvalidArgument( { 0, 5.9999, 1500., 3500. } ) &&
                             throwsInvalidArgument( { -1, 2, 1500., 3500. } ) &&
                             throwsInvalidArgument( { 0, std::numeric_limits< double >::quiet_NaN(), 1500., 3500. } ),
                           "non-integral base selector accepted in " + std::string( __PRETTY_FUNCTION__ ) );
  throwExceptionOnFailure( throwsInvalidArgument( { 0, 2, 1500., 3500., 1., 3. } ) &&
                             throwsInvalidArgument( { 0, 6, 1500., 3500., 1., 3., 2400., 1. } ),
                           "surplus material properties accepted in " + std::string( __PRETTY_FUNCTION__ ) );
  throwExceptionOnFailure( throwsInvalidArgument( { 0, 4, 1500., 3500., 0. } ),
                           "Ogden exponent theta = 0 accepted in " + std::string( __PRETTY_FUNCTION__ ) );
  throwExceptionOnFailure( throwsInvalidArgument( { 0, 6, 1500., 3500., 0., 3. } ) &&
                             throwsInvalidArgument( { 0, 6, 1500., 3500., 1., 1. } ),
                           "invalid Doll-Schweizerhof shape parameters accepted in " +
                             std::string( __PRETTY_FUNCTION__ ) );
  throwExceptionOnFailure( !throwsInvalidArgument( { 3, 2, 1500., 3.0, 3500. } ) &&
                             !throwsInvalidArgument( { 0, 6, 1500., 3500., 1., 3., 2400. } ),
                           "valid material properties rejected in " + std::string( __PRETTY_FUNCTION__ ) );

  const std::vector< double >       withDensity = { 1, 4, 700., -50., 20., 3500., -2., 2400. };
  const CompressibleHyperelasticity mat( withDensity.data(), static_cast< int >( withDensity.size() ), 1 );
  throwExceptionOnFailure( checkIfEqual( mat.getDensity( nullptr ), 2400., 1e-14 ),
                           "density is not the property following the potential parameters in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  const std::vector< double >       withoutDensity = { 1, 4, 700., -50., 20., 3500., -2. };
  const CompressibleHyperelasticity matWithoutDensity( withoutDensity.data(),
                                                       static_cast< int >( withoutDensity.size() ),
                                                       1 );
  bool                              threw = false;
  try {
    matWithoutDensity.getDensity( nullptr );
  }
  catch ( const std::runtime_error& ) {
    threw = true;
  }
  throwExceptionOnFailure( threw,
                           "getDensity() without density did not throw in " + std::string( __PRETTY_FUNCTION__ ) );
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
    testEveryBasePairMatchesSharedPotentials,
    testNeoHookeDollSchweizerhofIsPenceGouB,
    testMaterialPropertiesAreValidated,
  };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
