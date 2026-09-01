#include "Marmot/BergstromBoyce.h"
#include "Marmot/MarmotDeformationMeasures.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include "Marmot/MarmotMaterialPointSolverFiniteStrain.h"
#include "Marmot/MarmotTesting.h"
#include <cmath>
#include <cstring>
#include <string>
#include <vector>

using namespace Marmot::Testing;
using namespace Marmot::Materials;
using namespace Marmot::FastorStandardTensors;
using namespace Marmot::FastorIndices;

// Material properties: muA, kappaA, muB, kappaB, c1, c2, c3, implementationType
static const std::vector< double > PROPS = { 100.0, 1000.0, 50.0, 1000.0, 0.05, 1.0, 1.0, 0.0 };

BergstromBoyce makeMaterial()
{
  return BergstromBoyce( PROPS.data(), PROPS.size(), 1 );
}

std::array< double, 9 > initialState()
{
  std::array< double, 9 > stateVars{};
  BergstromBoyce          mat = makeMaterial();
  mat.initializeYourself( stateVars.data(), 9 );
  return stateVars;
}

Tensor33d combinedNeoHookeStress( const Tensor33d& F, double muSum, double kappaSum )
{
  BergstromBoyce mat = makeMaterial();
  Tensor33d      C   = Marmot::ContinuumMechanics::DeformationMeasures::rightCauchyGreen( F );

  double    psi;
  Tensor33d dPsi_dC;
  std::tie( psi, dPsi_dC ) = mat.neoHookePotential< double >( C, muSum, kappaSum );

  Tensor33d PK2 = 2. * dPsi_dC;
  Tensor33d tau = einsum< iI, IJ, jJ, to_ij >( F, PK2, F );
  return tau;
}

// Test 1: dt -> 0 limit reproduces the combined (unrelaxed) neo-Hookean response,
// with the error shrinking (roughly) linearly as dt shrinks.
void testInstantaneousLimit()
{
  Tensor33d F = Spatial3D::I;
  F( 0, 0 ) += 0.05;
  F( 1, 1 ) -= 0.01;
  F( 2, 2 ) -= 0.01;

  Tensor33d targetStress = combinedNeoHookeStress( F, PROPS[0] + PROPS[2], PROPS[1] + PROPS[3] );

  auto errorAtDt = []( const Tensor33d& F, const Tensor33d& target, double dt ) {
    BergstromBoyce                              mat = makeMaterial();
    std::array< double, 9 >                     sv  = initialState();
    BergstromBoyce::ConstitutiveResponse< 3 >   response( Tensor33d( 0.0 ), 0.0, 0.0, sv.data() );
    BergstromBoyce::AlgorithmicModuli< 3 >      tangent;
    BergstromBoyce::Deformation< 3 >            def{ F };
    BergstromBoyce::TimeIncrement                timeInc{ 0.0, dt };
    mat.computeStress( response, tangent, def, timeInc );
    Tensor33d diff = response.tau - target;
    return sqrt( Fastor::inner( diff, diff ) );
  };

  double e1 = errorAtDt( F, targetStress, 1e-4 );
  double e2 = errorAtDt( F, targetStress, 1e-5 );

  throwExceptionOnFailure( e1 < 1e-2, "I-1: dt->0 limit error at dt=1e-4 too large: " + std::to_string( e1 ) );
  throwExceptionOnFailure( e2 < e1, "I-1: dt->0 limit error did not shrink with smaller dt" );
}

// Test 2: exact degeneracy c1 = 0 -> no flow at all, ever. Stress must equal the
// combined neo-Hookean closed form at every increment of an arbitrary
// ramp+hold history, and Fv must stay at the identity throughout.
void testZeroFlowDegeneracy()
{
  std::vector< double > props0 = PROPS;
  props0[4]                    = 0.0; // c1 = 0
  BergstromBoyce mat( props0.data(), props0.size(), 1 );
  std::array< double, 9 > sv{};
  mat.initializeYourself( sv.data(), 9 );

  std::vector< Tensor33d > Fhistory;
  for ( int i = 1; i <= 5; ++i ) {
    Tensor33d F = Spatial3D::I;
    F( 0, 0 ) += 0.02 * i;
    F( 0, 1 ) = 0.01 * i;
    Fhistory.push_back( F );
  }
  // hold at the last deformation for a few more increments
  for ( int i = 0; i < 3; ++i )
    Fhistory.push_back( Fhistory.back() );

  for ( const auto& F : Fhistory ) {
    BergstromBoyce::ConstitutiveResponse< 3 > response( Tensor33d( 0.0 ), 0.0, 0.0, sv.data() );
    BergstromBoyce::AlgorithmicModuli< 3 >    tangent;
    BergstromBoyce::Deformation< 3 >          def{ F };
    BergstromBoyce::TimeIncrement              timeInc{ 0.0, 10.0 };
    mat.computeStress( response, tangent, def, timeInc );

    Tensor33d target = combinedNeoHookeStress( F, props0[0] + props0[2], props0[1] + props0[3] );

    throwExceptionOnFailure( checkIfEqual( response.tau, target, 1e-10 ),
                             "I-2: c1=0 degeneracy: stress mismatch vs. combined neo-Hookean closed form" );

    Tensor33d Fv( 0.0 );
    memcpy( Fv.data(), sv.data(), 9 * sizeof( double ) );
    throwExceptionOnFailure( checkIfEqual( Fv, Tensor33d( Spatial3D::I ), 1e-10 ),
                             "I-2: c1=0 degeneracy: Fv drifted away from identity" );
  }
}

// Test 3: det(Fe) == det(F) exactly at every increment (isochoric flow invariant),
// and during a long hold the deviatoric stress relaxes toward network A alone.
void testIsochoricFlowAndRelaxationSplit()
{
  BergstromBoyce mat = makeMaterial();
  std::array< double, 9 > sv{};
  mat.initializeYourself( sv.data(), 9 );

  auto step = [&]( const Tensor33d& F, double dt ) {
    BergstromBoyce::ConstitutiveResponse< 3 > response( Tensor33d( 0.0 ), 0.0, 0.0, sv.data() );
    BergstromBoyce::AlgorithmicModuli< 3 >    tangent;
    BergstromBoyce::Deformation< 3 >          def{ F };
    BergstromBoyce::TimeIncrement              timeInc{ 0.0, dt };
    mat.computeStress( response, tangent, def, timeInc );

    Tensor33d Fv( 0.0 );
    memcpy( Fv.data(), sv.data(), 9 * sizeof( double ) );
    Tensor33d Fe = F % Fastor::inverse( Fv );

    double detFe = Fastor::determinant( Fe );
    double detF  = Fastor::determinant( F );
    throwExceptionOnFailure( std::abs( detFe - detF ) < 1e-10,
                             "I-3: det(Fe) != det(F): " + std::to_string( detFe ) + " vs " +
                               std::to_string( detF ) );
    return response.tau;
  };

  // ramp
  Tensor33d Fend = Spatial3D::I;
  Fend( 0, 0 ) += 0.15;
  Fend( 1, 1 ) -= 0.03;
  Fend( 2, 2 ) -= 0.03;

  Tensor33d tauRampEnd( 0.0 );
  for ( int i = 1; i <= 10; ++i ) {
    Tensor33d F = Spatial3D::I + ( double( i ) / 10.0 ) * ( Fend - Spatial3D::I );
    tauRampEnd  = step( F, 1.0 );
  }

  ( void )tauRampEnd;

  // hold
  Tensor33d tauHold( 0.0 );
  for ( int i = 0; i < 30; ++i )
    tauHold = step( Fend, 20.0 );

  Tensor33d devHold           = Marmot::deviatoric( tauHold );
  Tensor33d tauA_alone        = combinedNeoHookeStress( Fend, PROPS[0], PROPS[1] );
  Tensor33d devA_alone        = Marmot::deviatoric( tauA_alone );
  Tensor33d devDiff           = devHold - devA_alone;
  double    relDiff = sqrt( Fastor::inner( devDiff, devDiff ) ) / sqrt( Fastor::inner( devA_alone, devA_alone ) );

  throwExceptionOnFailure( relDiff < 0.05,
                           "I-3: deviatoric stress did not relax toward network-A-alone after long hold, "
                           "relative difference = " +
                             std::to_string( relDiff ) );
}

// Test 4: objectivity under superposed rigid rotation of the current configuration.
void testObjectivity()
{
  BergstromBoyce mat = makeMaterial();

  Tensor33d F_unrotated = Spatial3D::I;
  F_unrotated( 0, 0 )   = 1.05;
  F_unrotated( 0, 1 )   = 0.04;
  F_unrotated( 1, 0 )   = 0.02;
  F_unrotated( 1, 1 )   = 0.98;
  F_unrotated( 2, 2 )   = 0.99;

  std::array< double, 9 > svUnrotated{};
  mat.initializeYourself( svUnrotated.data(), 9 );
  BergstromBoyce::ConstitutiveResponse< 3 > responseU( Tensor33d( 0.0 ), 0.0, 0.0, svUnrotated.data() );
  BergstromBoyce::AlgorithmicModuli< 3 >    tangentU;
  BergstromBoyce::Deformation< 3 >          defU{ F_unrotated };
  BergstromBoyce::TimeIncrement              timeIncU{ 0.0, 5.0 };
  mat.computeStress( responseU, tangentU, defU, timeIncU );
  Tensor33d stressUnrotated = responseU.tau;

  for ( int phi_deg = 0; phi_deg <= 180; phi_deg += 30 ) {
    std::array< double, 9 > sv{};
    mat.initializeYourself( sv.data(), 9 );

    double    phi = Marmot::Math::degToRad( phi_deg );
    Tensor33d Q( 0.0 );
    Q( 0, 0 ) = cos( phi );
    Q( 0, 1 ) = -sin( phi );
    Q( 1, 0 ) = sin( phi );
    Q( 1, 1 ) = cos( phi );
    Q( 2, 2 ) = 1;

    Tensor33d F_rotated = einsum< ik, kj, to_ij >( Q, F_unrotated );

    BergstromBoyce::ConstitutiveResponse< 3 > response( Tensor33d( 0.0 ), 0.0, 0.0, sv.data() );
    BergstromBoyce::AlgorithmicModuli< 3 >    tangent;
    BergstromBoyce::Deformation< 3 >          def{ F_rotated };
    BergstromBoyce::TimeIncrement              timeInc{ 0.0, 5.0 };
    mat.computeStress( response, tangent, def, timeInc );

    Tensor33d stressRotated = einsum< iI, IJ, jJ, to_ij >( Q, stressUnrotated, Q );

    throwExceptionOnFailure( checkIfEqual( response.tau, stressRotated, 1e-9 ),
                             "I-4: objectivity failed at phi_deg=" + std::to_string( phi_deg ) );

    Tensor33d FvUnrotated( 0.0 ), FvRotated( 0.0 );
    memcpy( FvUnrotated.data(), svUnrotated.data(), 9 * sizeof( double ) );
    memcpy( FvRotated.data(), sv.data(), 9 * sizeof( double ) );
    throwExceptionOnFailure( checkIfEqual( FvUnrotated, FvRotated, 1e-9 ),
                             "I-4: Fv is not invariant under superposed spatial rotation" );
  }
}

// Test 5: isotropy under a reference-configuration rotation (single first step
// from a fresh, isotropic initial state).
void testIsotropy()
{
  BergstromBoyce mat = makeMaterial();

  Tensor33d F_unrotated = Spatial3D::I;
  F_unrotated( 0, 0 )   = 1.05;
  F_unrotated( 0, 1 )   = 0.04;
  F_unrotated( 1, 0 )   = 0.02;
  F_unrotated( 1, 1 )   = 0.98;
  F_unrotated( 2, 2 )   = 0.99;

  std::array< double, 9 > svUnrotated{};
  mat.initializeYourself( svUnrotated.data(), 9 );
  BergstromBoyce::ConstitutiveResponse< 3 > responseU( Tensor33d( 0.0 ), 0.0, 0.0, svUnrotated.data() );
  BergstromBoyce::AlgorithmicModuli< 3 >    tangentU;
  BergstromBoyce::Deformation< 3 >          defU{ F_unrotated };
  BergstromBoyce::TimeIncrement              timeIncU{ 0.0, 5.0 };
  mat.computeStress( responseU, tangentU, defU, timeIncU );
  Tensor33d stressUnrotated = responseU.tau;

  for ( int phi_deg = 0; phi_deg <= 180; phi_deg += 30 ) {
    double    phi = Marmot::Math::degToRad( phi_deg );
    Tensor33d R( 0.0 );
    R( 0, 0 ) = cos( phi );
    R( 0, 1 ) = -sin( phi );
    R( 1, 0 ) = sin( phi );
    R( 1, 1 ) = cos( phi );
    R( 2, 2 ) = 1;

    Tensor33d F_rotated = einsum< ik, kj, to_ij >( F_unrotated, R );

    std::array< double, 9 > sv{};
    mat.initializeYourself( sv.data(), 9 );
    BergstromBoyce::ConstitutiveResponse< 3 > response( Tensor33d( 0.0 ), 0.0, 0.0, sv.data() );
    BergstromBoyce::AlgorithmicModuli< 3 >    tangent;
    BergstromBoyce::Deformation< 3 >          def{ F_rotated };
    BergstromBoyce::TimeIncrement              timeInc{ 0.0, 5.0 };
    mat.computeStress( response, tangent, def, timeInc );

    throwExceptionOnFailure( checkIfEqual( response.tau, stressUnrotated, 1e-9 ),
                             "I-5: isotropy failed at phi_deg=" + std::to_string( phi_deg ) );
  }
}

// Test 6: uniaxial tension/relaxation via the material-point solver -- axial
// stress must stay positive and non-decreasing during the ramp, and
// non-increasing (no overshoot) during the subsequent relaxation hold.
void testUniaxialRelaxationWithMPSolver()
{
  using namespace Marmot::Solvers;

  auto        solveropts = MarmotMaterialPointSolverFiniteStrain::SolverOptions();
  std::string matName    = "BERGSTROMBOYCE";
  auto        solver = MarmotMaterialPointSolverFiniteStrain( matName, PROPS.data(), PROPS.size(), solveropts );

  MarmotMaterialPointSolverFiniteStrain::Step rampStep;
  rampStep.gradUIncrementTarget         = Tensor33d( 0.0 );
  rampStep.gradUIncrementTarget( 0, 0 ) = 0.1;
  rampStep.stressIncrementTarget        = Tensor33d( 0.0 );
  rampStep.isGradUComponentControlled          = Tensor33t< bool >( false );
  rampStep.isGradUComponentControlled( 0, 0 )   = true;
  rampStep.isStressComponentControlled          = Tensor33t< bool >( true );
  rampStep.isStressComponentControlled( 0, 0 )  = false;
  rampStep.timeStart = 0.0;
  rampStep.timeEnd   = 1.0;
  rampStep.dTStart   = 0.1;
  rampStep.dTMax     = 0.1;
  rampStep.dTMin     = 0.01;
  solver.addStep( rampStep );

  MarmotMaterialPointSolverFiniteStrain::Step holdStep = rampStep;
  holdStep.gradUIncrementTarget( 0, 0 )                = 0.0;
  holdStep.timeStart                                   = 1.0;
  holdStep.timeEnd                                      = 100.0;
  holdStep.dTStart                                      = 1.0;
  holdStep.dTMax                                        = 10.0;
  holdStep.dTMin                                        = 0.1;
  solver.addStep( holdStep );

  solver.solve();
  auto history = solver.getHistory();

  throwExceptionOnFailure( history.size() > 2, "I-6: MP solver produced too few increments" );

  double prevTime   = -1.0;
  double prevStress = -1.0;
  bool   inHold      = false;
  for ( const auto& h : history ) {
    double axialStress = h.stress( 0, 0 );
    throwExceptionOnFailure( axialStress > -1e-8, "I-6: axial stress went negative in tension" );

    if ( h.time <= 1.0 + 1e-9 ) {
      if ( prevTime >= 0.0 && h.time > prevTime )
        throwExceptionOnFailure( axialStress >= prevStress - 1e-6, "I-6: axial stress decreased during the ramp" );
    }
    else {
      if ( !inHold ) {
        inHold = true;
      }
      else if ( h.time > prevTime ) {
        throwExceptionOnFailure( axialStress <= prevStress + 1e-6,
                                 "I-6: axial stress increased (overshoot) during the relaxation hold" );
      }
    }
    prevTime   = h.time;
    prevStress = axialStress;
  }
}

int main()
{
  auto tests = std::vector< std::function< void() > >{ testInstantaneousLimit,
                                                       testZeroFlowDegeneracy,
                                                       testIsochoricFlowAndRelaxationSplit,
                                                       testObjectivity,
                                                       testIsotropy,
                                                       testUniaxialRelaxationWithMPSolver };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
