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

// Three parameter sets, one per HyperelasticBase, all chosen to represent the
// SAME physical material (equivalent to a combined neo-Hookean network with
// muA=100, muB=50, kappaA=kappaB=1000) whenever the higher-order Yeoh/Mooney-
// Rivlin coefficients are zero -- this lets every base-agnostic physics check
// below be run identically across all three bases against one shared target.
//
// Properties: hyperelasticBase, kappaA, kappaB, A1, A2, A3, B1, B2, B3, c1, c2, c3,
// implementationType
struct PropsVariant {
  std::string          name;
  std::vector< double > props;
};

// ArrudaBoyce's lambdaL (A2/B2) is a modest, physically meaningful locking
// stretch (2.0), not a huge asymptotic-limit value: unlike the OTHER three
// bases here, ArrudaBoyce does not need to numerically approximate a
// DIFFERENT model to pass the shared I-1..I-6 checks below, because
// combinedStress() (see its own comment) evaluates each variant's "combined
// network" target using THAT SAME VARIANT's OWN base -- an EXACT identity
// for any lambdaL (Psi is exactly linear in mu for fixed lambdaL), not an
// approximation that only holds for lambdaL -> infinity. lambdaL=2.0 keeps
// chain stretch safely below the locking limit for every deformation used
// in I-1..I-6 (largest reaches lambdaChain ~= 1.03 in I-3) while still
// exercising the potential's real nonlinearity, unlike an enormous lambdaL
// which would make it numerically indistinguishable from a constant-modulus
// response. The reduction to NeoHooke as lambdaL -> infinity specifically is
// checked separately, at the scalar (isochoric-invariant) level, by
// testArrudaBoyceReducesToNeoHooke below.
static const std::vector< PropsVariant > VARIANTS = {
  { "NeoHooke", { 0, 1000.0, 1000.0, 100.0, 0.0, 0.0, 50.0, 0.0, 0.0, 0.05, 1.0, 1.0, 0.0 } },
  { "Yeoh", { 1, 1000.0, 1000.0, 50.0, 0.0, 0.0, 25.0, 0.0, 0.0, 0.05, 1.0, 1.0, 0.0 } },
  { "MooneyRivlin", { 2, 1000.0, 1000.0, 50.0, 0.0, 0.0, 25.0, 0.0, 0.0, 0.05, 1.0, 1.0, 0.0 } },
  { "ArrudaBoyce", { 3, 1000.0, 1000.0, 100.0, 2.0, 0.0, 50.0, 2.0, 0.0, 0.05, 1.0, 1.0, 0.0 } },
};

BergstromBoyce makeMaterial( const std::vector< double >& props )
{
  return BergstromBoyce( props.data(), props.size(), 1 );
}

std::array< double, 9 > initialState( const std::vector< double >& props )
{
  std::array< double, 9 > stateVars{};
  BergstromBoyce          mat = makeMaterial( props );
  mat.initializeYourself( stateVars.data(), 9 );
  return stateVars;
}

// Combined (single, unrelaxed) network response, evaluated with whichever
// HyperelasticBase the calling variant itself uses -- the physical target
// shared by all four variants, since each represents the same material
// AS SEEN BY ITS OWN BASE. Generalized (base-aware) rather than hardcoded to
// #neoHookePotential specifically: NeoHooke/Yeoh/MooneyRivlin/ArrudaBoyce are
// each linear in their LEADING coefficient (p1: C10, C10, C10, mu) for FIXED
// higher-order shape parameters (p2,p3: C20/C30, C01/unused, lambdaL/unused),
// so combining network A and B by summing p1 (and kappaA+kappaB) and
// evaluating that SAME base ONCE, at network A's own p2/p3, is exactly
// equivalent to evaluating each network separately and adding the stresses --
// PROVIDED network A and B share the same p2/p3 (true of every VARIANTS
// entry above, by construction: Yeoh/MooneyRivlin's higher-order slots are
// both zero, ArrudaBoyce's lambdaL_A==lambdaL_B==2.0). p2/p3 themselves must
// NOT be summed (unlike p1): they are nonlinear shape parameters, not moduli
// -- e.g. summing two equal lambdaL's would double the locking stretch, not
// preserve it. This matters specifically for ArrudaBoyce: its isochoric
// invariant Ibar1=I1*J^(-2/3) is a genuinely different tensor field from the
// other three bases' raw I1 away from J=1 (dIbar1/dC != dI1/dC in general,
// even though both potentials happen to coincide at C=I) -- comparing
// ArrudaBoyce's full compressible response against a raw-I1-based target
// would NOT converge as lambdaL grows, since the mismatch is structural, not
// a residual of the asymptotic reduction.
Tensor33d combinedStress( const Tensor33d& F, int base, double p1, double p2, double p3, double kappa )
{
  BergstromBoyce mat = makeMaterial( VARIANTS[0].props );
  Tensor33d      C   = Marmot::ContinuumMechanics::DeformationMeasures::rightCauchyGreen( F );

  double    psi;
  Tensor33d dPsi_dC;
  std::tie( psi, dPsi_dC ) = mat.hyperelasticPotential< double >( C, base, p1, p2, p3, kappa );

  Tensor33d PK2 = 2. * dPsi_dC;
  Tensor33d tau = einsum< iI, IJ, jJ, to_ij >( F, PK2, F );
  return tau;
}

// ---------------------------------------------------------------------------
// Potential-level checks: stress-free reference state, and the two exact
// algebraic reductions (Yeoh -> NeoHooke, MooneyRivlin -> NeoHooke) when the
// higher-order coefficients vanish.
// ---------------------------------------------------------------------------

std::vector< Tensor33d > randomSymmetricPositiveDefiniteCs()
{
  std::vector< Tensor33d > Cs;
  auto makeC = [&]( double a, double b, double c, double d, double e, double f ) {
    Tensor33d F = Spatial3D::I;
    F( 0, 0 ) += a;
    F( 1, 1 ) += b;
    F( 2, 2 ) += c;
    F( 0, 1 ) = d;
    F( 1, 2 ) = e;
    F( 0, 2 ) = f;
    return Marmot::ContinuumMechanics::DeformationMeasures::rightCauchyGreen( F );
  };
  Cs.push_back( makeC( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ) );        // C = I
  Cs.push_back( makeC( 0.1, -0.03, -0.03, 0.02, 0.01, -0.01 ) ); // moderate, mildly triaxial
  Cs.push_back( makeC( 0.3, 0.1, -0.2, 0.05, -0.03, 0.04 ) );    // larger, fully general
  return Cs;
}

void testPotentialStressFreeAtReference()
{
  BergstromBoyce mat = makeMaterial( VARIANTS[0].props );
  Tensor33d       C   = Spatial3D::I;

  for ( const auto& variant : VARIANTS ) {
    BergstromBoyce m2 = makeMaterial( variant.props );
    double         psi;
    Tensor33d      dPsi_dC;
    std::tie( psi, dPsi_dC ) = m2.hyperelasticPotential< double >( C, variant.props[0], 50.0, 5.0, 0.3, 1000.0 );
    throwExceptionOnFailure( checkIfEqual( dPsi_dC, Tensor33d( 0.0 ), 1e-10 ),
                             "P-1 [" + variant.name + "]: potential is not stress-free at C=I" );
  }
}

void testYeohReducesToNeoHooke()
{
  BergstromBoyce mat = makeMaterial( VARIANTS[0].props );
  const double    mu = 42.0, kappa = 777.0;

  for ( const auto& C : randomSymmetricPositiveDefiniteCs() ) {
    double    psiNH;
    Tensor33d dPsiNH_dC;
    std::tie( psiNH, dPsiNH_dC ) = mat.neoHookePotential< double >( C, mu, kappa );

    double    psiYeoh;
    Tensor33d dPsiYeoh_dC;
    std::tie( psiYeoh, dPsiYeoh_dC ) = mat.yeohPotential< double >( C, mu / 2., 0.0, 0.0, kappa );

    throwExceptionOnFailure( std::abs( psiYeoh - psiNH ) < 1e-10 * ( 1. + std::abs( psiNH ) ),
                             "P-2: Yeoh(C20=C30=0,C10=mu/2) energy does not match NeoHooke(mu)" );
    throwExceptionOnFailure( checkIfEqual( dPsiYeoh_dC, dPsiNH_dC, 1e-10 ),
                             "P-2: Yeoh(C20=C30=0,C10=mu/2) stress does not match NeoHooke(mu)" );
  }
}

void testMooneyRivlinReducesToNeoHooke()
{
  BergstromBoyce mat = makeMaterial( VARIANTS[0].props );
  const double    mu = 42.0, kappa = 777.0;

  for ( const auto& C : randomSymmetricPositiveDefiniteCs() ) {
    double    psiNH;
    Tensor33d dPsiNH_dC;
    std::tie( psiNH, dPsiNH_dC ) = mat.neoHookePotential< double >( C, mu, kappa );

    double    psiMR;
    Tensor33d dPsiMR_dC;
    std::tie( psiMR, dPsiMR_dC ) = mat.mooneyRivlinPotential< double >( C, mu / 2., 0.0, kappa );

    throwExceptionOnFailure( std::abs( psiMR - psiNH ) < 1e-10 * ( 1. + std::abs( psiNH ) ),
                             "P-3: MooneyRivlin(C01=0,C10=mu/2) energy does not match NeoHooke(mu)" );
    throwExceptionOnFailure( checkIfEqual( dPsiMR_dC, dPsiNH_dC, 1e-10 ),
                             "P-3: MooneyRivlin(C01=0,C10=mu/2) stress does not match NeoHooke(mu)" );
  }
}

// ArrudaBoyce reduces to NeoHooke only in the lambdaL -> infinity limit, and
// ONLY as a scalar function of the isochoric invariant Ibar1
// (Psi_iso(Ibar1) -> (mu/2)(Ibar1-3)) -- NOT as a full tensorial dPsi/dC
// match against BergstromBoyce's own raw-I1-based #neoHookePotential, which
// is a genuinely different tensor field away from C=I: dIbar1/dC != dI1/dC
// in general (they only coincide, trivially, at C=I itself), so a full
// compressible-gradient comparison against #neoHookePotential would NOT
// converge as lambdaL grows (confirmed empirically while developing this
// test -- the error stayed ~0.28 from lambdaL=1e2 all the way to lambdaL=1e6).
// This instead checks the actual, correct claim directly against the shared
// core-header function, at purely isochoric C (det C=1, so Ibar1=I1
// exactly), which isolates it from BergstromBoyce's own volumetric
// convention entirely.
void testArrudaBoyceReducesToNeoHooke()
{
  const double mu = 42.0;

  auto isochoricC = []( double a ) {
    Tensor33d F = Spatial3D::I;
    F( 0, 0 )   = a;
    F( 1, 1 )   = 1.0 / sqrt( a );
    F( 2, 2 )   = 1.0 / sqrt( a );
    return Marmot::ContinuumMechanics::DeformationMeasures::rightCauchyGreen( F );
  };

  // deliberately excludes a=1 (C=I): there, Psi_iso is identically zero for
  // ANY lambdaL (per the exact stress-free-at-reference proof behind P-1),
  // making an "error shrinks as lambdaL grows" comparison degenerate.
  for ( double a : { 1.2, 0.8, 1.6 } ) {
    Tensor33d    C           = isochoricC( a );
    const double I1         = Fastor::trace( C );
    const double psiTargetNH = mu / 2. * ( I1 - 3. );

    auto errorAtLambdaL = [&]( double lambdaL ) {
      double psiAB =
        Marmot::ContinuumMechanics::EnergyDensityFunctions::ArrudaBoyce8ChainPotential< double >( C, mu, lambdaL );
      return std::abs( psiAB - psiTargetNH );
    };

    // The true asymptotic error is O(1/lambdaL^2) (verified analytically and
    // by direct numerical experiment while developing this test): it shrinks
    // cleanly from lambdaL=1e2 through 1e4, but the naive log(w/w0) formula
    // suffers catastrophic cancellation once w and w0 both sit within ~1e-8
    // of 1 (lambdaL >~ 1e4 here), so the raw error GROWS again from 1e5 to
    // 1e6 -- a floating-point-precision artifact of this test's direct
    // double-precision evaluation, not a defect in the potential itself (no
    // physically meaningful locking stretch approaches lambdaL=1e6 in
    // practice). Compare 1e2 vs. 1e4, which stays inside the clean
    // convergence regime.
    double eSmall = errorAtLambdaL( 1e2 );
    double eLarge = errorAtLambdaL( 1e4 );

    throwExceptionOnFailure( eLarge < eSmall,
                             "P-4: ArrudaBoyce isochoric energy error vs. (mu/2)(Ibar1-3) did not shrink as lambdaL "
                             "grew (1e2: " +
                               std::to_string( eSmall ) + ", 1e4: " + std::to_string( eLarge ) + ")" );
    throwExceptionOnFailure( eLarge < 1e-3 * ( 1. + std::abs( psiTargetNH ) ),
                             "P-4: ArrudaBoyce(lambdaL=1e4) isochoric energy does not converge to (mu/2)(Ibar1-3): " +
                               std::to_string( eLarge ) );
  }
}

// Numerically verify dI2/dC = I1*I - C, via central finite differences of a
// standalone I2(C) = 1/2*(I1^2 - tr(C^2)) implementation -- do not just trust
// the closed-form derivative used inside mooneyRivlinPotential.
double I2standalone( const Tensor33d& C )
{
  const double I1  = Fastor::trace( C );
  Tensor33d    CSq = einsum< IK, KJ, to_ij >( C, C );
  const double trCSq = Fastor::trace( CSq );
  return 0.5 * ( I1 * I1 - trCSq );
}

void testI2DerivativeIdentity()
{
  const double h = 1e-6;

  for ( const auto& C : randomSymmetricPositiveDefiniteCs() ) {
    const double I1 = Fastor::trace( C );
    Tensor33d    analyticaldI2_dC = Marmot::multiplyFastorTensorWithScalar( Tensor33d( Spatial3D::I ), I1 ) - C;

    // NOTE: perturb only the (i,j) entry, leaving (j,i) at its original value -- this
    // reproduces the naive component-wise partial derivative that "I1*I - C" represents.
    // A *symmetrized* perturbation (also bumping (j,i) by the same amount) double-counts
    // the off-diagonal sensitivity and gives 2x the correct off-diagonal entries -- this
    // was caught by this very check on the first attempt (see git history / fork report).
    Tensor33d fdDI2_dC( 0.0 );
    for ( int i = 0; i < 3; ++i ) {
      for ( int j = 0; j < 3; ++j ) {
        Tensor33d Cp = C, Cm = C;
        Cp( i, j ) += h;
        Cm( i, j ) -= h;
        fdDI2_dC( i, j ) = ( I2standalone( Cp ) - I2standalone( Cm ) ) / ( 2. * h );
      }
    }

    throwExceptionOnFailure( checkIfEqual( fdDI2_dC, analyticaldI2_dC, 1e-4 ),
                             "P-4: dI2/dC = I1*I - C identity failed finite-difference check" );
  }
}

// ---------------------------------------------------------------------------
// Full-model, base-agnostic physics checks, run once per HyperelasticBase
// variant against the shared combined-neo-Hookean target.
// ---------------------------------------------------------------------------

// Test 1: dt -> 0 limit reproduces the combined (unrelaxed) neo-Hookean response,
// with the error shrinking (roughly) linearly as dt shrinks.
void testInstantaneousLimit( const PropsVariant& variant )
{
  Tensor33d F = Spatial3D::I;
  F( 0, 0 ) += 0.05;
  F( 1, 1 ) -= 0.01;
  F( 2, 2 ) -= 0.01;

  const auto& p            = variant.props;
  Tensor33d   targetStress = combinedStress( F, int( p[0] ), p[3] + p[6], p[4], p[5], p[1] + p[2] );

  auto errorAtDt = [&]( const Tensor33d& F, const Tensor33d& target, double dt ) {
    BergstromBoyce                              mat = makeMaterial( variant.props );
    std::array< double, 9 >                     sv  = initialState( variant.props );
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

  throwExceptionOnFailure( e1 < 1e-2,
                           "I-1 [" + variant.name + "]: dt->0 limit error at dt=1e-4 too large: " +
                             std::to_string( e1 ) );
  throwExceptionOnFailure( e2 < e1, "I-1 [" + variant.name + "]: dt->0 limit error did not shrink with smaller dt" );
}

// Test 2: exact degeneracy c1 = 0 -> no flow at all, ever. Stress must equal the
// combined neo-Hookean closed form at every increment of an arbitrary
// ramp+hold history, and Fv must stay at the identity throughout.
void testZeroFlowDegeneracy( const PropsVariant& variant )
{
  std::vector< double > props0 = variant.props;
  props0[9]                    = 0.0; // c1 = 0
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

    const auto& p      = variant.props;
    Tensor33d   target = combinedStress( F, int( p[0] ), p[3] + p[6], p[4], p[5], p[1] + p[2] );

    throwExceptionOnFailure( checkIfEqual( response.tau, target, 1e-10 ),
                             "I-2 [" + variant.name +
                               "]: c1=0 degeneracy: stress mismatch vs. combined neo-Hookean closed form" );

    Tensor33d Fv( 0.0 );
    memcpy( Fv.data(), sv.data(), 9 * sizeof( double ) );
    throwExceptionOnFailure( checkIfEqual( Fv, Tensor33d( Spatial3D::I ), 1e-10 ),
                             "I-2 [" + variant.name + "]: c1=0 degeneracy: Fv drifted away from identity" );
  }
}

// Test 3: det(Fe) == det(F) exactly at every increment (isochoric flow invariant),
// and during a long hold the deviatoric stress relaxes toward network A alone.
void testIsochoricFlowAndRelaxationSplit( const PropsVariant& variant )
{
  BergstromBoyce mat = makeMaterial( variant.props );
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
                             "I-3 [" + variant.name + "]: det(Fe) != det(F): " + std::to_string( detFe ) + " vs " +
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

  const auto& p = variant.props;

  Tensor33d devHold    = Marmot::deviatoric( tauHold );
  Tensor33d tauA_alone = combinedStress( Fend, int( p[0] ), p[3], p[4], p[5], p[1] );
  Tensor33d devA_alone = Marmot::deviatoric( tauA_alone );
  Tensor33d devDiff    = devHold - devA_alone;
  double    relDiff = sqrt( Fastor::inner( devDiff, devDiff ) ) / sqrt( Fastor::inner( devA_alone, devA_alone ) );

  throwExceptionOnFailure( relDiff < 0.05,
                           "I-3 [" + variant.name +
                             "]: deviatoric stress did not relax toward network-A-alone after long hold, "
                             "relative difference = " +
                             std::to_string( relDiff ) );
}

// Test 4: objectivity under superposed rigid rotation of the current configuration.
void testObjectivity( const PropsVariant& variant )
{
  BergstromBoyce mat = makeMaterial( variant.props );

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
                             "I-4 [" + variant.name + "]: objectivity failed at phi_deg=" +
                               std::to_string( phi_deg ) );

    Tensor33d FvUnrotated( 0.0 ), FvRotated( 0.0 );
    memcpy( FvUnrotated.data(), svUnrotated.data(), 9 * sizeof( double ) );
    memcpy( FvRotated.data(), sv.data(), 9 * sizeof( double ) );
    throwExceptionOnFailure( checkIfEqual( FvUnrotated, FvRotated, 1e-9 ),
                             "I-4 [" + variant.name + "]: Fv is not invariant under superposed spatial rotation" );
  }
}

// Test 5: isotropy under a reference-configuration rotation (single first step
// from a fresh, isotropic initial state).
void testIsotropy( const PropsVariant& variant )
{
  BergstromBoyce mat = makeMaterial( variant.props );

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
                             "I-5 [" + variant.name + "]: isotropy failed at phi_deg=" + std::to_string( phi_deg ) );
  }
}

// Test 6: uniaxial tension/relaxation via the material-point solver -- axial
// stress must stay positive and non-decreasing during the ramp, and
// non-increasing (no overshoot) during the subsequent relaxation hold.
void testUniaxialRelaxationWithMPSolver( const PropsVariant& variant )
{
  using namespace Marmot::Solvers;

  auto        solveropts = MarmotMaterialPointSolverFiniteStrain::SolverOptions();
  std::string matName    = "BERGSTROMBOYCE";
  auto solver = MarmotMaterialPointSolverFiniteStrain( matName, variant.props.data(), variant.props.size(),
                                                       solveropts );

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

  throwExceptionOnFailure( history.size() > 2, "I-6 [" + variant.name + "]: MP solver produced too few increments" );

  double prevTime   = -1.0;
  double prevStress = -1.0;
  bool   inHold      = false;
  for ( const auto& h : history ) {
    double axialStress = h.stress( 0, 0 );
    throwExceptionOnFailure( axialStress > -1e-8, "I-6 [" + variant.name + "]: axial stress went negative in tension" );

    if ( h.time <= 1.0 + 1e-9 ) {
      if ( prevTime >= 0.0 && h.time > prevTime )
        throwExceptionOnFailure( axialStress >= prevStress - 1e-6,
                                 "I-6 [" + variant.name + "]: axial stress decreased during the ramp" );
    }
    else {
      if ( !inHold ) {
        inHold = true;
      }
      else if ( h.time > prevTime ) {
        throwExceptionOnFailure( axialStress <= prevStress + 1e-6,
                                 "I-6 [" + variant.name +
                                   "]: axial stress increased (overshoot) during the relaxation hold" );
      }
    }
    prevTime   = h.time;
    prevStress = axialStress;
  }
}

int main()
{
  std::vector< std::function< void() > > tests = {
    testPotentialStressFreeAtReference,
    testYeohReducesToNeoHooke,
    testMooneyRivlinReducesToNeoHooke,
    testArrudaBoyceReducesToNeoHooke,
    testI2DerivativeIdentity,
  };

  for ( const auto& variant : VARIANTS ) {
    tests.push_back( [variant]() { testInstantaneousLimit( variant ); } );
    tests.push_back( [variant]() { testZeroFlowDegeneracy( variant ); } );
    tests.push_back( [variant]() { testIsochoricFlowAndRelaxationSplit( variant ); } );
    tests.push_back( [variant]() { testObjectivity( variant ); } );
    tests.push_back( [variant]() { testIsotropy( variant ); } );
    tests.push_back( [variant]() { testUniaxialRelaxationWithMPSolver( variant ); } );
  }

  executeTestsAndCollectExceptions( tests );

  return 0;
}
