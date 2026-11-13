#include "Fastor/Fastor.h"
#include "Marmot/MarmotEnergyDensityFunctions.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotTesting.h"
#include <complex>
#include <stdexcept>
#include <vector>

using namespace Marmot::Testing;
using namespace Marmot::Math::Tensors::FastorTensors;

std::tuple< Tensor33d, Tensor33d, double, double, double, double > computationParameters()
{
  // define material properties
  const double K = 10000.0;
  const double G = 1000.0;
  // define right Cauchy-Green tensor
  Tensor33d C;
  C( 0, 0 ) = 1.25;
  C( 1, 1 ) = 2.0;
  C( 2, 2 ) = 2.0;
  C( 0, 1 ) = 0.5;
  C( 0, 2 ) = 0.5;
  C( 1, 0 ) = 0.5;
  C( 1, 2 ) = 0.0;
  C( 2, 0 ) = 0.5;
  C( 2, 1 ) = 0.0;
  // define inverse C
  Tensor33d invC;
  invC( 0, 0 ) = 1.0;
  invC( 1, 1 ) = 0.5625;
  invC( 2, 2 ) = 0.5625;
  invC( 0, 1 ) = -0.25;
  invC( 0, 2 ) = -0.25;
  invC( 1, 0 ) = -0.25;
  invC( 1, 2 ) = 0.0625;
  invC( 2, 0 ) = -0.25;
  invC( 2, 1 ) = 0.0625;
  // define J and I1
  const double J  = 2.0;
  const double I1 = 5.25;
  return { C, invC, K, G, J, I1 };
}

auto testPenceGouPotentialA()
{
  using namespace Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions;
  // get computation parameters
  std::tuple< Tensor33d, Tensor33d, double, double, double, double > params = computationParameters();
  Tensor33d                                                          C      = get< 0 >( params );
  const double                                                       K      = get< 2 >( params );
  const double                                                       G      = get< 3 >( params );
  // compute expected and resulting value and compare
  double expected = 5098.519486106722;
  double res      = PenceGouPotentialA( C, K, G );
  throwExceptionOnFailure( checkIfEqual( res, expected, 1e-15 ),
                           Marmot::Utilities::MakeString()
                             << __PRETTY_FUNCTION__ << " energy density A function failed." );
}

auto testPenceGouPotentialB()
{
  using namespace Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions;
  // get computation parameters
  std::tuple< Tensor33d, Tensor33d, double, double, double, double > params = computationParameters();
  Tensor33d                                                          C      = get< 0 >( params );
  const double                                                       K      = get< 2 >( params );
  const double                                                       G      = get< 3 >( params );
  // compute expected and resulting value and compare
  double expected = 2966.146377987021;
  double res      = PenceGouPotentialB( C, K, G );
  throwExceptionOnFailure( checkIfEqual( res, expected, 1e-15 ),
                           Marmot::Utilities::MakeString()
                             << __PRETTY_FUNCTION__ << " energy density B function failed." );
}

auto testPenceGouPotentialC()
{
  using namespace Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions;
  // get computation parameters
  std::tuple< Tensor33d, Tensor33d, double, double, double, double > params = computationParameters();
  Tensor33d                                                          C      = get< 0 >( params );
  const double                                                       K      = get< 2 >( params );
  const double                                                       G      = get< 3 >( params );
  // compute expected and resulting value and compare
  double expected = 1018.0232353221228;
  double res      = PenceGouPotentialC( C, K, G );
  throwExceptionOnFailure( checkIfEqual( res, expected, 1e-15 ),
                           Marmot::Utilities::MakeString()
                             << __PRETTY_FUNCTION__ << " energy density C function failed." );
}

auto testFirstOrderDerivedB()
{
  using namespace Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions::FirstOrder;
  // get computation parameters
  std::tuple< Tensor33d, Tensor33d, double, double, double, double > params = computationParameters();
  Tensor33d                                                          C      = get< 0 >( params );
  Tensor33d                                                          invC   = get< 1 >( params );
  const double                                                       K      = get< 2 >( params );
  const double                                                       G      = get< 3 >( params );
  const double                                                       J      = get< 4 >( params );
  const double                                                       I1     = get< 5 >( params );
  // compute expected and resulting values
  double                          expectedED = 2966.146377987021;
  std::tuple< double, Tensor33d > res        = PenceGouPotentialB( C, K, G );
  double                          resED      = get< 0 >( res );
  Tensor33d                       resS       = get< 1 >( res );
  // compare energy density
  throwExceptionOnFailure( checkIfEqual( resED, expectedED, 1e-15 ),
                           Marmot::Utilities::MakeString()
                             << __PRETTY_FUNCTION__ << " energy density B function failed." );
  // constant parameters
  const double dW_dI1 = G / 2. * pow( J, -2. / 3 );
  const double dW_dJ  = K / 4. * ( J - pow( J, -3. ) ) - G * I1 / 3. * pow( J, -5. / 3 );
  double       expectedS;
  // compare stress
  for ( int i = 0; i < 3; i++ )
    for ( int j = 0; j < 3; j++ ) {
      expectedS = dW_dI1 * ( i == j ? 1 : 0 ) + dW_dJ * J / 2. * invC( i, j );
      throwExceptionOnFailure( checkIfEqual( resS( i, j ), expectedS, 1e-12 ),
                               Marmot::Utilities::MakeString()
                                 << __PRETTY_FUNCTION__ << " stress computation B failed for position (" << i << ", "
                                 << j << ")" );
    }
}

auto testSecondOrderDerivedB()
{
  using namespace Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions::SecondOrder;
  // get computation parameters
  std::tuple< Tensor33d, Tensor33d, double, double, double, double > params = computationParameters();
  Tensor33d                                                          C      = get< 0 >( params );
  Tensor33d                                                          invC   = get< 1 >( params );
  const double                                                       K      = get< 2 >( params );
  const double                                                       G      = get< 3 >( params );
  const double                                                       J      = get< 4 >( params );
  const double                                                       I1     = get< 5 >( params );
  // compute expected and resulting values
  double                                                 expectedED = 2966.146377987021;
  std::tuple< double, Tensor33d, Tensor3333t< double > > res        = PenceGouPotentialB( C, K, G );
  double                                                 resED      = get< 0 >( res );
  Tensor33d                                              resS       = get< 1 >( res );
  Tensor3333t< double >                                  resCSE     = get< 2 >( res );
  // compare energy density
  throwExceptionOnFailure( checkIfEqual( resED, expectedED, 1e-15 ),
                           Marmot::Utilities::MakeString()
                             << __PRETTY_FUNCTION__ << " energy density B function failed." );
  // constant parameters
  const double dW_dI1    = G / 2. * pow( J, -2. / 3 );
  const double dW_dJ     = K / 8. * ( pow( J, 2. ) - pow( J, -2. ) ) - G * I1 / 6. * pow( J, -2. / 3 );
  const double d2W_dI1dJ = -G / 3. * pow( J, -5. / 3 );
  const double d2W_dJ2   = K / 4. * ( 1. + 3. * pow( J, -4. ) ) + I1 * G * 5. / 9. * pow( J, -8. / 3 );
  double       expectedS;
  double       expectedCSE;
  // compare stress
  for ( int i = 0; i < 3; i++ )
    for ( int j = 0; j < 3; j++ ) {
      expectedS = dW_dI1 * ( i == j ? 1 : 0 ) + dW_dJ * J / 2. * invC( i, j );
      throwExceptionOnFailure( checkIfEqual( resS( i, j ), expectedS, 1e-12 ),
                               Marmot::Utilities::MakeString()
                                 << __PRETTY_FUNCTION__ << " stress computation B failed for position (" << i << ", "
                                 << j << ")" );
      for ( int k = 0; k < 3; k++ )
        for ( int l = 0; l < 3; l++ ) {
          expectedCSE = ( J * J ) / 4 * d2W_dJ2 * invC( j, i ) * invC( l, k ) +
                        dW_dJ * ( J / 4. * invC( j, i ) * invC( l, k ) - J / 2. * invC( j, k ) * invC( l, i ) ) +
                        d2W_dI1dJ * J / 2. * ( invC( j, i ) * ( k == l ? 1 : 0 ) + ( i == j ? 1 : 0 ) * invC( l, k ) );
          throwExceptionOnFailure( checkIfEqual( resCSE( i, j, k, l ), expectedCSE, 1e-12 ),
                                   Marmot::Utilities::MakeString()
                                     << __PRETTY_FUNCTION__ << " tangent computation B failed for position (" << i
                                     << ", " << j << ", " << k << ", " << l << ")" );
        }
    }
}

// Cross-check that the plain (energy-only) and FirstOrder (energy + analytic
// dPsi/dC) overloads of ArrudaBoyce8ChainPotential cannot drift out of sync: both
// funnel through the same internal `Detail::arrudaBoyce8ChainEnergyAndDerivative`
// helper, but this test verifies that end-to-end rather than assuming it. The
// gradient is cross-checked against a complex-step derivative of the plain overload
// itself (not against a hand-derived formula), so this specifically guards against
// the two overloads' shared helper silently diverging from what the plain overload
// actually computes.
auto testArrudaBoycePlainMatchesFirstOrder()
{
  using namespace Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions;

  Tensor33d    C       = get< 0 >( computationParameters() );
  const double mu      = 500.0;
  const double lambdaL = 2.5;

  const double psiPlain = ArrudaBoyce8ChainPotential< double >( C, mu, lambdaL );

  double    psiFOD;
  Tensor33d dPsiFOD_dC;
  std::tie( psiFOD, dPsiFOD_dC ) = FirstOrder::ArrudaBoyce8ChainPotential< double >( C, mu, lambdaL );

  throwExceptionOnFailure( checkIfEqual( psiPlain, psiFOD, 1e-14 ),
                           Marmot::Utilities::MakeString()
                             << __PRETTY_FUNCTION__ << " plain and FirstOrder energy densities disagree." );

  // complex-step derivative of the plain overload, entry by entry (C treated as a
  // general, not-necessarily-symmetric matrix argument, matching how dPsi/dC is
  // defined and used everywhere else in this codebase)
  const double h = 1e-20;
  Tensor33d    dPsiCSDA( 0.0 );
  for ( int i = 0; i < 3; i++ )
    for ( int j = 0; j < 3; j++ ) {
      Tensor33t< std::complex< double > >
        Cc = Marmot::Math::Tensors::fastorTensorFromDoubleTensor< std::complex< double > >( C );
      Cc( i, j ) += std::complex< double >( 0.0, h );
      std::complex< double > psiC = ArrudaBoyce8ChainPotential< std::complex< double > >( Cc, mu, lambdaL );
      dPsiCSDA( i, j )            = psiC.imag() / h;
    }

  throwExceptionOnFailure( checkIfEqual( dPsiFOD_dC, dPsiCSDA, 1e-9 ),
                           Marmot::Utilities::MakeString()
                             << __PRETTY_FUNCTION__
                             << " FirstOrder gradient disagrees with a complex-step derivative "
                                "of the plain overload." );
}

// Cross-check that Arruda-Boyce returns exactly zero energy and zero derivative at
// C=I even when the locking stretch lambdaL is arbitrarily close to 1.0 (where
// w0 = 1 - 1/lambdaL^2 is very small and an unscaled floor would clamp the reference state).
auto testArrudaBoyceZeroEnergyAtReferenceForNearOneLockingStretch()
{
  using namespace Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions;
  const Tensor33d C  = Marmot::Math::Tensors::fastorTensorFromDoubleTensor< double >( Spatial3D::I );
  const double    mu = 500.0;

  for ( double lambdaL : { 1.00001, 1.0001, 1.001, 1.005, 1.01, 2.5 } ) {
    const double psiPlain = ArrudaBoyce8ChainPotential< double >( C, mu, lambdaL );
    double       psiFOD;
    Tensor33d    dPsiFOD_dC;
    std::tie( psiFOD, dPsiFOD_dC ) = FirstOrder::ArrudaBoyce8ChainPotential< double >( C, mu, lambdaL );

    throwExceptionOnFailure( std::abs( psiPlain ) < 1e-12,
                             Marmot::Utilities::MakeString()
                               << __PRETTY_FUNCTION__
                               << " plain Arruda-Boyce energy non-zero at C=I for lambdaL=" << lambdaL );
    throwExceptionOnFailure( std::abs( psiFOD ) < 1e-12,
                             Marmot::Utilities::MakeString()
                               << __PRETTY_FUNCTION__
                               << " FirstOrder Arruda-Boyce energy non-zero at C=I for lambdaL=" << lambdaL );
    throwExceptionOnFailure( checkIfEqual( dPsiFOD_dC, Tensor33d( 0.0 ), 1e-10 ),
                             Marmot::Utilities::MakeString()
                               << __PRETTY_FUNCTION__
                               << " FirstOrder Arruda-Boyce gradient non-zero at C=I for lambdaL=" << lambdaL );
  }
}

// Generic plain-vs-FirstOrder cross-check, same methodology as
// testArrudaBoycePlainMatchesFirstOrder: the FirstOrder gradient
// is checked against a complex-step derivative of the plain overload itself,
// guarding against the two overloads' formulas silently drifting apart.
template < typename PlainFn, typename FirstOrderFn, typename... Args >
void checkPlainMatchesFirstOrder( const char*      testName,
                                  PlainFn          plainFn,
                                  FirstOrderFn     firstOrderFn,
                                  const Tensor33d& C,
                                  Args... args )
{
  const double psiPlain = plainFn( C, args... );

  double    psiFOD;
  Tensor33d dPsiFOD_dC;
  std::tie( psiFOD, dPsiFOD_dC ) = firstOrderFn( C, args... );

  throwExceptionOnFailure( checkIfEqual( psiPlain, psiFOD, 1e-14 ),
                           Marmot::Utilities::MakeString()
                             << testName << " plain and FirstOrder energy densities disagree." );

  const double h = 1e-20;
  Tensor33d    dPsiCSDA( 0.0 );
  for ( int i = 0; i < 3; i++ )
    for ( int j = 0; j < 3; j++ ) {
      Tensor33t< std::complex< double > >
        Cc = Marmot::Math::Tensors::fastorTensorFromDoubleTensor< std::complex< double > >( C );
      Cc( i, j ) += std::complex< double >( 0.0, h );
      std::complex< double > psiC = plainFn( Cc, args... );
      dPsiCSDA( i, j )            = psiC.imag() / h;
    }

  throwExceptionOnFailure( checkIfEqual( dPsiFOD_dC, dPsiCSDA, 1e-9 ),
                           Marmot::Utilities::MakeString()
                             << testName
                             << " FirstOrder gradient disagrees with a complex-step derivative "
                                "of the plain overload." );
}

auto testNeoHookePlainMatchesFirstOrder()
{
  using namespace Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions;
  Tensor33d C = get< 0 >( computationParameters() );
  checkPlainMatchesFirstOrder(
    __PRETTY_FUNCTION__,
    []( const auto& C, double mu ) { return NeoHookePotential( C, mu ); },
    []( const auto& C, double mu ) { return FirstOrder::NeoHookePotential( C, mu ); },
    C,
    750.0 );
}

auto testYeohPlainMatchesFirstOrder()
{
  using namespace Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions;
  Tensor33d C = get< 0 >( computationParameters() );
  checkPlainMatchesFirstOrder(
    __PRETTY_FUNCTION__,
    []( const auto& C, double c10, double c20, double c30 ) { return YeohPotential( C, c10, c20, c30 ); },
    []( const auto& C, double c10, double c20, double c30 ) { return FirstOrder::YeohPotential( C, c10, c20, c30 ); },
    C,
    500.0,
    50.0,
    5.0 );
}

auto testMooneyRivlinPlainMatchesFirstOrder()
{
  using namespace Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions;
  Tensor33d C = get< 0 >( computationParameters() );
  checkPlainMatchesFirstOrder(
    __PRETTY_FUNCTION__,
    []( const auto& C, double c10, double c01 ) { return MooneyRivlinPotential( C, c10, c01 ); },
    []( const auto& C, double c10, double c01 ) { return FirstOrder::MooneyRivlinPotential( C, c10, c01 ); },
    C,
    400.0,
    100.0 );
}

auto testVolumetricPenaltyPlainMatchesFirstOrder()
{
  using namespace Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions;
  Tensor33d C = get< 0 >( computationParameters() );
  checkPlainMatchesFirstOrder(
    __PRETTY_FUNCTION__,
    []( const auto& C, double kappa ) { return VolumetricPenaltyPotential( C, kappa ); },
    []( const auto& C, double kappa ) { return FirstOrder::VolumetricPenaltyPotential( C, kappa ); },
    C,
    10000.0 );
}

// Generic checks for a volumetric strain energy function U(J) acc. Doll & Schweizerhof (2000),
// https://doi.org/10.1115/1.321146:
//  - plain vs FirstOrder (complex-step derivative of the plain overload),
//  - closed-form value at the test state C of computationParameters() (J = 2),
//  - conditions (5)_2 and (6): U = 0 and dU/dC = 0 at C = I,
//  - condition (8)_2: d^2U/dJ^2 = kappa at J = 1 (central difference of the analytic dU/dJ along C = J^{2/3} I).
template < typename PlainFn, typename FirstOrderFn, typename... Args >
void checkVolumetricPotential( const char*  testName,
                               PlainFn      plainFn,
                               FirstOrderFn firstOrderFn,
                               const double expectedAtJ2,
                               Args... args )
{
  const Tensor33d C = get< 0 >( computationParameters() );
  checkPlainMatchesFirstOrder( testName, plainFn, firstOrderFn, C, args... );

  throwExceptionOnFailure( checkIfEqual( plainFn( C, args... ), expectedAtJ2, 1e-13 ),
                           Marmot::Utilities::MakeString()
                             << testName << " energy density does not match the closed form at J = 2." );

  const Tensor33d I = Marmot::Math::Tensors::fastorTensorFromDoubleTensor< double >( Spatial3D::I );
  double          U0;
  Tensor33d       dU0_dC;
  std::tie( U0, dU0_dC ) = firstOrderFn( I, args... );
  throwExceptionOnFailure( std::abs( U0 ) < 1e-12 && checkIfEqual( dU0_dC, Tensor33d( 0.0 ), 1e-12 ),
                           Marmot::Utilities::MakeString()
                             << testName << " energy or its derivative is non-zero at C = I." );

  // volumetric stress p = dU/dJ from the analytic gradient: along C = J^{2/3} I, dU/dC = p J/2 C^{-T} = p J^{1/3}/2 I
  const double kappa = get< 0 >( std::make_tuple( args... ) );
  const double h     = 1e-4;
  auto         pOfJ  = [&]( double J ) {
    const auto [U, dU_dC] = firstOrderFn( Tensor33d( std::pow( J, 2. / 3. ) * I ), args... );
    return 2. * dU_dC( 0, 0 ) / std::cbrt( J );
  };
  const double d2U_dJ = ( pOfJ( 1. + h ) - pOfJ( 1. - h ) ) / ( 2. * h );
  throwExceptionOnFailure( checkIfEqual( d2U_dJ / kappa, 1.0, 1e-6 ),
                           Marmot::Utilities::MakeString()
                             << testName << " d2U/dJ2 at J = 1 does not equal the bulk modulus." );
}

auto testVolumetricPotentialsDollSchweizerhof()
{
  using namespace Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions;
  const double k   = 10000.0;
  const double ln2 = std::log( 2.0 );

  // U_1 .. U_5, Table 1 (U_2 = VolumetricPenaltyPotential)
  checkVolumetricPotential(
    "SussmanBatheVolumetricPotential (U1)",
    []( const auto& C, double kappa ) { return SussmanBatheVolumetricPotential( C, kappa ); },
    []( const auto& C, double kappa ) { return FirstOrder::SussmanBatheVolumetricPotential( C, kappa ); },
    k / 2.,
    k );
  checkVolumetricPotential(
    "VolumetricPenaltyPotential (U2)",
    []( const auto& C, double kappa ) { return VolumetricPenaltyPotential( C, kappa ); },
    []( const auto& C, double kappa ) { return FirstOrder::VolumetricPenaltyPotential( C, kappa ); },
    k / 2. * ln2 * ln2,
    k );
  checkVolumetricPotential(
    "SimoTaylorVolumetricPotential (U3)",
    []( const auto& C, double kappa ) { return SimoTaylorVolumetricPotential( C, kappa ); },
    []( const auto& C, double kappa ) { return FirstOrder::SimoTaylorVolumetricPotential( C, kappa ); },
    k / 4. * ( 1. + ln2 * ln2 ),
    k );
  for ( const double theta : { -2.0, -5.0, 9.0 } )
    checkVolumetricPotential(
      "OgdenVolumetricPotential (U4)",
      []( const auto& C, double kappa, double th ) { return OgdenVolumetricPotential( C, kappa, th ); },
      []( const auto& C, double kappa, double th ) { return FirstOrder::OgdenVolumetricPotential( C, kappa, th ); },
      k / ( theta * theta ) * ( theta * ln2 + std::pow( 2., -theta ) - 1. ),
      k,
      theta );
  checkVolumetricPotential(
    "LiuHofstetterMangVolumetricPotential (U5)",
    []( const auto& C, double kappa ) { return LiuHofstetterMangVolumetricPotential( C, kappa ); },
    []( const auto& C, double kappa ) { return FirstOrder::LiuHofstetterMangVolumetricPotential( C, kappa ); },
    k * ( 2. * ln2 - 1. ),
    k );

  // U_6, Eq. (13), with the parameter sets U_6a .. U_6f of the paper
  for ( const auto& [a, b] : std::vector< std::pair< double, double > >{ { 2., 2. },
                                                                         { 0.001, 3.79 },
                                                                         { 1., 1.001 },
                                                                         { 0.001, 10. },
                                                                         { 2.3, 1.4 },
                                                                         { 0.45, 1.05 } } )
    checkVolumetricPotential(
      "DollSchweizerhofVolumetricPotential (U6)",
      []( const auto& C, double kappa, double alpha, double beta ) {
        return DollSchweizerhofVolumetricPotential( C, kappa, alpha, beta );
      },
      []( const auto& C, double kappa, double alpha, double beta ) {
        return FirstOrder::DollSchweizerhofVolumetricPotential( C, kappa, alpha, beta );
      },
      k * ( ( std::pow( 2., a + 1. ) / ( a + 1. ) + std::pow( 2., 1. - b ) / ( b - 1. ) ) / ( a + b ) -
            1. / ( ( a + 1. ) * ( b - 1. ) ) ),
      k,
      a,
      b );

  // U_6gen, Eq. (15)
  {
    const std::vector< double > alphas   = { 0.5, 2.0 };
    const std::vector< double > betas    = { 1.5, 3.0 };
    double                      expected = 0.;
    for ( std::size_t i = 0; i < alphas.size(); i++ )
      expected += ( std::pow( 2., alphas[i] + 1. ) - 1. ) / ( alphas[i] + 1. ) +
                  ( std::pow( 2., 1. - betas[i] ) - 1. ) / ( betas[i] - 1. );
    expected *= k / 7.0;
    checkVolumetricPotential(
      "GeneralizedDollSchweizerhofVolumetricPotential (U6gen)",
      [&]( const auto& C, double kappa ) {
        return GeneralizedDollSchweizerhofVolumetricPotential( C, kappa, alphas, betas );
      },
      [&]( const auto& C, double kappa ) {
        return FirstOrder::GeneralizedDollSchweizerhofVolumetricPotential( C, kappa, alphas, betas );
      },
      expected,
      k );
  }

  // U_7, U_8, Eq. (16)
  checkVolumetricPotential(
    "DollSchweizerhofExponentialVolumetricPotential (U7)",
    []( const auto& C, double kappa ) { return DollSchweizerhofExponentialVolumetricPotential( C, kappa ); },
    []( const auto& C, double kappa ) {
      return FirstOrder::DollSchweizerhofExponentialVolumetricPotential( C, kappa );
    },
    k / 2. * ( std::exp( 1. ) - ln2 - 1. ),
    k );
  checkVolumetricPotential(
    "DollSchweizerhofLogarithmicVolumetricPotential (U8)",
    []( const auto& C, double kappa ) { return DollSchweizerhofLogarithmicVolumetricPotential( C, kappa ); },
    []( const auto& C, double kappa ) {
      return FirstOrder::DollSchweizerhofLogarithmicVolumetricPotential( C, kappa );
    },
    k / 2. * ln2,
    k );
}

// Relations between the volumetric potentials stated in Doll & Schweizerhof (2000).
auto testVolumetricPotentialRelationsDollSchweizerhof()
{
  using namespace Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions;
  const double    k = 10000.0;
  const Tensor33d C = get< 0 >( computationParameters() );
  const Tensor33d CInv( Fastor::inverse( C ) );

  // U_3 = (U_1 + U_2) / 2
  throwExceptionOnFailure( checkIfEqual( SimoTaylorVolumetricPotential( C, k ),
                                         0.5 * ( SussmanBatheVolumetricPotential( C, k ) +
                                                 VolumetricPenaltyPotential( C, k ) ),
                                         1e-14 ),
                           Marmot::Utilities::MakeString() << __PRETTY_FUNCTION__ << " U3 != (U1 + U2)/2." );

  // U_4 with theta = -2 equals kappa/4 (J^2 - 1 - 2 ln J), the volumetric part of standardNeoHooke
  const double J = std::sqrt( Fastor::determinant( C ) );
  throwExceptionOnFailure( checkIfEqual( OgdenVolumetricPotential( C, k, -2.0 ),
                                         k / 4. * ( J * J - 1. - 2. * std::log( J ) ),
                                         1e-14 ),
                           Marmot::Utilities::MakeString()
                             << __PRETTY_FUNCTION__ << " U4(theta=-2) != kappa/4 (J^2 - 1 - 2 ln J)." );

  // beta = alpha + 2: U_6(J) = U_6(1/J)
  throwExceptionOnFailure( checkIfEqual( DollSchweizerhofVolumetricPotential( C, k, 1.5, 3.5 ),
                                         DollSchweizerhofVolumetricPotential( CInv, k, 1.5, 3.5 ),
                                         1e-13 ),
                           Marmot::Utilities::MakeString()
                             << __PRETTY_FUNCTION__ << " U6(J) != U6(1/J) for beta = alpha + 2." );

  // beta = alpha: p_6(J) = -p_6(1/J), with p = dU/dJ = 2/J dU/dC : C / 3 ... evaluated via dU/dC = p J/2 C^{-T}
  {
    double    U, UInv;
    Tensor33d dU_dC, dUInv_dC;
    std::tie( U, dU_dC )       = FirstOrder::DollSchweizerhofVolumetricPotential( C, k, 2.0, 2.0 );
    std::tie( UInv, dUInv_dC ) = FirstOrder::DollSchweizerhofVolumetricPotential( CInv, k, 2.0, 2.0 );
    // p = 2/J (dU/dC : C) / 3, since C^{-T} : C = 3
    const double p    = 2. / J * Fastor::inner( dU_dC, C ) / 3.;
    const double pInv = 2. * J * Fastor::inner( dUInv_dC, CInv ) / 3.;
    throwExceptionOnFailure( checkIfEqual( p, -pInv, 1e-13 ),
                             Marmot::Utilities::MakeString()
                               << __PRETTY_FUNCTION__ << " p6(J) != -p6(1/J) for beta = alpha." );
  }

  // U_6gen with a duplicated parameter pair reduces to U_6
  {
    const std::vector< double > alphas = { 0.7, 0.7 };
    const std::vector< double > betas  = { 2.5, 2.5 };
    throwExceptionOnFailure( checkIfEqual( GeneralizedDollSchweizerhofVolumetricPotential( C, k, alphas, betas ),
                                           DollSchweizerhofVolumetricPotential( C, k, 0.7, 2.5 ),
                                           1e-14 ),
                             Marmot::Utilities::MakeString()
                               << __PRETTY_FUNCTION__ << " U6gen with duplicated pairs != U6." );
  }

  // mismatched parameter arrays are rejected
  {
    const std::vector< double > alphas = { 0.7, 1.0 };
    const std::vector< double > betas  = { 2.5 };
    bool                        thrown = false;
    try {
      GeneralizedDollSchweizerhofVolumetricPotential( C, k, alphas, betas );
    }
    catch ( const std::invalid_argument& ) {
      thrown = true;
    }
    throwExceptionOnFailure( thrown,
                             Marmot::Utilities::MakeString()
                               << __PRETTY_FUNCTION__ << " U6gen accepted mismatched alphas/betas." );
  }
}

// Yeoh and Mooney-Rivlin must reduce exactly to NeoHooke's isochoric part when
// the extra shape parameters vanish and C10=mu/2, since all three now share
// the same isochoric invariant Ibar1 (no linear-shift terms to cause drift).
auto testYeohMooneyRivlinReduceToNeoHooke()
{
  using namespace Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions;
  Tensor33d    C  = get< 0 >( computationParameters() );
  const double mu = 750.0;

  const double psiNeoHooke = NeoHookePotential( C, mu );
  const double psiYeoh     = YeohPotential( C, mu / 2., 0.0, 0.0 );
  const double psiMR       = MooneyRivlinPotential( C, mu / 2., 0.0 );

  throwExceptionOnFailure( checkIfEqual( psiYeoh, psiNeoHooke, 1e-13 ),
                           Marmot::Utilities::MakeString()
                             << __PRETTY_FUNCTION__ << " Yeoh does not reduce to NeoHooke." );
  throwExceptionOnFailure( checkIfEqual( psiMR, psiNeoHooke, 1e-13 ),
                           Marmot::Utilities::MakeString()
                             << __PRETTY_FUNCTION__ << " Mooney-Rivlin does not reduce to NeoHooke." );
}

auto testStandardNeoHookeEnergyMatchesClosedForm()
{
  using namespace Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions::ThirdOrder;
  std::tuple< Tensor33d, Tensor33d, double, double, double, double > params = computationParameters();
  Tensor33d                                                          C      = get< 0 >( params );
  const double                                                       K      = get< 2 >( params );
  const double                                                       G      = get< 3 >( params );

  const double lambda   = K - 2.0 / 3.0 * G;
  const double trC      = Fastor::trace( C );
  const double detC     = Fastor::determinant( C );
  const double expected = G / 2. * ( trC - 3.0 - log( detC ) ) + lambda / 4. * ( detC - 1.0 - log( detC ) );

  auto [psi, dPsi_dC, d2Psi_dC2, d3Psi_dC3] = standardNeoHooke( C, K, G );

  throwExceptionOnFailure( checkIfEqual( psi, expected, 1e-10 ),
                           Marmot::Utilities::MakeString()
                             << __PRETTY_FUNCTION__ << " energy density does not match the closed form" );
}

auto testStandardNeoHookeFirstDerivativeMatchesNumericalDifferentiation()
{
  using namespace Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions::ThirdOrder;
  std::tuple< Tensor33d, Tensor33d, double, double, double, double > params = computationParameters();
  Tensor33d                                                          C      = get< 0 >( params );
  const double                                                       K      = get< 2 >( params );
  const double                                                       G      = get< 3 >( params );

  auto [psi0, dPsi_dC, d2Psi_dC2, d3Psi_dC3] = standardNeoHooke( C, K, G );

  const double h = 1e-6;
  for ( int i = 0; i < 3; i++ )
    for ( int j = 0; j < 3; j++ ) {
      Tensor33d Cp = C;
      Tensor33d Cm = C;
      Cp( i, j ) += h;
      Cm( i, j ) -= h;
      const double psiPlus  = std::get< 0 >( standardNeoHooke( Cp, K, G ) );
      const double psiMinus = std::get< 0 >( standardNeoHooke( Cm, K, G ) );
      const double numDeriv = ( psiPlus - psiMinus ) / ( 2. * h );
      throwExceptionOnFailure( checkIfEqual( dPsi_dC( i, j ), numDeriv, 1e-6 ),
                               Marmot::Utilities::MakeString()
                                 << __PRETTY_FUNCTION__ << " dPsi_dC(" << i << "," << j
                                 << ") does not match central-difference numerical differentiation" );
    }
}

// The 6 independent symmetric "unit" directions of a symmetric 3x3 tensor (00, 11, 22, and the
// 3 off-diagonal directions with both (a,b) and (b,a) set to 1). d2Psi_dC2/d3Psi_dC3 are only
// meaningful when contracted against directions like these -- C is always the (symmetric) right
// Cauchy-Green tensor, and the returned tensors internally symmetrize over each differentiated
// index pair to match that (e.g. d2Psi_dC2(i,j,k,l) == d2Psi_dC2(i,j,l,k) by construction), so
// isolating a single (k,l) component via an asymmetric single-entry perturbation does not recover
// it. Verifying directional bilinear/trilinear forms instead sidesteps that convention entirely.
std::vector< Tensor33d > symmetricUnitDirections()
{
  std::vector< Tensor33d > directions;
  for ( int a = 0; a < 3; a++ )
    for ( int b = a; b < 3; b++ ) {
      Tensor33d dC( 0.0 );
      dC( a, b ) = 1.0;
      dC( b, a ) = 1.0;
      directions.push_back( dC );
    }
  return directions;
}

auto testStandardNeoHookeSecondDerivativeMatchesNumericalDifferentiation()
{
  using namespace Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions::ThirdOrder;
  std::tuple< Tensor33d, Tensor33d, double, double, double, double > params = computationParameters();
  Tensor33d                                                          C      = get< 0 >( params );
  const double                                                       K      = get< 2 >( params );
  const double                                                       G      = get< 3 >( params );

  auto [psi0, dPsi_dC0, d2Psi_dC2, d3Psi_dC3] = standardNeoHooke( C, K, G );

  const double h = 1e-5;
  for ( const Tensor33d& dC : symmetricUnitDirections() ) {
    Tensor33d Cp = C + h * dC;
    Tensor33d Cm = C - h * dC;

    const Tensor33d dPsi_dC_plus  = std::get< 1 >( standardNeoHooke( Cp, K, G ) );
    const Tensor33d dPsi_dC_minus = std::get< 1 >( standardNeoHooke( Cm, K, G ) );

    for ( int i = 0; i < 3; i++ )
      for ( int j = 0; j < 3; j++ ) {
        const double numDirectional = ( dPsi_dC_plus( i, j ) - dPsi_dC_minus( i, j ) ) / ( 2. * h );

        double analyticalDirectional = 0.0;
        for ( int k = 0; k < 3; k++ )
          for ( int l = 0; l < 3; l++ )
            analyticalDirectional += d2Psi_dC2( i, j, k, l ) * dC( k, l );

        throwExceptionOnFailure( checkIfEqual( analyticalDirectional, numDirectional, 1e-4 ),
                                 Marmot::Utilities::MakeString()
                                   << __PRETTY_FUNCTION__ << " d2Psi_dC2 contracted with a symmetric "
                                   << "direction does not match numerical differentiation at (" << i << "," << j
                                   << ")" );
      }
  }
}

auto testStandardNeoHookeThirdDerivativeMatchesNumericalDifferentiation()
{
  using namespace Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions::ThirdOrder;
  std::tuple< Tensor33d, Tensor33d, double, double, double, double > params = computationParameters();
  Tensor33d                                                          C      = get< 0 >( params );
  const double                                                       K      = get< 2 >( params );
  const double                                                       G      = get< 3 >( params );

  auto [psi0, dPsi_dC0, d2Psi_dC2_0, d3Psi_dC3] = standardNeoHooke( C, K, G );

  const std::vector< Tensor33d > directions = symmetricUnitDirections();

  const double h = 1e-4;
  for ( const Tensor33d& dCOuter : directions ) {
    Tensor33d Cp = C + h * dCOuter;
    Tensor33d Cm = C - h * dCOuter;

    const Tensor3333d d2Psi_dC2_plus  = std::get< 2 >( standardNeoHooke( Cp, K, G ) );
    const Tensor3333d d2Psi_dC2_minus = std::get< 2 >( standardNeoHooke( Cm, K, G ) );

    for ( const Tensor33d& dCInner : directions )
      for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ ) {
          double d2Plus = 0.0, d2Minus = 0.0;
          for ( int k = 0; k < 3; k++ )
            for ( int l = 0; l < 3; l++ ) {
              d2Plus += d2Psi_dC2_plus( i, j, k, l ) * dCInner( k, l );
              d2Minus += d2Psi_dC2_minus( i, j, k, l ) * dCInner( k, l );
            }
          const double numDirectional = ( d2Plus - d2Minus ) / ( 2. * h );

          double analyticalDirectional = 0.0;
          for ( int k = 0; k < 3; k++ )
            for ( int l = 0; l < 3; l++ )
              for ( int m = 0; m < 3; m++ )
                for ( int n = 0; n < 3; n++ )
                  analyticalDirectional += d3Psi_dC3( i, j, k, l, m, n ) * dCInner( k, l ) * dCOuter( m, n );

          throwExceptionOnFailure( checkIfEqual( analyticalDirectional, numDirectional, 1e-2 ),
                                   Marmot::Utilities::MakeString()
                                     << __PRETTY_FUNCTION__ << " d3Psi_dC3 contracted with two symmetric "
                                     << "directions does not match numerical differentiation at (" << i << "," << j
                                     << ")" );
        }
  }
}

int main()
{
  auto tests = std::vector<
    std::function< void() > >{ testPenceGouPotentialA,
                               testPenceGouPotentialB,
                               testPenceGouPotentialC,
                               testFirstOrderDerivedB,
                               testSecondOrderDerivedB,
                               testArrudaBoycePlainMatchesFirstOrder,
                               testArrudaBoyceZeroEnergyAtReferenceForNearOneLockingStretch,
                               testNeoHookePlainMatchesFirstOrder,
                               testYeohPlainMatchesFirstOrder,
                               testMooneyRivlinPlainMatchesFirstOrder,
                               testVolumetricPenaltyPlainMatchesFirstOrder,
                               testVolumetricPotentialsDollSchweizerhof,
                               testVolumetricPotentialRelationsDollSchweizerhof,
                               testYeohMooneyRivlinReduceToNeoHooke,
                               testStandardNeoHookeEnergyMatchesClosedForm,
                               testStandardNeoHookeFirstDerivativeMatchesNumericalDifferentiation,
                               testStandardNeoHookeSecondDerivativeMatchesNumericalDifferentiation,
                               testStandardNeoHookeThirdDerivativeMatchesNumericalDifferentiation };

  executeTestsAndCollectExceptions( tests );
  return 0;
}
