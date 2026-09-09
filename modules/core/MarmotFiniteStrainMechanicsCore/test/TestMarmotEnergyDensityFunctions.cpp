#include "Fastor/Fastor.h"
#include "Marmot/MarmotEnergyDensityFunctions.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotTesting.h"
#include <complex>

using namespace Marmot::Testing;
using namespace Marmot::FastorStandardTensors;

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
  using namespace Marmot::ContinuumMechanics::EnergyDensityFunctions;
  // get computation parameters
  std::tuple< Tensor33d, Tensor33d, double, double, double, double > params = computationParameters();
  Tensor33d                                                          C      = get< 0 >( params );
  const double                                                       K      = get< 2 >( params );
  const double                                                       G      = get< 3 >( params );
  // compute expected and resulting value and compare
  double expected = 5098.519486106722;
  double res      = PenceGouPotentialA( C, K, G );
  throwExceptionOnFailure( checkIfEqual( res, expected, 1e-15 ),
                           MakeString() << __PRETTY_FUNCTION__ << " energy density A function failed." );
}

auto testPenceGouPotentialB()
{
  using namespace Marmot::ContinuumMechanics::EnergyDensityFunctions;
  // get computation parameters
  std::tuple< Tensor33d, Tensor33d, double, double, double, double > params = computationParameters();
  Tensor33d                                                          C      = get< 0 >( params );
  const double                                                       K      = get< 2 >( params );
  const double                                                       G      = get< 3 >( params );
  // compute expected and resulting value and compare
  double expected = 2966.146377987021;
  double res      = PenceGouPotentialB( C, K, G );
  throwExceptionOnFailure( checkIfEqual( res, expected, 1e-15 ),
                           MakeString() << __PRETTY_FUNCTION__ << " energy density B function failed." );
}

auto testPenceGouPotentialC()
{
  using namespace Marmot::ContinuumMechanics::EnergyDensityFunctions;
  // get computation parameters
  std::tuple< Tensor33d, Tensor33d, double, double, double, double > params = computationParameters();
  Tensor33d                                                          C      = get< 0 >( params );
  const double                                                       K      = get< 2 >( params );
  const double                                                       G      = get< 3 >( params );
  // compute expected and resulting value and compare
  double expected = 1018.0232353221228;
  double res      = PenceGouPotentialC( C, K, G );
  throwExceptionOnFailure( checkIfEqual( res, expected, 1e-15 ),
                           MakeString() << __PRETTY_FUNCTION__ << " energy density C function failed." );
}

auto testFirstOrderDerivedB()
{
  using namespace Marmot::ContinuumMechanics::EnergyDensityFunctions::FirstOrderDerived;
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
                           MakeString() << __PRETTY_FUNCTION__ << " energy density B function failed." );
  // constant parameters
  const double dW_dI1 = G / 2. * pow( J, -2. / 3 );
  const double dW_dJ  = K / 4. * ( J - pow( J, -3. ) ) - G * I1 / 3. * pow( J, -5. / 3 );
  double       expectedS;
  // compare stress
  for ( int i = 0; i < 3; i++ )
    for ( int j = 0; j < 3; j++ ) {
      expectedS = dW_dI1 * ( i == j ? 1 : 0 ) + dW_dJ * J / 2. * invC( i, j );
      throwExceptionOnFailure( checkIfEqual( resS( i, j ), expectedS, 1e-12 ),
                               MakeString() << __PRETTY_FUNCTION__ << " stress computation B failed for position (" << i
                                            << ", " << j << ")" );
    }
}

auto testSecondOrderDerivedB()
{
  using namespace Marmot::ContinuumMechanics::EnergyDensityFunctions::SecondOrderDerived;
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
                           MakeString() << __PRETTY_FUNCTION__ << " energy density B function failed." );
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
                               MakeString() << __PRETTY_FUNCTION__ << " stress computation B failed for position (" << i
                                            << ", " << j << ")" );
      for ( int k = 0; k < 3; k++ )
        for ( int l = 0; l < 3; l++ ) {
          expectedCSE = ( J * J ) / 4 * d2W_dJ2 * invC( j, i ) * invC( l, k ) +
                        dW_dJ * ( J / 4. * invC( j, i ) * invC( l, k ) - J / 2. * invC( j, k ) * invC( l, i ) ) +
                        d2W_dI1dJ * J / 2. * ( invC( j, i ) * ( k == l ? 1 : 0 ) + ( i == j ? 1 : 0 ) * invC( l, k ) );
          throwExceptionOnFailure( checkIfEqual( resCSE( i, j, k, l ), expectedCSE, 1e-12 ),
                                   MakeString() << __PRETTY_FUNCTION__ << " tangent computation B failed for position ("
                                                << i << ", " << j << ", " << k << ", " << l << ")" );
        }
    }
}

// Cross-check that the plain (energy-only) and FirstOrderDerived (energy + analytic
// dPsi/dC) overloads of ArrudaBoyce8ChainPotential cannot drift out of sync: both
// funnel through the same internal `detail::arrudaBoyce8ChainEnergyAndDerivative`
// helper, but this test verifies that end-to-end rather than assuming it. The
// gradient is cross-checked against a complex-step derivative of the plain overload
// itself (not against a hand-derived formula), so this specifically guards against
// the two overloads' shared helper silently diverging from what the plain overload
// actually computes.
auto testArrudaBoycePlainMatchesFirstOrderDerived()
{
  using namespace Marmot::ContinuumMechanics::EnergyDensityFunctions;

  Tensor33d    C       = get< 0 >( computationParameters() );
  const double mu      = 500.0;
  const double lambdaL = 2.5;

  const double psiPlain = ArrudaBoyce8ChainPotential< double >( C, mu, lambdaL );

  double    psiFOD;
  Tensor33d dPsiFOD_dC;
  std::tie( psiFOD, dPsiFOD_dC ) = FirstOrderDerived::ArrudaBoyce8ChainPotential< double >( C, mu, lambdaL );

  throwExceptionOnFailure( checkIfEqual( psiPlain, psiFOD, 1e-14 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " plain and FirstOrderDerived energy densities disagree." );

  // complex-step derivative of the plain overload, entry by entry (C treated as a
  // general, not-necessarily-symmetric matrix argument, matching how dPsi/dC is
  // defined and used everywhere else in this codebase)
  const double h = 1e-20;
  Tensor33d    dPsiCSDA( 0.0 );
  for ( int i = 0; i < 3; i++ )
    for ( int j = 0; j < 3; j++ ) {
      Tensor33t< std::complex< double > > Cc = Marmot::fastorTensorFromDoubleTensor< std::complex< double > >( C );
      Cc( i, j ) += std::complex< double >( 0.0, h );
      std::complex< double > psiC = ArrudaBoyce8ChainPotential< std::complex< double > >( Cc, mu, lambdaL );
      dPsiCSDA( i, j )             = psiC.imag() / h;
    }

  throwExceptionOnFailure( checkIfEqual( dPsiFOD_dC, dPsiCSDA, 1e-9 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " FirstOrderDerived gradient disagrees with a complex-step derivative "
                                           "of the plain overload." );
}

// Generic plain-vs-FirstOrderDerived cross-check, same methodology as
// testArrudaBoycePlainMatchesFirstOrderDerived: the FirstOrderDerived gradient
// is checked against a complex-step derivative of the plain overload itself,
// guarding against the two overloads' formulas silently drifting apart.
template < typename PlainFn, typename FirstOrderFn, typename... Args >
void checkPlainMatchesFirstOrderDerived( const char*   testName,
                                         PlainFn       plainFn,
                                         FirstOrderFn  firstOrderFn,
                                         const Tensor33d& C,
                                         Args... args )
{
  const double psiPlain = plainFn( C, args... );

  double    psiFOD;
  Tensor33d dPsiFOD_dC;
  std::tie( psiFOD, dPsiFOD_dC ) = firstOrderFn( C, args... );

  throwExceptionOnFailure( checkIfEqual( psiPlain, psiFOD, 1e-14 ),
                           MakeString() << testName << " plain and FirstOrderDerived energy densities disagree." );

  const double h = 1e-20;
  Tensor33d    dPsiCSDA( 0.0 );
  for ( int i = 0; i < 3; i++ )
    for ( int j = 0; j < 3; j++ ) {
      Tensor33t< std::complex< double > > Cc = Marmot::fastorTensorFromDoubleTensor< std::complex< double > >( C );
      Cc( i, j ) += std::complex< double >( 0.0, h );
      std::complex< double > psiC = plainFn( Cc, args... );
      dPsiCSDA( i, j )             = psiC.imag() / h;
    }

  throwExceptionOnFailure( checkIfEqual( dPsiFOD_dC, dPsiCSDA, 1e-9 ),
                           MakeString() << testName
                                        << " FirstOrderDerived gradient disagrees with a complex-step derivative "
                                           "of the plain overload." );
}

auto testNeoHookePlainMatchesFirstOrderDerived()
{
  using namespace Marmot::ContinuumMechanics::EnergyDensityFunctions;
  Tensor33d C = get< 0 >( computationParameters() );
  checkPlainMatchesFirstOrderDerived(
    __PRETTY_FUNCTION__,
    []( const auto& C, double mu ) { return NeoHookePotential( C, mu ); },
    []( const auto& C, double mu ) { return FirstOrderDerived::NeoHookePotential( C, mu ); },
    C,
    750.0 );
}

auto testYeohPlainMatchesFirstOrderDerived()
{
  using namespace Marmot::ContinuumMechanics::EnergyDensityFunctions;
  Tensor33d C = get< 0 >( computationParameters() );
  checkPlainMatchesFirstOrderDerived(
    __PRETTY_FUNCTION__,
    []( const auto& C, double c10, double c20, double c30 ) { return YeohPotential( C, c10, c20, c30 ); },
    []( const auto& C, double c10, double c20, double c30 ) { return FirstOrderDerived::YeohPotential( C, c10, c20, c30 ); },
    C,
    500.0,
    50.0,
    5.0 );
}

auto testMooneyRivlinPlainMatchesFirstOrderDerived()
{
  using namespace Marmot::ContinuumMechanics::EnergyDensityFunctions;
  Tensor33d C = get< 0 >( computationParameters() );
  checkPlainMatchesFirstOrderDerived(
    __PRETTY_FUNCTION__,
    []( const auto& C, double c10, double c01 ) { return MooneyRivlinPotential( C, c10, c01 ); },
    []( const auto& C, double c10, double c01 ) { return FirstOrderDerived::MooneyRivlinPotential( C, c10, c01 ); },
    C,
    400.0,
    100.0 );
}

auto testVolumetricPenaltyPlainMatchesFirstOrderDerived()
{
  using namespace Marmot::ContinuumMechanics::EnergyDensityFunctions;
  Tensor33d C = get< 0 >( computationParameters() );
  checkPlainMatchesFirstOrderDerived(
    __PRETTY_FUNCTION__,
    []( const auto& C, double kappa ) { return VolumetricPenaltyPotential( C, kappa ); },
    []( const auto& C, double kappa ) { return FirstOrderDerived::VolumetricPenaltyPotential( C, kappa ); },
    C,
    10000.0 );
}

// Yeoh and Mooney-Rivlin must reduce exactly to NeoHooke's isochoric part when
// the extra shape parameters vanish and C10=mu/2, since all three now share
// the same isochoric invariant Ibar1 (no linear-shift terms to cause drift).
auto testYeohMooneyRivlinReduceToNeoHooke()
{
  using namespace Marmot::ContinuumMechanics::EnergyDensityFunctions;
  Tensor33d    C  = get< 0 >( computationParameters() );
  const double mu = 750.0;

  const double psiNeoHooke = NeoHookePotential( C, mu );
  const double psiYeoh     = YeohPotential( C, mu / 2., 0.0, 0.0 );
  const double psiMR       = MooneyRivlinPotential( C, mu / 2., 0.0 );

  throwExceptionOnFailure( checkIfEqual( psiYeoh, psiNeoHooke, 1e-13 ),
                           MakeString() << __PRETTY_FUNCTION__ << " Yeoh does not reduce to NeoHooke." );
  throwExceptionOnFailure( checkIfEqual( psiMR, psiNeoHooke, 1e-13 ),
                           MakeString() << __PRETTY_FUNCTION__ << " Mooney-Rivlin does not reduce to NeoHooke." );
}

int main()
{
  auto tests = std::vector< std::function< void() > >{ testPenceGouPotentialA,
                                                       testPenceGouPotentialB,
                                                       testPenceGouPotentialC,
                                                       testFirstOrderDerivedB,
                                                       testSecondOrderDerivedB,
                                                       testArrudaBoycePlainMatchesFirstOrderDerived,
                                                       testNeoHookePlainMatchesFirstOrderDerived,
                                                       testYeohPlainMatchesFirstOrderDerived,
                                                       testMooneyRivlinPlainMatchesFirstOrderDerived,
                                                       testVolumetricPenaltyPlainMatchesFirstOrderDerived,
                                                       testYeohMooneyRivlinReduceToNeoHooke };

  executeTestsAndCollectExceptions( tests );
  return 0;
}
