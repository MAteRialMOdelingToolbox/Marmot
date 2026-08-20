#include "Fastor/Fastor.h"
#include "Marmot/MarmotEnergyDensityFunctions.h"
#include "Marmot/MarmotTesting.h"
#include <array>

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
  using namespace Marmot::ContinuumMechanics::EnergyDensityFunctions::Compressible;
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
  using namespace Marmot::ContinuumMechanics::EnergyDensityFunctions::Compressible;
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
  using namespace Marmot::ContinuumMechanics::EnergyDensityFunctions::Compressible;
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
  using namespace Marmot::ContinuumMechanics::EnergyDensityFunctions::Compressible::FirstOrderDerived;
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
  using namespace Marmot::ContinuumMechanics::EnergyDensityFunctions::Compressible::SecondOrderDerived;
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

namespace {

  // Six independent symmetric directions spanning the space of symmetric perturbations of C. C is
  // always symmetric in practice (C=F^T F), so a rigorous finite-difference check of a
  // dPsi/dC-type derivative must probe it only along symmetric directions - perturbing a single
  // off-diagonal component C(i,j) alone (leaving C(j,i) untouched) produces an unphysical
  // asymmetric C that the analytic (symmetric-domain) derivative formulas are not meant to match.
  std::array< Tensor33d, 6 > symmetricBasisDirections()
  {
    std::array< Tensor33d, 6 > D;
    for ( auto& d : D )
      d = Tensor33d( 0.0 );
    D[0]( 0, 0 ) = 1.;
    D[1]( 1, 1 ) = 1.;
    D[2]( 2, 2 ) = 1.;
    D[3]( 0, 1 ) = D[3]( 1, 0 ) = 1.;
    D[4]( 0, 2 ) = D[4]( 2, 0 ) = 1.;
    D[5]( 1, 2 ) = D[5]( 2, 1 ) = 1.;
    return D;
  }

  // Verifies {psi, dPsi_dC, d2Psi_dCdC} = potentialFn(C) against central-difference directional
  // derivatives along every pair of symmetric basis directions.
  template < typename PotentialFn >
  void checkSecondOrderDerivedAgainstFD( const std::string& testName, const Tensor33d& C, PotentialFn potentialFn )
  {
    const auto [psi, dPsi_dC, d2Psi_dCdC] = potentialFn( C );

    static const std::array< Tensor33d, 6 > D = symmetricBasisDirections();
    const double                            h = 1e-6;

    for ( int a = 0; a < 6; a++ ) {
      const Tensor33d Cp = C + h * D[a];
      const Tensor33d Cm = C - h * D[a];

      const auto [psiP, dPsi_dCp, d2Psi_dCdCp] = potentialFn( Cp );
      const auto [psiM, dPsi_dCm, d2Psi_dCdCm] = potentialFn( Cm );

      const double directionalFD       = ( psiP - psiM ) / ( 2 * h );
      double       directionalAnalytic = 0.;
      for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ )
          directionalAnalytic += dPsi_dC( i, j ) * D[a]( i, j );

      throwExceptionOnFailure( checkIfEqual( directionalAnalytic, directionalFD, 1e-6 ),
                               MakeString() << testName << " first derivative mismatch along direction " << a );

      for ( int b = 0; b < 6; b++ ) {
        double dirP = 0., dirM = 0.;
        for ( int k = 0; k < 3; k++ )
          for ( int l = 0; l < 3; l++ ) {
            dirP += dPsi_dCp( k, l ) * D[b]( k, l );
            dirM += dPsi_dCm( k, l ) * D[b]( k, l );
          }
        const double d2DirectionalFD = ( dirP - dirM ) / ( 2 * h );

        double d2DirectionalAnalytic = 0.;
        for ( int i = 0; i < 3; i++ )
          for ( int j = 0; j < 3; j++ )
            for ( int k = 0; k < 3; k++ )
              for ( int l = 0; l < 3; l++ )
                d2DirectionalAnalytic += d2Psi_dCdC( i, j, k, l ) * D[a]( i, j ) * D[b]( k, l );

        throwExceptionOnFailure( checkIfEqual( d2DirectionalAnalytic, d2DirectionalFD, 1e-3 ),
                                 MakeString() << testName << " second derivative mismatch along directions (" << a
                                              << ", " << b << ")" );
      }
    }
  }

} // namespace

// Verifies the analytic first and second derivative of the (purely isochoric)
// Incompressible::NeoHookePotential against central-difference directional derivatives.
auto testIncompressibleNeoHookeSecondOrderDerived()
{
  using namespace Marmot::ContinuumMechanics::EnergyDensityFunctions::Incompressible;

  Tensor33d    C  = get< 0 >( computationParameters() );
  const double mu = 1000.;

  throwExceptionOnFailure( checkIfEqual( get< 0 >( SecondOrderDerived::NeoHookePotential( C, mu ) ),
                                         NeoHookePotential( C, mu ),
                                         1e-14 ),
                           MakeString() << __PRETTY_FUNCTION__ << " energy density mismatch between the bare "
                                        << "potential and its SecondOrderDerived variant." );

  checkSecondOrderDerivedAgainstFD( __PRETTY_FUNCTION__, C, [mu]( const Tensor33d& C_ ) {
    return SecondOrderDerived::NeoHookePotential( C_, mu );
  } );
}

// Verifies the analytic first and second derivative of the (purely isochoric)
// Incompressible::MooneyRivlinPotential against central-difference directional derivatives.
auto testIncompressibleMooneyRivlinSecondOrderDerived()
{
  using namespace Marmot::ContinuumMechanics::EnergyDensityFunctions::Incompressible;

  Tensor33d    C  = get< 0 >( computationParameters() );
  const double C1 = 800.;
  const double C2 = 200.;

  throwExceptionOnFailure( checkIfEqual( get< 0 >( SecondOrderDerived::MooneyRivlinPotential( C, C1, C2 ) ),
                                         MooneyRivlinPotential( C, C1, C2 ),
                                         1e-14 ),
                           MakeString() << __PRETTY_FUNCTION__ << " energy density mismatch between the bare "
                                        << "potential and its SecondOrderDerived variant." );

  checkSecondOrderDerivedAgainstFD( __PRETTY_FUNCTION__, C, [C1, C2]( const Tensor33d& C_ ) {
    return SecondOrderDerived::MooneyRivlinPotential( C_, C1, C2 );
  } );
}

int main()
{
  auto tests = std::vector< std::function< void() > >{ testPenceGouPotentialA,
                                                       testPenceGouPotentialB,
                                                       testPenceGouPotentialC,
                                                       testFirstOrderDerivedB,
                                                       testSecondOrderDerivedB,
                                                       testIncompressibleNeoHookeSecondOrderDerived,
                                                       testIncompressibleMooneyRivlinSecondOrderDerived };

  executeTestsAndCollectExceptions( tests );
  return 0;
}
