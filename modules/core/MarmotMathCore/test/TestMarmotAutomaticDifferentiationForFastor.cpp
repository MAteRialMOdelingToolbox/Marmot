#include "Marmot/MarmotAutomaticDifferentiationForFastor.h"
#include "Marmot/MarmotDeformationMeasures.h"
#include "Marmot/MarmotEnergyDensityFunctions.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotTesting.h"

using namespace Fastor;
using namespace Marmot::AutomaticDifferentiation;
using namespace Marmot::Testing;
using namespace Marmot::FastorStandardTensors;
using namespace Marmot::ContinuumMechanics;

void testTensorToScalar()
{

  // initialise F with identity tensor
  Tensor33d F;
  F.eye();
  Tensor33d C = DeformationMeasures::CauchyGreen( F );

  // set material parameters
  const double K = 3500;
  const double G = 1000;

  std::function< autodiff::dual( const Fastor::Tensor< autodiff::dual, 3, 3 >& ) > psi =
    [&]( const Fastor::Tensor< autodiff::dual, 3, 3 >& Ce_ ) {
      return EnergyDensityFunctions::PenceGouPotentialB( Ce_, K, G );
    };

  // autodiff solution
  Tensor33d dPsi_dC = df_dT( psi, C );

  // check results for dPsi_dC
  throwExceptionOnFailure( checkIfEqual< double >( dPsi_dC, Tensor33d( 0.0 ) ),
                           MakeString() << __PRETTY_FUNCTION__ << "dPsi_dC should be zero for the F = I" );

  // extension in first direction
  F( 0, 0 ) += 1e-3;
  C = DeformationMeasures::CauchyGreen( F );

  // analytical solution
  auto [psi_, dPsi_dC_analytical] = EnergyDensityFunctions::FirstOrderDerived::PenceGouPotentialB( C, K, G );

  // autodiff solution
  dPsi_dC = df_dT( psi, C );

  // check results for dPsi_dC
  throwExceptionOnFailure( checkIfEqual< double >( dPsi_dC, dPsi_dC_analytical, 1e-12 ),
                           MakeString() << __PRETTY_FUNCTION__ << "dPsi_dC for odometric extension failed" );

  // shear deformation
  F( 0, 1 ) += 1e-3;
  C = DeformationMeasures::CauchyGreen( F );

  // analytical solution
  std::tie( psi_, dPsi_dC_analytical ) = EnergyDensityFunctions::FirstOrderDerived::PenceGouPotentialB( C, K, G );

  // autodiff solution
  dPsi_dC = df_dT( psi, C );

  // check results for dPsi_dC
  throwExceptionOnFailure( checkIfEqual< double >( dPsi_dC, dPsi_dC_analytical, 1e-12 ),
                           MakeString() << __PRETTY_FUNCTION__ << "dPsi_dC for mixed deformation failed" );
}

void testTensorToScalarWith2ndOrderDuals()
{

  // initialise F with identity tensor
  Tensor33d F;
  F.eye();
  Tensor33t< autodiff::dual > F_dual = Marmot::makeDual( F );

  // seed F_dual to check if shifting to higher order duals works
  seed< 1 >( F_dual( 0, 0 ), 1 );

  // compute right Cauchy-Green tensor
  Tensor33t< autodiff::dual > C_dual = DeformationMeasures::CauchyGreen( F_dual );

  // set material parameters
  const double K = 3500;
  const double G = 1000;

  std::function< autodiff::dual2nd( const Tensor33t< autodiff::dual2nd >& ) > psi =
    [&]( const Tensor33t< autodiff::dual2nd >& Ce_ ) {
      return EnergyDensityFunctions::PenceGouPotentialB( Ce_, K, G );
    };

  // analytical solution
  auto [psi_, dPsi_dC_analytical] = EnergyDensityFunctions::FirstOrderDerived::PenceGouPotentialB( C_dual, K, G );

  // autodiff solution
  Tensor33t< autodiff::dual > dPsi_dC = df_dT< 1 >( psi, C_dual );

  // check results for dPsi_dC
  throwExceptionOnFailure( checkIfEqual< autodiff::dual >( dPsi_dC, dPsi_dC_analytical, 1e-12 ),
                           MakeString() << __PRETTY_FUNCTION__ << "dPsi_dPsi for F = I failed" );

  // extension in first direction
  F_dual( 0, 0 ) += 1e-3;
  C_dual = DeformationMeasures::CauchyGreen( F_dual );

  // analytical solution
  std::tie( psi_, dPsi_dC_analytical ) = EnergyDensityFunctions::FirstOrderDerived::PenceGouPotentialB( C_dual, K, G );

  // autodiff solution
  dPsi_dC = df_dT< 1 >( psi, C_dual );

  // check results for dPsi_dC
  throwExceptionOnFailure( checkIfEqual< autodiff::dual >( dPsi_dC, dPsi_dC_analytical, 1e-12 ),
                           MakeString() << __PRETTY_FUNCTION__ << "dPsi_dC for odometric extension failed" );

  // shear deformation
  F_dual( 0, 1 ) += 1e-3;
  C_dual = DeformationMeasures::CauchyGreen( F_dual );

  // analytical solution
  std::tie( psi_, dPsi_dC_analytical ) = EnergyDensityFunctions::FirstOrderDerived::PenceGouPotentialB( C_dual, K, G );

  // autodiff solution
  dPsi_dC = df_dT< 1 >( psi, C_dual );

  throwExceptionOnFailure( checkIfEqual< autodiff::dual >( dPsi_dC, dPsi_dC_analytical, 1e-12 ),
                           MakeString() << __PRETTY_FUNCTION__ << "dPsi_dC for mixed deformation failed" );
}

void testTensorToTensor()
{
  Tensor33d F;
  F.random();

  std::function< Tensor33t< autodiff::dual >( const Tensor33t< autodiff::dual >& ) > f =
    [&]( const Fastor::Tensor< autodiff::dual, 3, 3 >& F_ ) {
      Tensor33t< autodiff::dual > C = DeformationMeasures::CauchyGreen( F_ );
      return C;
    };

  // autodiff solution
  auto [C, dC_dF] = dF_dT( f, F );

  // analytical solution
  auto [C_analytical, dC_dF_analytical] = DeformationMeasures::FirstOrderDerived::CauchyGreen( F );

  // check results for C
  throwExceptionOnFailure( checkIfEqual< double >( C, C_analytical, 1e-12 ),
                           MakeString() << __PRETTY_FUNCTION__ << "computation of C failed" );

  // check results for dC_dF
  throwExceptionOnFailure( checkIfEqual< double >( dC_dF, dC_dF_analytical, 1e-12 ),
                           MakeString() << __PRETTY_FUNCTION__ << "computation of dC_dF failed" );
}

void testTensorToScalarSecondOrder()
{

  Tensor33d F;
  F.eye();
  Tensor33d C = DeformationMeasures::CauchyGreen( F );

  const double K = 3500;
  const double G = 1000;

  std::function< autodiff::dual2nd( const Tensor33t< autodiff::dual2nd >& ) > f =
    [&]( const Fastor::Tensor< autodiff::dual2nd, 3, 3 >& C_ ) {
      return EnergyDensityFunctions::PenceGouPotentialB( C_, K, G );
    };

  auto [psi, dPsi_dC, d2Psi_dC2] = SecondOrder::d2f_dT2( f, C );

  auto [psi_,
        dPsi_dC_analytical,
        d2Psi_dC2_analytical] = EnergyDensityFunctions::SecondOrderDerived::PenceGouPotentialB( C, K, G );

  throwExceptionOnFailure( checkIfEqual( psi, psi_, 1e-12 ),
                           MakeString() << __PRETTY_FUNCTION__ << "computation of psi failed" );

  throwExceptionOnFailure( checkIfEqual< double >( dPsi_dC, dPsi_dC_analytical, 1e-12 ),
                           MakeString() << __PRETTY_FUNCTION__ << "computation of dPsi_dC failed" );

  throwExceptionOnFailure( checkIfEqual< double >( d2Psi_dC2, d2Psi_dC2_analytical, 1e-12 ),
                           MakeString() << __PRETTY_FUNCTION__ << "computation of d2Psi_dC2 failed" );
}

void testTensorToScalarSecondOrderMixed()
{

  Tensor33d F;
  F.eye();
  F( 1, 2 ) += 1e-4;

  Tensor33d C = DeformationMeasures::CauchyGreen( F );

  const double K = 3500;
  const double G = 1000;

  std::function< autodiff::dual2nd( const Tensor33t< autodiff::dual2nd >&, const autodiff::dual2nd ) > f =
    [&]( const Fastor::Tensor< autodiff::dual2nd, 3, 3 >& C_, const autodiff::dual2nd omega_ ) {
      const dual2nd psi = EnergyDensityFunctions::PenceGouPotentialB( C_, K, G );
      const dual2nd res = ( -pow( omega_, 2. ) + 1. ) * psi;
      return res;
    };

  const double omega = 0.5;

  // autodiff solution
  auto d2Psi_dCdOmega = SecondOrder::d2f_dTensor_dScalar( f, C, omega );

  // analytical solution
  auto [psi, dPsi_dC_analytical] = EnergyDensityFunctions::FirstOrderDerived::PenceGouPotentialB( C, K, G );

  Tensor33d d2Psi_dCdOmega_analytical = -2. * omega * dPsi_dC_analytical;

  throwExceptionOnFailure( checkIfEqual< double >( d2Psi_dCdOmega, d2Psi_dCdOmega_analytical, 1e-12 ),
                           MakeString() << __PRETTY_FUNCTION__ << "computation of d2Psi_dC_dOmega failed" );
}

int main()
{

  auto tests = std::vector< std::function< void() > >{ testTensorToScalar,
                                                       testTensorToScalarWith2ndOrderDuals,
                                                       testTensorToTensor,
                                                       testTensorToScalarSecondOrder,
                                                       testTensorToScalarSecondOrderMixed };

  executeTestsAndCollectExceptions( tests );
}
