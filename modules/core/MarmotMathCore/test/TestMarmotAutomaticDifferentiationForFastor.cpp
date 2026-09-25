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
  Tensor33d C = DeformationMeasures::rightCauchyGreen( F );

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
  C = DeformationMeasures::rightCauchyGreen( F );

  // analytical solution
  auto [psi_, dPsi_dC_analytical] = EnergyDensityFunctions::FirstOrderDerived::PenceGouPotentialB( C, K, G );

  // autodiff solution
  dPsi_dC = df_dT( psi, C );

  // check results for dPsi_dC
  throwExceptionOnFailure( checkIfEqual< double >( dPsi_dC, dPsi_dC_analytical, 1e-12 ),
                           MakeString() << __PRETTY_FUNCTION__ << "dPsi_dC for odometric extension failed" );

  // shear deformation
  F( 0, 1 ) += 1e-3;
  C = DeformationMeasures::rightCauchyGreen( F );

  // analytical solution
  std::tie( psi_, dPsi_dC_analytical ) = EnergyDensityFunctions::FirstOrderDerived::PenceGouPotentialB( C, K, G );

  // autodiff solution
  dPsi_dC = df_dT( psi, C );

  // check results for dPsi_dC
  throwExceptionOnFailure( checkIfEqual< double >( dPsi_dC, dPsi_dC_analytical, 1e-12 ),
                           MakeString() << __PRETTY_FUNCTION__ << "dPsi_dC for mixed deformation failed" );
}

void testTensorToScalarSymmetric()
{
  Tensor33d F;
  F.eye();
  F( 0, 0 ) += 1e-3;
  F( 0, 1 ) += 2e-3;
  Tensor33d C = DeformationMeasures::rightCauchyGreen( F );

  const double K = 3500;
  const double G = 1000;

  std::function< autodiff::dual( const Fastor::Tensor< autodiff::dual, 3, 3 >& ) > psi =
    [&]( const Fastor::Tensor< autodiff::dual, 3, 3 >& Ce_ ) {
      return EnergyDensityFunctions::PenceGouPotentialB( Ce_, K, G );
    };

  Tensor33d dPsi_dC_full = df_dT( psi, C );
  Tensor33d dPsi_dC_sym  = df_dT( psi, C, true );

  throwExceptionOnFailure( checkIfEqual< double >( dPsi_dC_sym, dPsi_dC_full, 1e-12 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << "symmetric df_dT doesn't match the dense computation" );
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
  Tensor33t< autodiff::dual > C_dual = DeformationMeasures::rightCauchyGreen( F_dual );

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
  Tensor33t< autodiff::dual > dPsi_dC = df_dT< 1 >( psi, C_dual ).second;

  // check results for dPsi_dC
  throwExceptionOnFailure( checkIfEqual< autodiff::dual >( dPsi_dC, dPsi_dC_analytical, 1e-12 ),
                           MakeString() << __PRETTY_FUNCTION__ << "dPsi_dPsi for F = I failed" );

  // extension in first direction
  F_dual( 0, 0 ) += 1e-3;
  C_dual = DeformationMeasures::rightCauchyGreen( F_dual );

  // analytical solution
  std::tie( psi_, dPsi_dC_analytical ) = EnergyDensityFunctions::FirstOrderDerived::PenceGouPotentialB( C_dual, K, G );

  // autodiff solution
  dPsi_dC = df_dT< 1 >( psi, C_dual ).second;

  // check results for dPsi_dC
  throwExceptionOnFailure( checkIfEqual< autodiff::dual >( dPsi_dC, dPsi_dC_analytical, 1e-12 ),
                           MakeString() << __PRETTY_FUNCTION__ << "dPsi_dC for odometric extension failed" );

  // shear deformation
  F_dual( 0, 1 ) += 1e-3;
  C_dual = DeformationMeasures::rightCauchyGreen( F_dual );

  // analytical solution
  std::tie( psi_, dPsi_dC_analytical ) = EnergyDensityFunctions::FirstOrderDerived::PenceGouPotentialB( C_dual, K, G );

  // autodiff solution
  dPsi_dC = df_dT< 1 >( psi, C_dual ).second;

  throwExceptionOnFailure( checkIfEqual< autodiff::dual >( dPsi_dC, dPsi_dC_analytical, 1e-12 ),
                           MakeString() << __PRETTY_FUNCTION__ << "dPsi_dC for mixed deformation failed" );
}

void testTensorToTensor()
{
  Tensor33d F;
  F.random();

  std::function< Tensor33t< autodiff::dual >( const Tensor33t< autodiff::dual >& ) > f =
    [&]( const Fastor::Tensor< autodiff::dual, 3, 3 >& F_ ) {
      Tensor33t< autodiff::dual > C = DeformationMeasures::rightCauchyGreen( F_ );
      return C;
    };

  // autodiff solution
  auto [C, dC_dF] = dF_dT( f, F );

  // analytical solution
  auto [C_analytical, dC_dF_analytical] = DeformationMeasures::FirstOrderDerived::rightCauchyGreen( F );

  // check results for C
  throwExceptionOnFailure( checkIfEqual< double >( C, C_analytical, 1e-12 ),
                           MakeString() << __PRETTY_FUNCTION__ << "computation of C failed" );

  // check results for dC_dF
  throwExceptionOnFailure( checkIfEqual< double >( dC_dF, dC_dF_analytical, 1e-12 ),
                           MakeString() << __PRETTY_FUNCTION__ << "computation of dC_dF failed" );
}

void testTensorToTensorSymmetric()
{
  Tensor33d F;
  F.eye();
  F( 0, 0 ) += 1e-3;
  F( 0, 1 ) += 2e-3;
  Tensor33d C = DeformationMeasures::rightCauchyGreen( F );

  // f(C) = C*C: transpose-equivariant (f(A^T) = f(A)^T for any A), same-dimension square output, so the
  // symmetry-exploiting path (which relies on that property) is applicable
  std::function< Tensor33t< autodiff::dual >( const Tensor33t< autodiff::dual >& ) > f =
    []( const Fastor::Tensor< autodiff::dual, 3, 3 >& C_ ) {
      return Fastor::einsum< Marmot::FastorIndices::ij, Marmot::FastorIndices::jk >( C_, C_ );
    };

  auto [Fval_full, dF_dC_full] = dF_dT( f, C );
  auto [Fval_sym, dF_dC_sym]   = dF_dT( f, C, true );

  throwExceptionOnFailure( checkIfEqual< double >( Fval_sym, Fval_full, 1e-12 ),
                           MakeString() << __PRETTY_FUNCTION__ << "function value mismatch" );
  throwExceptionOnFailure( checkIfEqual< double >( dF_dC_sym, dF_dC_full, 1e-12 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << "symmetric dF_dT doesn't match the dense computation" );
}

void testTensorToScalarSecondOrder()
{

  Tensor33d F;
  F.eye();
  Tensor33d C = DeformationMeasures::rightCauchyGreen( F );

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

void testTensorToScalarSecondOrderSymmetric()
{
  Tensor33d F;
  F.eye();
  F( 0, 0 ) += 1e-3;
  F( 1, 2 ) += 2e-3;
  Tensor33d C = DeformationMeasures::rightCauchyGreen( F );

  const double K = 3500;
  const double G = 1000;

  std::function< autodiff::dual2nd( const Tensor33t< autodiff::dual2nd >& ) > f =
    [&]( const Fastor::Tensor< autodiff::dual2nd, 3, 3 >& C_ ) {
      return EnergyDensityFunctions::PenceGouPotentialB( C_, K, G );
    };

  auto [psi_full, dPsi_dC_full, d2Psi_dC2_full] = SecondOrder::d2f_dT2( f, C );
  auto [psi_sym, dPsi_dC_sym, d2Psi_dC2_sym]    = SecondOrder::d2f_dT2( f, C, true );

  throwExceptionOnFailure( checkIfEqual( psi_sym, psi_full, 1e-12 ),
                           MakeString() << __PRETTY_FUNCTION__ << "psi mismatch" );
  throwExceptionOnFailure( checkIfEqual< double >( dPsi_dC_sym, dPsi_dC_full, 1e-12 ),
                           MakeString() << __PRETTY_FUNCTION__ << "dPsi_dC mismatch" );
  throwExceptionOnFailure( checkIfEqual< double >( d2Psi_dC2_sym, d2Psi_dC2_full, 1e-12 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << "symmetric d2f_dT2 doesn't match the dense computation" );
}

void testTensorToScalarSecondOrderMixed()
{

  Tensor33d F;
  F.eye();
  F( 1, 2 ) += 1e-4;

  Tensor33d C = DeformationMeasures::rightCauchyGreen( F );

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

void testTensorToScalarSecondOrderMixedSymmetric()
{
  Tensor33d F;
  F.eye();
  F( 1, 2 ) += 1e-4;

  Tensor33d C = DeformationMeasures::rightCauchyGreen( F );

  const double K = 3500;
  const double G = 1000;

  std::function< autodiff::dual2nd( const Tensor33t< autodiff::dual2nd >&, const autodiff::dual2nd ) > f =
    [&]( const Fastor::Tensor< autodiff::dual2nd, 3, 3 >& C_, const autodiff::dual2nd omega_ ) {
      const dual2nd psi = EnergyDensityFunctions::PenceGouPotentialB( C_, K, G );
      const dual2nd res = ( -pow( omega_, 2. ) + 1. ) * psi;
      return res;
    };

  const double omega = 0.5;

  auto d2Psi_dCdOmega_full = SecondOrder::d2f_dTensor_dScalar( f, C, omega );
  auto d2Psi_dCdOmega_sym  = SecondOrder::d2f_dTensor_dScalar( f, C, omega, true );

  throwExceptionOnFailure( checkIfEqual< double >( d2Psi_dCdOmega_sym, d2Psi_dCdOmega_full, 1e-12 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << "symmetric d2f_dTensor_dScalar doesn't match the dense computation" );
}

void testTensorToScalarThirdOrder()
{

  Tensor33d F;
  F.eye();
  Tensor33d C = DeformationMeasures::rightCauchyGreen( F );

  // Use a simpler test function: f(C) = sum(C^2)
  std::function< autodiff::dual3rd( const Tensor33t< autodiff::dual3rd >& ) > f =
    [&]( const Fastor::Tensor< autodiff::dual3rd, 3, 3 >& C_ ) {
      autodiff::dual3rd result = 0.0;
      for ( size_t i = 0; i < 3; i++ ) {
        for ( size_t j = 0; j < 3; j++ ) {
          result += C_( i, j ) * C_( i, j );
        }
      }
      return result;
    };

  auto [psi, dPsi_dC, d2Psi_dC2, d3Psi_dC3] = ThirdOrder::d3f_dT3( f, C );

  double expected_psi = 3.0;

  throwExceptionOnFailure( std::isfinite( psi ),
                           MakeString() << __PRETTY_FUNCTION__ << "psi is not finite (NaN or Inf)" );

  throwExceptionOnFailure( checkIfEqual( psi, expected_psi, 1e-10 ),
                           MakeString() << __PRETTY_FUNCTION__ << "computation of psi failed" );

  // check that the third derivative tensor has the correct shape (3x3x3x3x3x3)
  throwExceptionOnFailure( d3Psi_dC3.dimension( 0 ) == 3 && d3Psi_dC3.dimension( 1 ) == 3 &&
                             d3Psi_dC3.dimension( 2 ) == 3 && d3Psi_dC3.dimension( 3 ) == 3 &&
                             d3Psi_dC3.dimension( 4 ) == 3 && d3Psi_dC3.dimension( 5 ) == 3,
                           MakeString() << __PRETTY_FUNCTION__ << "d3Psi_dC3 has incorrect dimensions" );

  throwExceptionOnFailure( checkIfEqual( d3Psi_dC3, Tensor333333d( 0.0 ) ),
                           MakeString() << __PRETTY_FUNCTION__ << "d3Psi_dC3 is not equal to zero tensor as expected" );
}

void testTensorToScalarThirdOrderSymmetric()
{
  Tensor33d F;
  F.eye();
  F( 0, 0 ) += 1e-3;
  F( 1, 2 ) += 2e-3;
  Tensor33d C = DeformationMeasures::rightCauchyGreen( F );

  // det(C): a transpose-invariant scalar function with genuinely nonzero third derivatives, unlike the trivial
  // sum(C^2) used above, so the symmetry-exploiting mirroring is actually exercised
  std::function< autodiff::dual3rd( const Tensor33t< autodiff::dual3rd >& ) > f =
    []( const Fastor::Tensor< autodiff::dual3rd, 3, 3 >& C_ ) {
      return C_( 0, 0 ) * ( C_( 1, 1 ) * C_( 2, 2 ) - C_( 1, 2 ) * C_( 2, 1 ) ) -
             C_( 0, 1 ) * ( C_( 1, 0 ) * C_( 2, 2 ) - C_( 1, 2 ) * C_( 2, 0 ) ) +
             C_( 0, 2 ) * ( C_( 1, 0 ) * C_( 2, 1 ) - C_( 1, 1 ) * C_( 2, 0 ) );
    };

  auto [psi_full, dPsi_dC_full, d2Psi_dC2_full, d3Psi_dC3_full] = ThirdOrder::d3f_dT3( f, C );
  auto [psi_sym, dPsi_dC_sym, d2Psi_dC2_sym, d3Psi_dC3_sym]     = ThirdOrder::d3f_dT3( f, C, true );

  throwExceptionOnFailure( checkIfEqual( psi_sym, psi_full, 1e-10 ),
                           MakeString() << __PRETTY_FUNCTION__ << "psi mismatch" );
  throwExceptionOnFailure( checkIfEqual< double >( dPsi_dC_sym, dPsi_dC_full, 1e-10 ),
                           MakeString() << __PRETTY_FUNCTION__ << "dPsi_dC mismatch" );
  throwExceptionOnFailure( checkIfEqual< double >( d2Psi_dC2_sym, d2Psi_dC2_full, 1e-10 ),
                           MakeString() << __PRETTY_FUNCTION__ << "d2Psi_dC2 mismatch" );
  throwExceptionOnFailure( checkIfEqual< double >( d3Psi_dC3_sym, d3Psi_dC3_full, 1e-10 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << "symmetric d3f_dT3 doesn't match the dense computation" );
}

int main()
{

  auto tests = std::vector< std::function< void() > >{ testTensorToScalar,
                                                       testTensorToScalarSymmetric,
                                                       testTensorToScalarWith2ndOrderDuals,
                                                       testTensorToTensor,
                                                       testTensorToTensorSymmetric,
                                                       testTensorToScalarSecondOrder,
                                                       testTensorToScalarSecondOrderSymmetric,
                                                       testTensorToScalarSecondOrderMixed,
                                                       testTensorToScalarSecondOrderMixedSymmetric,
                                                       testTensorToScalarThirdOrder,
                                                       testTensorToScalarThirdOrderSymmetric };

  executeTestsAndCollectExceptions( tests );
}
