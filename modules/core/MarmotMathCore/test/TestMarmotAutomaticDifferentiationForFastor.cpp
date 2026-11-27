#include "Marmot/MarmotAutomaticDifferentiationForFastor.h"
#include "Marmot/MarmotDeformationMeasures.h"
#include "Marmot/MarmotEnergyDensityFunctions.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotTesting.h"

using namespace Fastor;
using namespace Marmot::Math::Differentiation::Automatic;
using namespace Marmot::Testing;
using namespace Marmot::Math::Tensors::FastorTensors;
using namespace Marmot::ContinuumMechanics;

void testTensorToScalar()
{

  // initialise F with identity tensor
  Tensor33d F;
  F.eye();
  Tensor33d C = Marmot::ContinuumMechanics::Kinematics::DeformationMeasures::rightCauchyGreen( F );

  // set material parameters
  const double K = 3500;
  const double G = 1000;

  std::function< JacobianDual< 9 >( const Fastor::Tensor< JacobianDual< 9 >, 3, 3 >& ) > psi =
    [&]( const Fastor::Tensor< JacobianDual< 9 >, 3, 3 >& Ce_ ) {
      return Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions::PenceGouPotentialB( Ce_, K, G );
    };

  // autodiff solution
  Tensor33d dPsi_dC = df_dT( psi, C );

  // check results for dPsi_dC
  throwExceptionOnFailure( checkIfEqual< double >( dPsi_dC, Tensor33d( 0.0 ) ),
                           Marmot::Utilities::MakeString()
                             << __PRETTY_FUNCTION__ << "dPsi_dC should be zero for the F = I" );

  // extension in first direction
  F( 0, 0 ) += 1e-3;
  C = Marmot::ContinuumMechanics::Kinematics::DeformationMeasures::rightCauchyGreen( F );

  // analytical solution
  auto [psi_, dPsi_dC_analytical] = Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions::FirstOrder::
    PenceGouPotentialB( C, K, G );

  // autodiff solution
  dPsi_dC = df_dT( psi, C );

  // check results for dPsi_dC
  throwExceptionOnFailure( checkIfEqual< double >( dPsi_dC, dPsi_dC_analytical, 1e-12 ),
                           Marmot::Utilities::MakeString()
                             << __PRETTY_FUNCTION__ << "dPsi_dC for odometric extension failed" );

  // shear deformation
  F( 0, 1 ) += 1e-3;
  C = Marmot::ContinuumMechanics::Kinematics::DeformationMeasures::rightCauchyGreen( F );

  // analytical solution
  std::tie( psi_, dPsi_dC_analytical ) = Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions::FirstOrder::
    PenceGouPotentialB( C, K, G );

  // autodiff solution
  dPsi_dC = df_dT( psi, C );

  // check results for dPsi_dC
  throwExceptionOnFailure( checkIfEqual< double >( dPsi_dC, dPsi_dC_analytical, 1e-12 ),
                           Marmot::Utilities::MakeString()
                             << __PRETTY_FUNCTION__ << "dPsi_dC for mixed deformation failed" );
}

void testTensorToScalarSymmetric()
{
  Tensor33d F;
  F.eye();
  F( 0, 0 ) += 1e-3;
  F( 0, 1 ) += 2e-3;
  Tensor33d C = Marmot::ContinuumMechanics::Kinematics::DeformationMeasures::rightCauchyGreen( F );

  const double K = 3500;
  const double G = 1000;

  std::function< JacobianDual< 9 >( const Fastor::Tensor< JacobianDual< 9 >, 3, 3 >& ) > psi =
    [&]( const Fastor::Tensor< JacobianDual< 9 >, 3, 3 >& Ce_ ) {
      return Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions::PenceGouPotentialB( Ce_, K, G );
    };

  Tensor33d dPsi_dC_full = df_dT( psi, C );
  Tensor33d dPsi_dC_sym  = df_dT( psi, C, true );

  throwExceptionOnFailure( checkIfEqual< double >( dPsi_dC_sym, dPsi_dC_full, 1e-12 ),
                           Marmot::Utilities::MakeString()
                             << __PRETTY_FUNCTION__ << "symmetric df_dT doesn't match the dense computation" );
}

void testTensorToScalarWith2ndOrderDuals()
{

  // initialise F with identity tensor
  Tensor33d F;
  F.eye();
  Tensor33t< dual > F_dual = Marmot::Math::Tensors::makeDual( F );

  // seed F_dual to check if shifting to higher order duals works
  seed< 1 >( F_dual( 0, 0 ), 1 );

  // compute right Cauchy-Green tensor
  Tensor33t< dual > C_dual = Marmot::ContinuumMechanics::Kinematics::DeformationMeasures::rightCauchyGreen( F_dual );

  // set material parameters
  const double K = 3500;
  const double G = 1000;

  std::function< dual2nd( const Tensor33t< dual2nd >& ) > psi = [&]( const Tensor33t< dual2nd >& Ce_ ) {
    return Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions::PenceGouPotentialB( Ce_, K, G );
  };

  // analytical solution
  auto [psi_, dPsi_dC_analytical] = Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions::FirstOrder::
    PenceGouPotentialB( C_dual, K, G );

  // autodiff solution
  Tensor33t< dual > dPsi_dC = df_dT< 1 >( psi, C_dual ).second;

  // check results for dPsi_dC
  throwExceptionOnFailure( checkIfEqual< dual >( dPsi_dC, dPsi_dC_analytical, 1e-12 ),
                           Marmot::Utilities::MakeString() << __PRETTY_FUNCTION__ << "dPsi_dPsi for F = I failed" );

  // extension in first direction
  F_dual( 0, 0 ) += 1e-3;
  C_dual = Marmot::ContinuumMechanics::Kinematics::DeformationMeasures::rightCauchyGreen( F_dual );

  // analytical solution
  std::tie( psi_, dPsi_dC_analytical ) = Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions::FirstOrder::
    PenceGouPotentialB( C_dual, K, G );

  // autodiff solution
  dPsi_dC = df_dT< 1 >( psi, C_dual ).second;

  // check results for dPsi_dC
  throwExceptionOnFailure( checkIfEqual< dual >( dPsi_dC, dPsi_dC_analytical, 1e-12 ),
                           Marmot::Utilities::MakeString()
                             << __PRETTY_FUNCTION__ << "dPsi_dC for odometric extension failed" );

  // shear deformation
  F_dual( 0, 1 ) += 1e-3;
  C_dual = Marmot::ContinuumMechanics::Kinematics::DeformationMeasures::rightCauchyGreen( F_dual );

  // analytical solution
  std::tie( psi_, dPsi_dC_analytical ) = Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions::FirstOrder::
    PenceGouPotentialB( C_dual, K, G );

  // autodiff solution
  dPsi_dC = df_dT< 1 >( psi, C_dual ).second;

  throwExceptionOnFailure( checkIfEqual< dual >( dPsi_dC, dPsi_dC_analytical, 1e-12 ),
                           Marmot::Utilities::MakeString()
                             << __PRETTY_FUNCTION__ << "dPsi_dC for mixed deformation failed" );
}

void testTensorToScalarWith2ndOrderDualsSymmetric()
{
  Tensor33d F;
  F.eye();
  F( 0, 0 ) += 1e-3;
  F( 0, 1 ) += 2e-3;
  Tensor33d         C      = Marmot::ContinuumMechanics::Kinematics::DeformationMeasures::rightCauchyGreen( F );
  Tensor33t< dual > C_dual = Marmot::Math::Tensors::makeDual( C );

  const double K = 3500;
  const double G = 1000;

  std::function< dual2nd( const Tensor33t< dual2nd >& ) > psi = [&]( const Tensor33t< dual2nd >& Ce_ ) {
    return Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions::PenceGouPotentialB( Ce_, K, G );
  };

  auto [val_full, dPsi_dC_full] = df_dT< 1 >( psi, C_dual );
  auto [val_sym, dPsi_dC_sym]   = df_dT< 1 >( psi, C_dual, true );

  throwExceptionOnFailure( checkIfEqual( double( val_sym ), double( val_full ), 1e-12 ),
                           Marmot::Utilities::MakeString() << __PRETTY_FUNCTION__ << "function value mismatch" );
  throwExceptionOnFailure( checkIfEqual< dual >( dPsi_dC_sym, dPsi_dC_full, 1e-12 ),
                           Marmot::Utilities::MakeString()
                             << __PRETTY_FUNCTION__
                             << "symmetric higher-order df_dT doesn't match the dense computation" );
}

void testTensorToTensor()
{
  Tensor33d F;
  F.random();

  std::function< Tensor33t< JacobianDual< 9 > >( const Tensor33t< JacobianDual< 9 > >& ) > f =
    [&]( const Fastor::Tensor< JacobianDual< 9 >, 3, 3 >& F_ ) {
      Tensor33t< JacobianDual< 9 > > C = Marmot::ContinuumMechanics::Kinematics::DeformationMeasures::rightCauchyGreen(
        F_ );
      return C;
    };

  // autodiff solution
  auto [C, dC_dF] = dF_dT( f, F );

  // analytical solution
  auto [C_analytical, dC_dF_analytical] = Marmot::ContinuumMechanics::Kinematics::DeformationMeasures::FirstOrder::
    rightCauchyGreen( F );

  // check results for C
  throwExceptionOnFailure( checkIfEqual< double >( C, C_analytical, 1e-12 ),
                           Marmot::Utilities::MakeString() << __PRETTY_FUNCTION__ << "computation of C failed" );

  // check results for dC_dF
  throwExceptionOnFailure( checkIfEqual< double >( dC_dF, dC_dF_analytical, 1e-12 ),
                           Marmot::Utilities::MakeString() << __PRETTY_FUNCTION__ << "computation of dC_dF failed" );
}

void testTensorToTensorSymmetric()
{
  Tensor33d F;
  F.eye();
  F( 0, 0 ) += 1e-3;
  F( 0, 1 ) += 2e-3;
  Tensor33d C = Marmot::ContinuumMechanics::Kinematics::DeformationMeasures::rightCauchyGreen( F );

  // f(C) = C*C: transpose-equivariant (f(A^T) = f(A)^T for any A), same-dimension square output, so the
  // symmetry-exploiting path (which relies on that property) is applicable
  std::function< Tensor33t< JacobianDual< 9 > >( const Tensor33t< JacobianDual< 9 > >& ) > f =
    []( const Fastor::Tensor< JacobianDual< 9 >, 3, 3 >& C_ ) {
      return Fastor::einsum< Marmot::Math::Tensors::FastorIndices::ij, Marmot::Math::Tensors::FastorIndices::jk >( C_,
                                                                                                                   C_ );
    };

  auto [Fval_full, dF_dC_full] = dF_dT( f, C );
  auto [Fval_sym, dF_dC_sym]   = dF_dT( f, C, true );

  throwExceptionOnFailure( checkIfEqual< double >( Fval_sym, Fval_full, 1e-12 ),
                           Marmot::Utilities::MakeString() << __PRETTY_FUNCTION__ << "function value mismatch" );
  throwExceptionOnFailure( checkIfEqual< double >( dF_dC_sym, dF_dC_full, 1e-12 ),
                           Marmot::Utilities::MakeString()
                             << __PRETTY_FUNCTION__ << "symmetric dF_dT doesn't match the dense computation" );
}

// dF_dSymmetricT returns the gradient along symmetric directions, i.e. the gradient of the dense computation
// symmetrized in its last index pair, for an output of any shape and a function that is not
// transpose-equivariant.
void testTensorToTensorWrtSymmetricTensor()
{
  Tensor33d F;
  F.eye();
  F( 0, 0 ) += 1e-1;
  F( 0, 1 ) += 2e-1;
  F( 1, 2 ) -= 1e-1;
  const Tensor33d C = Marmot::ContinuumMechanics::Kinematics::DeformationMeasures::rightCauchyGreen( F );

  // rank-4 output C (x) (B C) with a non-symmetric B: not transpose-equivariant
  Tensor33d B( 0.0 );
  B( 0, 0 )    = 1.0;
  B( 0, 1 )    = 2.0;
  B( 1, 2 )    = -1.0;
  B( 2, 0 )    = 0.5;
  B( 2, 2 )    = 3.0;
  const auto f = [&B]< typename Scalar >( const Tensor33t< Scalar >& C_ ) {
    Fastor::Tensor< Scalar, 3, 3, 3, 3 > out;
    for ( int i = 0; i < 3; i++ )
      for ( int j = 0; j < 3; j++ )
        for ( int k = 0; k < 3; k++ )
          for ( int l = 0; l < 3; l++ ) {
            Scalar BC_kl = 0.0;
            for ( int m = 0; m < 3; m++ )
              BC_kl += B( k, m ) * C_( m, l );
            out( i, j, k, l ) = C_( i, j ) * BC_kl;
          }
    return out;
  };

  // the dense reference by dF_dT, whose number type is JacobianDual< 9 >
  std::function< Fastor::Tensor< JacobianDual< 9 >, 3, 3, 3, 3 >( const Tensor33t< JacobianDual< 9 > >& ) > fDense = f;
  const auto [value_dense, gradient_dense] = dF_dT( fDense, C );
  const auto [value, gradient] = dF_dSymmetricT( [&f]( const Tensor33t< dual >& C_ ) { return f( C_ ); }, C );

  Fastor::Tensor< double, 3, 3, 3, 3, 3, 3 > gradient_expected;
  for ( int i = 0; i < 3; i++ )
    for ( int j = 0; j < 3; j++ )
      for ( int k = 0; k < 3; k++ )
        for ( int l = 0; l < 3; l++ )
          for ( int m = 0; m < 3; m++ )
            for ( int n = 0; n < 3; n++ )
              gradient_expected( i, j, k, l, m, n ) = 0.5 * ( gradient_dense( i, j, k, l, m, n ) +
                                                              gradient_dense( i, j, k, l, n, m ) );

  throwExceptionOnFailure( checkIfEqual< double >( value, value_dense, 1e-12 ),
                           Marmot::Utilities::MakeString() << __PRETTY_FUNCTION__ << "function value mismatch" );
  throwExceptionOnFailure( checkIfEqual< double >( gradient, gradient_expected, 1e-12 ),
                           Marmot::Utilities::MakeString()
                             << __PRETTY_FUNCTION__
                             << "dF_dSymmetricT doesn't match the symmetrized dense computation" );

  // the symmetrization matters for this function: the dense gradient is not symmetric in its last index pair
  throwExceptionOnFailure( !checkIfEqual< double >( gradient, gradient_dense, 1e-6 ),
                           Marmot::Utilities::MakeString()
                             << __PRETTY_FUNCTION__ << "test function does not distinguish the two gradients" );
}

void testTensorToScalarSecondOrder()
{

  Tensor33d F;
  F.eye();
  Tensor33d C = Marmot::ContinuumMechanics::Kinematics::DeformationMeasures::rightCauchyGreen( F );

  const double K = 3500;
  const double G = 1000;

  SecondOrder::tensor_to_scalar_function_type< 3 > f = [&]( const Tensor33t< HessianDual< 9 > >& C_ ) {
    return Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions::PenceGouPotentialB( C_, K, G );
  };

  auto [psi, dPsi_dC, d2Psi_dC2] = SecondOrder::d2f_dT2( f, C );

  auto [psi_, dPsi_dC_analytical, d2Psi_dC2_analytical] = Marmot::ContinuumMechanics::Elasticity::
    EnergyDensityFunctions::SecondOrder::PenceGouPotentialB( C, K, G );

  throwExceptionOnFailure( checkIfEqual( psi, psi_, 1e-12 ),
                           Marmot::Utilities::MakeString() << __PRETTY_FUNCTION__ << "computation of psi failed" );

  throwExceptionOnFailure( checkIfEqual< double >( dPsi_dC, dPsi_dC_analytical, 1e-12 ),
                           Marmot::Utilities::MakeString() << __PRETTY_FUNCTION__ << "computation of dPsi_dC failed" );

  throwExceptionOnFailure( checkIfEqual< double >( d2Psi_dC2, d2Psi_dC2_analytical, 1e-12 ),
                           Marmot::Utilities::MakeString()
                             << __PRETTY_FUNCTION__ << "computation of d2Psi_dC2 failed" );
}

void testTensorToScalarSecondOrderSymmetric()
{
  Tensor33d F;
  F.eye();
  F( 0, 0 ) += 1e-3;
  F( 1, 2 ) += 2e-3;
  Tensor33d C = Marmot::ContinuumMechanics::Kinematics::DeformationMeasures::rightCauchyGreen( F );

  const double K = 3500;
  const double G = 1000;

  SecondOrder::tensor_to_scalar_function_type< 3 > f = [&]( const Tensor33t< HessianDual< 9 > >& C_ ) {
    return Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions::PenceGouPotentialB( C_, K, G );
  };

  auto [psi_full, dPsi_dC_full, d2Psi_dC2_full] = SecondOrder::d2f_dT2( f, C );
  auto [psi_sym, dPsi_dC_sym, d2Psi_dC2_sym]    = SecondOrder::d2f_dT2( f, C, true );

  throwExceptionOnFailure( checkIfEqual( psi_sym, psi_full, 1e-12 ),
                           Marmot::Utilities::MakeString() << __PRETTY_FUNCTION__ << "psi mismatch" );
  throwExceptionOnFailure( checkIfEqual< double >( dPsi_dC_sym, dPsi_dC_full, 1e-12 ),
                           Marmot::Utilities::MakeString() << __PRETTY_FUNCTION__ << "dPsi_dC mismatch" );
  throwExceptionOnFailure( checkIfEqual< double >( d2Psi_dC2_sym, d2Psi_dC2_full, 1e-12 ),
                           Marmot::Utilities::MakeString()
                             << __PRETTY_FUNCTION__ << "symmetric d2f_dT2 doesn't match the dense computation" );

  auto [psi_, dPsi_dC_analytical, d2Psi_dC2_analytical] = Marmot::ContinuumMechanics::Elasticity::
    EnergyDensityFunctions::SecondOrder::PenceGouPotentialB( C, K, G );
  throwExceptionOnFailure( checkIfEqual( psi_full, psi_, 1e-12 ) &&
                             checkIfEqual< double >( dPsi_dC_full, dPsi_dC_analytical, 1e-12 ) &&
                             checkIfEqual< double >( d2Psi_dC2_full, d2Psi_dC2_analytical, 1e-10 ),
                           Marmot::Utilities::MakeString()
                             << __PRETTY_FUNCTION__ << "d2f_dT2 doesn't match the analytical derivatives" );
}

void testTensorToScalarSecondOrderMixed()
{

  Tensor33d F;
  F.eye();
  F( 1, 2 ) += 1e-4;

  Tensor33d C = Marmot::ContinuumMechanics::Kinematics::DeformationMeasures::rightCauchyGreen( F );

  const double K = 3500;
  const double G = 1000;

  std::function< dual2nd( const Tensor33t< dual2nd >&, const dual2nd ) > f =
    [&]( const Fastor::Tensor< dual2nd, 3, 3 >& C_, const dual2nd omega_ ) {
      const dual2nd psi = Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions::PenceGouPotentialB( C_,
                                                                                                              K,
                                                                                                              G );
      const dual2nd res = ( -pow( omega_, 2. ) + 1. ) * psi;
      return res;
    };

  const double omega = 0.5;

  // autodiff solution
  auto d2Psi_dCdOmega = SecondOrder::d2f_dTensor_dScalar( f, C, omega );

  // analytical solution
  auto [psi, dPsi_dC_analytical] = Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions::FirstOrder::
    PenceGouPotentialB( C, K, G );

  Tensor33d d2Psi_dCdOmega_analytical = -2. * omega * dPsi_dC_analytical;

  throwExceptionOnFailure( checkIfEqual< double >( d2Psi_dCdOmega, d2Psi_dCdOmega_analytical, 1e-12 ),
                           Marmot::Utilities::MakeString()
                             << __PRETTY_FUNCTION__ << "computation of d2Psi_dC_dOmega failed" );
}

void testTensorToScalarSecondOrderMixedSymmetric()
{
  Tensor33d F;
  F.eye();
  F( 1, 2 ) += 1e-4;

  Tensor33d C = Marmot::ContinuumMechanics::Kinematics::DeformationMeasures::rightCauchyGreen( F );

  const double K = 3500;
  const double G = 1000;

  std::function< dual2nd( const Tensor33t< dual2nd >&, const dual2nd ) > f =
    [&]( const Fastor::Tensor< dual2nd, 3, 3 >& C_, const dual2nd omega_ ) {
      const dual2nd psi = Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions::PenceGouPotentialB( C_,
                                                                                                              K,
                                                                                                              G );
      const dual2nd res = ( -pow( omega_, 2. ) + 1. ) * psi;
      return res;
    };

  const double omega = 0.5;

  auto d2Psi_dCdOmega_full = SecondOrder::d2f_dTensor_dScalar( f, C, omega );
  auto d2Psi_dCdOmega_sym  = SecondOrder::d2f_dTensor_dScalar( f, C, omega, true );

  throwExceptionOnFailure( checkIfEqual< double >( d2Psi_dCdOmega_sym, d2Psi_dCdOmega_full, 1e-12 ),
                           Marmot::Utilities::MakeString()
                             << __PRETTY_FUNCTION__
                             << "symmetric d2f_dTensor_dScalar doesn't match the dense computation" );
}

void testTensorToScalarThirdOrder()
{

  Tensor33d F;
  F.eye();
  Tensor33d C = Marmot::ContinuumMechanics::Kinematics::DeformationMeasures::rightCauchyGreen( F );

  // Use a simpler test function: f(C) = sum(C^2)
  std::function< dual3rd( const Tensor33t< dual3rd >& ) > f = [&]( const Fastor::Tensor< dual3rd, 3, 3 >& C_ ) {
    dual3rd result = 0.0;
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
                           Marmot::Utilities::MakeString() << __PRETTY_FUNCTION__ << "psi is not finite (NaN or Inf)" );

  throwExceptionOnFailure( checkIfEqual( psi, expected_psi, 1e-10 ),
                           Marmot::Utilities::MakeString() << __PRETTY_FUNCTION__ << "computation of psi failed" );

  // check that the third derivative tensor has the correct shape (3x3x3x3x3x3)
  throwExceptionOnFailure( d3Psi_dC3.dimension( 0 ) == 3 && d3Psi_dC3.dimension( 1 ) == 3 &&
                             d3Psi_dC3.dimension( 2 ) == 3 && d3Psi_dC3.dimension( 3 ) == 3 &&
                             d3Psi_dC3.dimension( 4 ) == 3 && d3Psi_dC3.dimension( 5 ) == 3,
                           Marmot::Utilities::MakeString()
                             << __PRETTY_FUNCTION__ << "d3Psi_dC3 has incorrect dimensions" );

  throwExceptionOnFailure( checkIfEqual( d3Psi_dC3, Tensor333333d( 0.0 ) ),
                           Marmot::Utilities::MakeString()
                             << __PRETTY_FUNCTION__ << "d3Psi_dC3 is not equal to zero tensor as expected" );
}

namespace {
  /** All derivatives up to third order of F at T, one evaluation per ordered triple of directions (729 for 3x3) */
  std::tuple< double, Tensor33d, Tensor3333d, Tensor333333d > thirdOrderReference(
    const ThirdOrder::tensor_to_scalar_function_type< 3 >& F,
    const Tensor33d&                                       T )
  {
    double               f = 0;
    Tensor33d            d1;
    Tensor3333d          d2;
    Tensor333333d        d3;
    Tensor33t< dual3rd > T_ = Marmot::Math::Tensors::makeHigherOrderDual< 3 >( T );
    dual3rd*             t  = T_.data();
    for ( size_t a = 0; a < 9; a++ )
      for ( size_t b = 0; b < 9; b++ )
        for ( size_t c = 0; c < 9; c++ ) {
          seed< 1 >( t[a], 1.0 );
          seed< 2 >( t[b], 1.0 );
          seed< 3 >( t[c], 1.0 );
          const dual3rd y = F( T_ );
          seed< 3 >( t[c], 0.0 );
          seed< 2 >( t[b], 0.0 );
          seed< 1 >( t[a], 0.0 );
          f                                              = value( y );
          d1( a / 3, a % 3 )                             = derivative< 1 >( y );
          d2( a / 3, a % 3, b / 3, b % 3 )               = derivative< 2 >( y );
          d3( a / 3, a % 3, b / 3, b % 3, c / 3, c % 3 ) = derivative< 3 >( y );
        }
    return { f, d1, d2, d3 };
  }
} // namespace

void testTensorToScalarThirdOrderAgainstReference()
{
  // a function that is not transpose-invariant, at a non-symmetric point: the full driver must match the reference
  Tensor33d T;
  T                                                 = { { 1.1, 0.2, -0.3 }, { 0.15, 0.9, 0.25 }, { -0.1, 0.35, 1.2 } };
  ThirdOrder::tensor_to_scalar_function_type< 3 > f = []( const Tensor33t< dual3rd >& C_ ) {
    return C_( 0, 0 ) * C_( 0, 1 ) * C_( 1, 2 ) + exp( C_( 1, 0 ) ) * C_( 0, 2 ) * C_( 0, 2 ) +
           C_( 2, 1 ) * pow( C_( 1, 1 ), 3.5 ) + log( C_( 2, 2 ) ) * C_( 2, 0 ) * C_( 1, 0 );
  };
  auto [f0, f1, f2, f3] = thirdOrderReference( f, T );
  auto [g0, g1, g2, g3] = ThirdOrder::d3f_dT3( f, T );
  throwExceptionOnFailure( checkIfEqual( g0, f0, 1e-14 ) && checkIfEqual< double >( g1, f1, 1e-13 ) &&
                             checkIfEqual< double >( g2, f2, 1e-13 ) && checkIfEqual< double >( g3, f3, 1e-12 ),
                           Marmot::Utilities::MakeString()
                             << __PRETTY_FUNCTION__ << "d3f_dT3 doesn't match the reference (general function)" );

  // a transpose-invariant energy at a symmetric point: the full and the symmetric driver must match the reference
  Tensor33d F;
  F.eye();
  F( 0, 0 ) += 0.02;
  F( 0, 1 ) += 0.01;
  F( 1, 2 ) -= 0.015;
  const Tensor33d C = Marmot::ContinuumMechanics::Kinematics::DeformationMeasures::rightCauchyGreen( F );
  ThirdOrder::tensor_to_scalar_function_type< 3 > psi = []( const Tensor33t< dual3rd >& C_ ) {
    return Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions::YeohPotential< dual3rd >( C_,
                                                                                                     3500.,
                                                                                                     500.,
                                                                                                     -10.,
                                                                                                     1. );
  };
  auto [r0, r1, r2, r3] = thirdOrderReference( psi, C );
  for ( bool isSymmetric : { false, true } ) {
    auto [h0, h1, h2, h3] = ThirdOrder::d3f_dT3( psi, C, isSymmetric );
    throwExceptionOnFailure( checkIfEqual( h0, r0, 1e-12 ) && checkIfEqual< double >( h1, r1, 1e-10 ) &&
                               checkIfEqual< double >( h2, r2, 1e-10 ) && checkIfEqual< double >( h3, r3, 1e-9 ),
                             Marmot::Utilities::MakeString()
                               << __PRETTY_FUNCTION__
                               << "d3f_dT3 doesn't match the reference (Yeoh energy, isSymmetric " << isSymmetric
                               << ")" );
  }
}

void testTensorToScalarThirdOrderSymmetric()
{
  Tensor33d F;
  F.eye();
  F( 0, 0 ) += 1e-3;
  F( 1, 2 ) += 2e-3;
  Tensor33d C = Marmot::ContinuumMechanics::Kinematics::DeformationMeasures::rightCauchyGreen( F );

  // det(C): a transpose-invariant scalar function with genuinely nonzero third derivatives, unlike the trivial
  // sum(C^2) used above, so the symmetry-exploiting mirroring is actually exercised
  std::function< dual3rd( const Tensor33t< dual3rd >& ) > f = []( const Fastor::Tensor< dual3rd, 3, 3 >& C_ ) {
    return C_( 0, 0 ) * ( C_( 1, 1 ) * C_( 2, 2 ) - C_( 1, 2 ) * C_( 2, 1 ) ) -
           C_( 0, 1 ) * ( C_( 1, 0 ) * C_( 2, 2 ) - C_( 1, 2 ) * C_( 2, 0 ) ) +
           C_( 0, 2 ) * ( C_( 1, 0 ) * C_( 2, 1 ) - C_( 1, 1 ) * C_( 2, 0 ) );
  };

  auto [psi_full, dPsi_dC_full, d2Psi_dC2_full, d3Psi_dC3_full] = ThirdOrder::d3f_dT3( f, C );
  auto [psi_sym, dPsi_dC_sym, d2Psi_dC2_sym, d3Psi_dC3_sym]     = ThirdOrder::d3f_dT3( f, C, true );

  throwExceptionOnFailure( checkIfEqual( psi_sym, psi_full, 1e-10 ),
                           Marmot::Utilities::MakeString() << __PRETTY_FUNCTION__ << "psi mismatch" );
  throwExceptionOnFailure( checkIfEqual< double >( dPsi_dC_sym, dPsi_dC_full, 1e-10 ),
                           Marmot::Utilities::MakeString() << __PRETTY_FUNCTION__ << "dPsi_dC mismatch" );
  throwExceptionOnFailure( checkIfEqual< double >( d2Psi_dC2_sym, d2Psi_dC2_full, 1e-10 ),
                           Marmot::Utilities::MakeString() << __PRETTY_FUNCTION__ << "d2Psi_dC2 mismatch" );
  throwExceptionOnFailure( checkIfEqual< double >( d3Psi_dC3_sym, d3Psi_dC3_full, 1e-10 ),
                           Marmot::Utilities::MakeString()
                             << __PRETTY_FUNCTION__ << "symmetric d3f_dT3 doesn't match the dense computation" );
}

void testNonSquareRank2Shapes()
{
  // regression test: df_dT and dF_dT must still compile and work for tensors that are not a square rank-2 tensor
  // (here a rank-1 vector input and, separately, a rank-1 vector output), i.e. isSymmetric=false (the default)
  // must be usable without triggering the isSymmetric-only code path's compile-time shape requirements

  // df_dT with a non-square-rank2 (rank-1 vector) input: f(T) = sum(T_i^2), gradient = 2*T
  std::function< JacobianDual< 4 >( const Fastor::Tensor< JacobianDual< 4 >, 4 >& ) > f =
    []( const Fastor::Tensor< JacobianDual< 4 >, 4 >& T ) {
      JacobianDual< 4 > result = 0.0;
      for ( size_t i = 0; i < 4; i++ )
        result += T( i ) * T( i );
      return result;
    };

  Fastor::Tensor< double, 4 > T{ 1., 2., 3., 4. };
  Fastor::Tensor< double, 4 > dF_dT_result = df_dT( f, T );
  Fastor::Tensor< double, 4 > dF_dT_target = 2. * T;

  throwExceptionOnFailure( checkIfEqual< double >( dF_dT_result, dF_dT_target, 1e-12 ),
                           Marmot::Utilities::MakeString()
                             << __PRETTY_FUNCTION__ << "df_dT with a rank-1 vector input failed" );

  // dF_dT with a non-square-rank2 (rank-1 vector) output: g(C)_0 = C_00^2, g(C)_1 = C_11 + C_22
  std::function< Fastor::Tensor< JacobianDual< 9 >, 2 >( const Fastor::Tensor< JacobianDual< 9 >, 3, 3 >& ) > g =
    []( const Fastor::Tensor< JacobianDual< 9 >, 3, 3 >& C_ ) {
      Fastor::Tensor< JacobianDual< 9 >, 2 > out( 0.0 );
      out( 0 ) = C_( 0, 0 ) * C_( 0, 0 );
      out( 1 ) = C_( 1, 1 ) + C_( 2, 2 );
      return out;
    };

  Tensor33d C{ { 1., 2., 3. }, { 4., 5., 6. }, { 7., 8., 9. } };
  auto [Gval, dG_dC] = dF_dT( g, C );

  Fastor::Tensor< double, 2 >       Gval_target( 0.0 );
  Fastor::Tensor< double, 2, 3, 3 > dG_dC_target( 0.0 );
  Gval_target( 0 )        = C( 0, 0 ) * C( 0, 0 );
  Gval_target( 1 )        = C( 1, 1 ) + C( 2, 2 );
  dG_dC_target( 0, 0, 0 ) = 2. * C( 0, 0 );
  dG_dC_target( 1, 1, 1 ) = 1.;
  dG_dC_target( 1, 2, 2 ) = 1.;

  throwExceptionOnFailure( checkIfEqual< double >( Gval, Gval_target, 1e-12 ),
                           Marmot::Utilities::MakeString()
                             << __PRETTY_FUNCTION__ << "dF_dT with a rank-1 vector output: value failed" );
  throwExceptionOnFailure( checkIfEqual< double >( dG_dC, dG_dC_target, 1e-12 ),
                           Marmot::Utilities::MakeString()
                             << __PRETTY_FUNCTION__ << "dF_dT with a rank-1 vector output: gradient failed" );
}

int main()
{

  auto tests = std::vector< std::function< void() > >{ testTensorToScalar,
                                                       testTensorToScalarSymmetric,
                                                       testTensorToScalarWith2ndOrderDuals,
                                                       testTensorToScalarWith2ndOrderDualsSymmetric,
                                                       testTensorToTensor,
                                                       testTensorToTensorSymmetric,
                                                       testTensorToTensorWrtSymmetricTensor,
                                                       testTensorToScalarSecondOrder,
                                                       testTensorToScalarSecondOrderSymmetric,
                                                       testTensorToScalarSecondOrderMixed,
                                                       testTensorToScalarSecondOrderMixedSymmetric,
                                                       testTensorToScalarThirdOrder,
                                                       testTensorToScalarThirdOrderSymmetric,
                                                       testTensorToScalarThirdOrderAgainstReference,
                                                       testNonSquareRank2Shapes };

  executeTestsAndCollectExceptions( tests );
}
