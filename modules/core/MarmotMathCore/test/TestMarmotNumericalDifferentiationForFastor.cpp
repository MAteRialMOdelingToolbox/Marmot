// #include "Marmot/MarmotJournal.h"
#include "Fastor/Fastor.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotNumericalDifferentiationForFastor.h"
#include "Marmot/MarmotTesting.h"
// #include "Marmot/MarmotConstants.h"

using namespace Marmot;
using namespace Marmot::Testing;
using namespace Marmot::NumericalAlgorithms::Differentiation;

void testScalarToTensor()
{
  std::function< Fastor::Tensor< double, 2, 2 >( double ) > scalarToTensor_func = []( const double x ) {
    Fastor::Tensor< double, 2, 2 > res;
    res( 0, 0 ) = 3.4 * x * x * exp( x / 5. ) + 10. * std::cos( x );
    res( 0, 1 ) = 2. * x * x * x;
    res( 1, 0 ) = std::sin( 2. * x );
    res( 1, 1 ) = 43.0;
    return res;
  };

  const double x_ = 4.5;

  Fastor::Tensor< double, 2, 2 > dF_dX_forward = ScalarToTensor::forwardDifference( scalarToTensor_func, x_ );
  Fastor::Tensor< double, 2, 2 > dF_dX_central = ScalarToTensor::centralDifference( scalarToTensor_func, x_ );
  Fastor::Tensor< double, 2, 2 > dF_dX_target;
  dF_dX_target( 0, 0 ) = 3.4 * 2. * x_ * exp( x_ / 5. ) + 3.4 * x_ * x_ * exp( x_ / 5. ) / 5. - 10. * std::sin( x_ );
  dF_dX_target( 0, 1 ) = 2. * 3. * x_ * x_;
  dF_dX_target( 1, 0 ) = std::cos( 2. * x_ ) * 2.;
  dF_dX_target( 1, 1 ) = 0.0;

  throwExceptionOnFailure( Fastor::isequal( dF_dX_forward, dF_dX_target, 1e-5 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " failed: forward difference doesn't yield the right result" );

  throwExceptionOnFailure( Fastor::isequal( dF_dX_central, dF_dX_target, 1e-9 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " failed: central difference doesn't yield the right result" );
}

void testTensorToScalar()
{
  std::function< double( const Fastor::Tensor< double, 2, 2 >& ) > tensorToScalar_func =
    []( const Fastor::Tensor< double, 2, 2 >& x ) {
      // auto   frobeniusNorm_squared = Fastor::einsum< Fastor::Index< 0, 1 >, Fastor::Index< 0, 1 > >( x, x );
      auto   frobeniusNorm_squared = Fastor::einsum< FastorIndices::ij, FastorIndices::ij >( x, x );
      double res                   = frobeniusNorm_squared.toscalar();
      return res;
    };

  const Fastor::Tensor< double, 2, 2 > x_{ { 1., 2. }, { 3., 4. } };

  Fastor::Tensor< double, 2, 2 > dF_dX_forward = TensorToScalar::forwardDifference( tensorToScalar_func, x_ );
  Fastor::Tensor< double, 2, 2 > dF_dX_central = TensorToScalar::centralDifference( tensorToScalar_func, x_ );
  Fastor::Tensor< double, 2, 2 > dF_dX_target  = 2. * x_;

  throwExceptionOnFailure( Fastor::isequal( dF_dX_forward, dF_dX_target, 1e-7 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " failed: forward difference doesn't yield the right result" );

  throwExceptionOnFailure( Fastor::isequal( dF_dX_central, dF_dX_target, 1e-9 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " failed: central difference doesn't yield the right result" );
}

void testTensorToScalarSymmetric()
{
  std::function< double( const Fastor::Tensor< double, 3, 3 >& ) > tensorToScalar_func =
    []( const Fastor::Tensor< double, 3, 3 >& x ) {
      auto   frobeniusNorm_squared = Fastor::einsum< FastorIndices::ij, FastorIndices::ij >( x, x );
      double res                   = frobeniusNorm_squared.toscalar();
      return res;
    };

  // a genuinely symmetric tensor, e.g. akin to a right Cauchy-Green tensor
  const Fastor::Tensor< double, 3, 3 > C{ { 4., 1., 2. }, { 1., 5., 3. }, { 2., 3., 6. } };

  Fastor::Tensor< double, 3, 3 > dF_dC_forward_full = TensorToScalar::forwardDifference( tensorToScalar_func, C );
  Fastor::Tensor< double, 3, 3 > dF_dC_forward_sym  = TensorToScalar::forwardDifference( tensorToScalar_func, C, true );
  Fastor::Tensor< double, 3, 3 > dF_dC_central_full = TensorToScalar::centralDifference( tensorToScalar_func, C );
  Fastor::Tensor< double, 3, 3 > dF_dC_central_sym  = TensorToScalar::centralDifference( tensorToScalar_func, C, true );

  // note: the two computations perturb different (but analytically equivalent) tensor entries, so floating-point
  // summation order can differ slightly; tolerances are chosen a couple orders of magnitude above the expected
  // truncation error of each method rather than machine precision
  throwExceptionOnFailure( Fastor::isequal( dF_dC_forward_sym, dF_dC_forward_full, 1e-7 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " failed: symmetric forward difference doesn't match the dense result" );

  throwExceptionOnFailure( Fastor::isequal( dF_dC_central_sym, dF_dC_central_full, 1e-7 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " failed: symmetric central difference doesn't match the dense result" );
}

void testTensorToTensor()
{
  std::function< Fastor::Tensor< double, 3, 3 >( const Fastor::Tensor< double, 3, 3 >& ) > tensorToTensor_func =
    []( const Fastor::Tensor< double, 3, 3 >& x ) {
      auto matrixproduct = Fastor::einsum< FastorIndices::ij, FastorIndices::jk >( x, x );
      return matrixproduct;
    };

  const Fastor::Tensor< double, 3, 3 > x_{ { 1., 2., 3. }, { 4., 5., 6. }, { 7., 8., 9. } };

  Fastor::Tensor< double, 3, 3, 3, 3 > dF_dX_forward = TensorToTensor::forwardDifference( tensorToTensor_func, x_ );
  Fastor::Tensor< double, 3, 3, 3, 3 > dF_dX_central = TensorToTensor::centralDifference( tensorToTensor_func, x_ );
  Fastor::Tensor< double, 3, 3, 3, 3 >
    dF_dX_target = Fastor::einsum< Fastor::Index< 0, 2 >,
                                   Fastor::Index< 3, 1 >,
                                   Fastor::OIndex< 0, 1, 2, 3 > >( FastorStandardTensors::Spatial3D::I, x_ ) +
                   Fastor::einsum< Fastor::Index< 0, 2 >,
                                   Fastor::Index< 1, 3 >,
                                   Fastor::OIndex< 0, 1, 2, 3 > >( x_, FastorStandardTensors::Spatial3D::I );

  throwExceptionOnFailure( Fastor::isequal( dF_dX_forward, dF_dX_target, 1e-6 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " failed: forward difference doesn't yield the right result" );

  throwExceptionOnFailure( Fastor::isequal( dF_dX_central, dF_dX_target, 1e-8 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " failed: central difference doesn't yield the right result" );
}

Fastor::Tensor< double, 3, 3, 3, 3 > createRank4StiffnessTensor( double E, double nu )
{
  /* lamé parameters */
  double                               lambda = nu * E / ( ( 1. + nu ) * ( 1. - 2. * nu ) );
  double                               G      = E / ( 2. * ( 1. + nu ) );
  Fastor::Tensor< double, 3, 3, 3, 3 > C;
  for ( int i = 0; i < 3; ++i ) {
    for ( int j = 0; j < 3; ++j ) {
      for ( int k = 0; k < 3; ++k ) {
        for ( int l = 0; l < 3; ++l ) {
          C( i, j, k, l ) = lambda * ( i == j ) * ( k == l ) +
                            G * ( ( i == k ) * ( j == l ) + ( i == l ) * ( j == k ) );
        }
      }
    }
  }
  return C;
}

Fastor::Tensor< double, 3, 3 > hookes_law( const Fastor::Tensor< double, 3, 3 >& strain )
{
  Fastor::Tensor< double, 3, 3, 3, 3 > C      = createRank4StiffnessTensor( 30000., 0.2 );
  Fastor::Tensor< double, 3, 3 >       stress = Fastor::einsum< FastorIndices::ijkl, FastorIndices::kl >( C, strain );
  return stress;
}

void testTensorToTensor2()
{
  std::function< Fastor::Tensor< double, 3, 3 >( const Fastor::Tensor< double, 3, 3 >& ) >
    tensorToTensor_func = hookes_law;

  const Fastor::Tensor< double, 3, 3 > strain{ { 1., 2., 3. }, { 4., 5., 6. }, { 7., 8., 9. } };

  const Fastor::Tensor< double, 3, 3, 3, 3 > dF_dX_forward = TensorToTensor::forwardDifference( tensorToTensor_func,
                                                                                                strain );
  const Fastor::Tensor< double, 3, 3, 3, 3 > dF_dX_central = TensorToTensor::centralDifference( tensorToTensor_func,
                                                                                                strain );
  const Fastor::Tensor< double, 3, 3, 3, 3 > dF_dX_target  = createRank4StiffnessTensor( 30000., 0.2 );

  throwExceptionOnFailure( Fastor::isequal( dF_dX_forward, dF_dX_target, 3e-3 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " failed: forward difference doesn't yield the right result" );

  throwExceptionOnFailure( Fastor::isequal( dF_dX_central, dF_dX_target, 1e-5 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " failed: central difference doesn't yield the right result" );
}

void testTensorToTensorSymmetric()
{
  std::function< Fastor::Tensor< double, 3, 3 >( const Fastor::Tensor< double, 3, 3 >& ) >
    tensorToTensor_func = hookes_law;

  // a symmetric strain tensor: the resulting stiffness (4th order derivative) is symmetric wrt swapping the
  // input tensor's index pair since the underlying stiffness tensor has minor symmetry C_ijkl = C_ijlk
  const Fastor::Tensor< double, 3, 3 > strain{ { 1., 2., 3. }, { 2., 4., 5. }, { 3., 5., 6. } };

  const Fastor::Tensor< double, 3, 3, 3, 3 >
    dF_dX_forward_full = TensorToTensor::forwardDifference( tensorToTensor_func, strain );
  const Fastor::Tensor< double, 3, 3, 3, 3 > dF_dX_forward_sym = TensorToTensor::forwardDifference( tensorToTensor_func,
                                                                                                    strain,
                                                                                                    true );
  const Fastor::Tensor< double, 3, 3, 3, 3 >
    dF_dX_central_full = TensorToTensor::centralDifference( tensorToTensor_func, strain );
  const Fastor::Tensor< double, 3, 3, 3, 3 > dF_dX_central_sym = TensorToTensor::centralDifference( tensorToTensor_func,
                                                                                                    strain,
                                                                                                    true );

  throwExceptionOnFailure( Fastor::isequal( dF_dX_forward_sym, dF_dX_forward_full, 1e-6 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " failed: symmetric forward difference doesn't match the dense result" );

  throwExceptionOnFailure( Fastor::isequal( dF_dX_central_sym, dF_dX_central_full, 1e-6 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " failed: symmetric central difference doesn't match the dense result" );
}

void testScalarToTensorComplex()
{
  std::function< Fastor::Tensor< std::complex< double >, 2, 2 >( std::complex< double > ) > scalarToTensor_func =
    []( const std::complex< double > x ) {
      Fastor::Tensor< std::complex< double >, 2, 2 > res;
      res( 0, 0 ) = 3.4 * x * x * exp( x / 5. ) + 10. * std::cos( x );
      res( 0, 1 ) = 2. * x * x * x;
      res( 1, 0 ) = std::sin( 2. * x );
      res( 1, 1 ) = 43.0;
      return res;
    };

  const double x_ = 4.5;

  auto [F_forward, dF_dX_forward] = Complex::ScalarToTensor::forwardDifference( scalarToTensor_func, x_ );
  Fastor::Tensor< double, 2, 2 > dF_dX_target;
  dF_dX_target( 0, 0 ) = 3.4 * 2. * x_ * exp( x_ / 5. ) + 3.4 * x_ * x_ * exp( x_ / 5. ) / 5. - 10. * std::sin( x_ );
  dF_dX_target( 0, 1 ) = 2. * 3. * x_ * x_;
  dF_dX_target( 1, 0 ) = std::cos( 2. * x_ ) * 2.;
  dF_dX_target( 1, 1 ) = 0.0;
  Fastor::Tensor< std::complex< double >, 2, 2 > F_target = scalarToTensor_func( x_ );

  auto convert_to_real =
    []( const Fastor::Tensor< std::complex< double >, 2, 2 >& complex_tensor ) -> Fastor::Tensor< double, 2, 2 > {
    Fastor::Tensor< double, 2, 2 > real_tensor;
    for ( Fastor::FASTOR_INDEX j = 0; j < complex_tensor.size(); ++j ) {
      real_tensor.data()[real_tensor.get_mem_index( j )] = complex_tensor.data()[complex_tensor.get_mem_index( j )]
                                                             .real();
    }
    return real_tensor;
  };

  Fastor::Tensor< double, 2, 2 > F_target_real = convert_to_real( F_target );

  throwExceptionOnFailure( Fastor::isequal( dF_dX_forward, dF_dX_target, 1e-10 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " failed: complex forward difference doesn't yield the right result" );

  throwExceptionOnFailure( Fastor::isequal( F_forward, F_target_real, 1e-10 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " failed: complex forward difference doesn't yield the right result" );
}

void testTensorToScalarComplex()
{
  std::function< std::complex< double >( const Fastor::Tensor< std::complex< double >, 2, 2 >& ) > tensorToScalar_func =
    []( const Fastor::Tensor< std::complex< double >, 2, 2 >& x ) {
      // auto   frobeniusNorm_squared = Fastor::einsum< Fastor::Index< 0, 1 >, Fastor::Index< 0, 1 > >( x, x );
      auto                   frobeniusNorm_squared = Fastor::einsum< FastorIndices::ij, FastorIndices::ij >( x, x );
      std::complex< double > res                   = frobeniusNorm_squared.toscalar();
      return res;
    };

  const Fastor::Tensor< double, 2, 2 > x_{ { 1., 2. }, { 3., 4. } };

  Fastor::Tensor< double, 2, 2 > dF_dX_forward  = Complex::TensorToScalar::forwardDifference( tensorToScalar_func, x_ );
  Fastor::Tensor< double, 2, 2 > dF_dX_forward2 = Complex::forwardDifference( tensorToScalar_func, x_ );
  Fastor::Tensor< double, 2, 2 > dF_dX_target   = 2. * x_;

  throwExceptionOnFailure( Fastor::isequal( dF_dX_forward, dF_dX_target, 1e-10 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " failed: forward difference doesn't yield the right result" );

  throwExceptionOnFailure( Fastor::isequal( dF_dX_forward2, dF_dX_target, 1e-10 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " failed: forward difference (function outside of namespace TensorToScalar) "
                                           "doesn't yield the right result" );
}

void testTensorToScalarComplexSymmetric()
{
  std::function< std::complex< double >( const Fastor::Tensor< std::complex< double >, 3, 3 >& ) > tensorToScalar_func =
    []( const Fastor::Tensor< std::complex< double >, 3, 3 >& x ) {
      auto                   frobeniusNorm_squared = Fastor::einsum< FastorIndices::ij, FastorIndices::ij >( x, x );
      std::complex< double > res                   = frobeniusNorm_squared.toscalar();
      return res;
    };

  const Fastor::Tensor< double, 3, 3 > C{ { 4., 1., 2. }, { 1., 5., 3. }, { 2., 3., 6. } };

  Fastor::Tensor< double, 3, 3 > dF_dC_full = Complex::TensorToScalar::forwardDifference( tensorToScalar_func, C );
  Fastor::Tensor< double, 3, 3 > dF_dC_sym = Complex::TensorToScalar::forwardDifference( tensorToScalar_func, C, true );

  throwExceptionOnFailure( Fastor::isequal( dF_dC_sym, dF_dC_full, 1e-10 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " failed: symmetric complex-step derivative doesn't match the dense "
                                           "result" );
}

void testTensorToTensorComplex()
{
  std::function< Fastor::Tensor< std::complex< double >, 3, 3 >(
    const Fastor::Tensor< std::complex< double >, 3, 3 >& ) >
    tensorToTensor_func = []( const Fastor::Tensor< std::complex< double >, 3, 3 >& x ) {
      auto matrixproduct = Fastor::einsum< FastorIndices::ij, FastorIndices::jk >( x, x );
      return matrixproduct;
    };

  const Fastor::Tensor< double, 3, 3 > x_{ { 1., 2., 3. }, { 4., 5., 6. }, { 7., 8., 9. } };

  Fastor::Tensor< double, 3, 3, 3, 3 > dF_dX_forward = Complex::TensorToTensor::forwardDifference( tensorToTensor_func,
                                                                                                   x_ );
  Fastor::Tensor< double, 3, 3, 3, 3 >
    dF_dX_target = Fastor::einsum< Fastor::Index< 0, 2 >,
                                   Fastor::Index< 3, 1 >,
                                   Fastor::OIndex< 0, 1, 2, 3 > >( FastorStandardTensors::Spatial3D::I, x_ ) +
                   Fastor::einsum< Fastor::Index< 0, 2 >,
                                   Fastor::Index< 1, 3 >,
                                   Fastor::OIndex< 0, 1, 2, 3 > >( x_, FastorStandardTensors::Spatial3D::I );

  throwExceptionOnFailure( Fastor::isequal( dF_dX_forward, dF_dX_target, 1e-10 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " failed: forward difference doesn't yield the right result" );
}

Fastor::Tensor< std::complex< double >, 3, 3 > hookes_law(
  const Fastor::Tensor< std::complex< double >, 3, 3 >& strain )
{
  Fastor::Tensor< double, 3, 3, 3, 3 > C_double = createRank4StiffnessTensor( 30000., 0.2 );
  Fastor::Tensor< std::complex< double >, 3, 3, 3, 3 >
    C_complex = fastorTensorFromDoubleTensor< std::complex< double > >( C_double );
  Fastor::Tensor< std::complex< double >, 3, 3 >
    stress = Fastor::einsum< FastorIndices::ijkl, FastorIndices::kl >( C_complex, strain );
  return stress;
}

void testTensorToTensorComplex2()
{
  // cast hookes_law to a function pointer with the correct signature
  std::function< Fastor::Tensor< std::complex< double >, 3, 3 >(
    const Fastor::Tensor< std::complex< double >, 3, 3 >& ) >
    tensorToTensor_func = static_cast< Fastor::Tensor< std::complex< double >, 3, 3 > ( * )(
      const Fastor::Tensor< std::complex< double >, 3, 3 >& ) >( hookes_law );

  const Fastor::Tensor< double, 3, 3 > strain{ { 1., 2., 3. }, { 4., 5., 6. }, { 7., 8., 9. } };

  const Fastor::Tensor< double, 3, 3, 3, 3 >
    dF_dX_forward = Complex::TensorToTensor::forwardDifference( tensorToTensor_func, strain );
  const Fastor::Tensor< double, 3, 3, 3, 3 > dF_dX_target = createRank4StiffnessTensor( 30000., 0.2 );

  throwExceptionOnFailure( Fastor::isequal( dF_dX_forward, dF_dX_target, 3e-3 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " failed: forward difference doesn't yield the right result" );
}

void testTensorToTensorComplexSymmetric()
{
  std::function< Fastor::Tensor< std::complex< double >, 3, 3 >(
    const Fastor::Tensor< std::complex< double >, 3, 3 >& ) >
    tensorToTensor_func = static_cast< Fastor::Tensor< std::complex< double >, 3, 3 > ( * )(
      const Fastor::Tensor< std::complex< double >, 3, 3 >& ) >( hookes_law );

  const Fastor::Tensor< double, 3, 3 > strain{ { 1., 2., 3. }, { 2., 4., 5. }, { 3., 5., 6. } };

  const Fastor::Tensor< double, 3, 3, 3, 3 >
    dF_dX_full = Complex::TensorToTensor::forwardDifference( tensorToTensor_func, strain );
  const Fastor::Tensor< double, 3, 3, 3, 3 >
    dF_dX_sym = Complex::TensorToTensor::forwardDifference( tensorToTensor_func, strain, true );

  throwExceptionOnFailure( Fastor::isequal( dF_dX_sym, dF_dX_full, 1e-10 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " failed: symmetric complex-step derivative doesn't match the dense "
                                           "result" );
}

int main()
{

  auto tests = std::vector< std::function< void() > >{ testScalarToTensor,
                                                       testTensorToScalar,
                                                       testTensorToScalarSymmetric,
                                                       testTensorToTensor,
                                                       testTensorToTensor2,
                                                       testTensorToTensorSymmetric,
                                                       testScalarToTensorComplex,
                                                       testTensorToScalarComplex,
                                                       testTensorToScalarComplexSymmetric,
                                                       testTensorToTensorComplex,
                                                       testTensorToTensorComplex2,
                                                       testTensorToTensorComplexSymmetric };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
