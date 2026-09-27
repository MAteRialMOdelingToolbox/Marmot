#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotTesting.h"
#include <complex>

using namespace Marmot::Testing;
using namespace Marmot::FastorStandardTensors;
using namespace Marmot::FastorStandardTensors::Spatial3D;

void testInvertMinorSymmetricFourthOrderTensorMatchesIsotropicCompliance()
{
  // isotropic elasticity tensor C = lambda * I2xI2 + 2*mu*Isym
  const double lambda = 1.0;
  const double mu     = 0.5;

  const Tensor3333d C = lambda * IHyd + 2. * mu * ISymm;

  const Tensor3333d Cinv = Marmot::invertMinorSymmetricFourthOrderTensor( C );

  // analytical isotropic compliance tensor: S = 1/(2*mu) * Isym - lambda / ( 2*mu*(3*lambda + 2*mu) ) * I2xI2
  const double      nuOverE   = lambda / ( 2. * mu * ( 3. * lambda + 2. * mu ) );
  const Tensor3333d Sexpected = 1. / ( 2. * mu ) * ISymm - nuOverE * IHyd;

  throwExceptionOnFailure( checkIfEqual( Cinv, Sexpected, 1e-10 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " invertMinorSymmetricFourthOrderTensor does not match the analytical "
                                           "isotropic compliance tensor" );
}

void testInvertMinorSymmetricFourthOrderTensorIsInvolutive()
{
  // for a different set of isotropic parameters, inverting twice must recover the original tensor
  const double lambda = 2.3;
  const double mu     = 0.7;

  const Tensor3333d C       = lambda * IHyd + 2. * mu * ISymm;
  const Tensor3333d Cinv    = Marmot::invertMinorSymmetricFourthOrderTensor( C );
  const Tensor3333d CinvInv = Marmot::invertMinorSymmetricFourthOrderTensor( Cinv );

  throwExceptionOnFailure( checkIfEqual( CinvInv, C, 1e-8 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " inverting the inverse of a minor-symmetric fourth order tensor must "
                                           "recover the original tensor" );
}

void testIdentity3D()
{
  using namespace Marmot;
  const auto I  = identity3D< double >();
  const auto Ic = identity3D< std::complex< double > >();
  for ( int i = 0; i < 3; i++ )
    for ( int j = 0; j < 3; j++ ) {
      throwExceptionOnFailure( I( i, j ) == ( i == j ? 1.0 : 0.0 ) &&
                                 I( i, j ) == FastorStandardTensors::Spatial3D::I( i, j ),
                               "identity3D< double >" );
      throwExceptionOnFailure( Ic( i, j ) == std::complex< double >( i == j ? 1.0 : 0.0, 0.0 ),
                               "identity3D< complex >" );
    }
}

void testDeviatoricTransposeIsTransposeOfDeviatoric()
{
  // DeviatoricTranspose used to self-referentially read its own (uninitialized) value in its
  // initializer instead of transposing Deviatoric; guard against a regression to that bug.
  const Tensor3333d expected = Fastor::transpose( Deviatoric );

  throwExceptionOnFailure( checkIfEqual( DeviatoricTranspose, expected ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " DeviatoricTranspose must equal transpose(Deviatoric)" );

  const Tensor3333d zero( 0.0 );
  throwExceptionOnFailure( !checkIfEqual( DeviatoricTranspose, zero, 1e-12 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " DeviatoricTranspose must not degenerate to the all-zero tensor "
                                           "produced by the old self-referential initialization bug" );
}

int main()
{
  auto
    tests = std::vector< std::function< void() > >{ testInvertMinorSymmetricFourthOrderTensorMatchesIsotropicCompliance,
                                                    testInvertMinorSymmetricFourthOrderTensorIsInvolutive,
                                                    testDeviatoricTransposeIsTransposeOfDeviatoric,
                                                    testIdentity3D };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
