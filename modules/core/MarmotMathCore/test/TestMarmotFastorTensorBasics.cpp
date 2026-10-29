#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotTesting.h"

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

void testMapEigenToFastorThirdRankTensor()
{
  // non-cubic sizes, so that a wrong index convention cannot pass by coincidence
  Fastor::Tensor< double, 2, 3, 4 > tensor;
  for ( size_t i = 0; i < 2; ++i )
    for ( size_t j = 0; j < 3; ++j )
      for ( size_t k = 0; k < 4; ++k )
        tensor( i, j, k ) = 100. * i + 10. * j + k;

  const auto matrix = Marmot::mapEigenToFastor( tensor );
  throwExceptionOnFailure( matrix.rows() == 2 && matrix.cols() == 12, "Wrong shape of the third rank tensor map." );

  for ( size_t i = 0; i < 2; ++i )
    for ( size_t j = 0; j < 3; ++j )
      for ( size_t k = 0; k < 4; ++k )
        throwExceptionOnFailure( matrix( i, j * 4 + k ) == tensor( i, j, k ),
                                 MakeString() << __PRETTY_FUNCTION__ << " map(i, j * n3 + k) != tensor(i, j, k)" );

  // the non-const overload writes through to the tensor
  Marmot::mapEigenToFastor( tensor )( 1, 2 * 4 + 3 ) = -1.;
  throwExceptionOnFailure( tensor( 1, 2, 3 ) == -1., "Writing through the third rank tensor map failed." );
}

void testMapEigenToFastorFourthRankTensor()
{
  Fastor::Tensor< double, 2, 3, 4, 5 > tensor;
  for ( size_t i = 0; i < 2; ++i )
    for ( size_t j = 0; j < 3; ++j )
      for ( size_t k = 0; k < 4; ++k )
        for ( size_t l = 0; l < 5; ++l )
          tensor( i, j, k, l ) = 1000. * i + 100. * j + 10. * k + l;

  const auto matrix = Marmot::mapEigenToFastor( tensor );
  throwExceptionOnFailure( matrix.rows() == 6 && matrix.cols() == 20, "Wrong shape of the fourth rank tensor map." );

  for ( size_t i = 0; i < 2; ++i )
    for ( size_t j = 0; j < 3; ++j )
      for ( size_t k = 0; k < 4; ++k )
        for ( size_t l = 0; l < 5; ++l )
          throwExceptionOnFailure( matrix( i * 3 + j, k * 5 + l ) == tensor( i, j, k, l ),
                                   MakeString()
                                     << __PRETTY_FUNCTION__ << " map(i * n2 + j, k * n4 + l) != tensor(i, j, k, l)" );

  // the non-const overload writes through to the tensor
  Marmot::mapEigenToFastor( tensor )( 1 * 3 + 2, 3 * 5 + 4 ) = -1.;
  throwExceptionOnFailure( tensor( 1, 2, 3, 4 ) == -1., "Writing through the fourth rank tensor map failed." );
}

int main()
{
  auto
    tests = std::vector< std::function< void() > >{ testInvertMinorSymmetricFourthOrderTensorMatchesIsotropicCompliance,
                                                    testInvertMinorSymmetricFourthOrderTensorIsInvolutive,
                                                    testDeviatoricTransposeIsTransposeOfDeviatoric,
                                                    testMapEigenToFastorThirdRankTensor,
                                                    testMapEigenToFastorFourthRankTensor };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
