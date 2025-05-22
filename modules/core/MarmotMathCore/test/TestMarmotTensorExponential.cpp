#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotTensorExponential.h"
#include "Marmot/MarmotTesting.h"
#include <Fastor/Fastor.h>
#include <cmath>

using namespace Marmot::Testing;

using namespace Marmot;

void test_tensorExponential()
{

  Fastor::Tensor< double, 2, 2 > tensor2d = { { .1, 0. }, { 0., .2 } };

  const auto
    tExp2D = ContinuumMechanics::TensorUtility::TensorExponential::computeTensorExponential< double, 2 >( tensor2d,
                                                                                                          2000,
                                                                                                          1e-8,
                                                                                                          1e-8 );

  throwExceptionOnFailure( checkIfEqual( ( tExp2D( 0, 0 ) ), std::exp( 0.1 ), 1e-8 ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
  throwExceptionOnFailure( checkIfEqual( ( tExp2D( 1, 1 ) ), std::exp( 0.2 ), 1e-8 ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );

  const auto [tExp2D_2, tExp2D_2_d_dExp2D] = ContinuumMechanics::TensorUtility::TensorExponential::FirstOrderDerived::
    computeTensorExponential( tensor2d, 1000, 1e-10 );

  throwExceptionOnFailure( checkIfEqual( ( tExp2D_2_d_dExp2D( 0, 0, 0, 0 ) ), std::exp( 0.1 ), 1e-8 ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
  throwExceptionOnFailure( checkIfEqual( ( tExp2D_2_d_dExp2D( 1, 1, 1, 1 ) ), std::exp( 0.2 ), 1e-8 ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

int main()
{

  test_tensorExponential();

  return 0;
}
