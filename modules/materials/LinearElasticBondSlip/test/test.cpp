#include "Marmot/MarmotBondSlipLawFactory.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotTesting.h"
#include <memory>

using namespace Marmot;
using namespace Marmot::Testing;
using namespace Marmot::Math;

void testLinearElasticBondSlip()
{
  const double props[] = { 200., 5000. };
  auto         law     = std::unique_ptr< Marmot::Materials::BondSlipLaw >(
    Marmot::Registry::BondSlipLawFactory::createBondSlipLaw( "LINEARELASTICBONDSLIP", props, 2, 1 ) );

  Marmot::Materials::BondSlipLaw::State state{ Vector3d::Zero(), 0.0, 0.0, nullptr };
  Matrix3d                              C;
  const Vector3d                        slip( 0.1, -0.02, 0.03 );
  law->computeBondStress( state, C, slip, slip, { 0.0, 1.0 } );

  throwExceptionOnFailure( checkIfEqual< double >( state.bondStress, Vector3d( 20., -100., 150. ), 1e-12 ),
                           "linear bond stress" );
  throwExceptionOnFailure( checkIfEqual< double >( C, Matrix3d( Vector3d( 200., 5000., 5000. ).asDiagonal() ), 1e-12 ),
                           "linear bond tangent" );
  throwExceptionOnFailure( checkIfEqual( state.elasticEnergyDensity, 0.5 * slip.dot( state.bondStress ), 1e-12 ),
                           "elastic energy" );
  throwExceptionOnFailure( law->getNumberOfRequiredStateVars() == 0, "no state variables" );
}

void testUnknownBondSlipLawThrows()
{
  const double props[] = { 1., 1. };
  bool         threw   = false;
  try {
    Marmot::Registry::BondSlipLawFactory::createBondSlipLaw( "NOTABONDSLIPLAW", props, 2, 1 );
  }
  catch ( const std::invalid_argument& ) {
    threw = true;
  }
  throwExceptionOnFailure( threw, "an unknown bond-slip law must be rejected" );
}

int main()
{
  executeTestsAndCollectExceptions( { testLinearElasticBondSlip, testUnknownBondSlipLawThrows } );
  return 0;
}
