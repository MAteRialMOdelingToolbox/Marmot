#include "Marmot/LinearViscoElasticWiechert.h"
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotTesting.h"

#include <Eigen/Dense>

#include <functional>
#include <vector>

using namespace Marmot::Testing;

namespace {

  Eigen::Vector< double, 7 > getMaterialPropertiesLinearViscoElasticWiechert()
  {
    Eigen::Vector< double, 7 > materialProperties;
    materialProperties << 2e5, 0.2, 0.5, 0.1, 10., 0.0001, 1.;
    return materialProperties;
  }

  void testLinearViscoElasticWiechertReferenceStress()
  {
    auto materialProperties = getMaterialPropertiesLinearViscoElasticWiechert();
    Marmot::Materials::LinearViscoElasticWiechert material( materialProperties.data(), materialProperties.size(), 1 );

    Eigen::VectorXd stateVars( material.getNumberOfRequiredStateVars() );
    material.initializeYourself( stateVars.data(), stateVars.size() );

    MarmotMaterialHypoElastic::state3D state{ Marmot::Vector6d::Zero(), 0.0, 0.0, stateVars.data() };
    Marmot::Matrix6d                   tangent = Marmot::Matrix6d::Zero();

    Marmot::Vector6d dStrain = Marmot::Vector6d::Zero();
    material.computeStress( state, tangent, dStrain, { 0.0, 28.0 } );

    dStrain << 10e-6, 10e-6, 0., 0., 0., 0.;
    material.computeStress( state, tangent, dStrain, { 28.0, 0.01 } );

    dStrain.setZero();
    material.computeStress( state, tangent, dStrain, { 28.01, 100.0 } );

    Marmot::Vector6d stressTarget;
    stressTarget << 2.77778377468268234, 2.77778377468268234, 1.11111350987307289, 0., 0., 0.;

    throwExceptionOnFailure( checkIfEqual< double >( state.stress, stressTarget, 1e-10 ),
                             "LinearViscoElasticWiechert reference stress mismatch." );
  }

  void testDensity()
  {
    const double                                  properties[8] = { 1e8, 0.3, 2e7, 0.25, 4., 1e-4, 1., 2400. };
    Marmot::Materials::LinearViscoElasticWiechert material( properties, 8, 1 );
    throwExceptionOnFailure( checkIfEqual( material.getDensity( nullptr ), properties[7] ),
                             "LinearViscoElasticWiechert density is incorrect." );
  }

} // namespace

int main()
{
  std::vector< std::function< void() > > tests = {
    testLinearViscoElasticWiechertReferenceStress,
    testDensity,
  };
  executeTestsAndCollectExceptions( tests );
  return 0;
}
