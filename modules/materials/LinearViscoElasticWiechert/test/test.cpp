#include "Marmot/LinearViscoElasticWiechert.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotTesting.h"

#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <functional>
#include <stdexcept>
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

  template < typename Callable >
  bool throwsInvalidArgument( Callable&& call )
  {
    try {
      call();
    }
    catch ( const std::invalid_argument& ) {
      return true;
    }
    return false;
  }

  void testConstructorRejectsInvalidProperties()
  {
    using Marmot::Materials::LinearViscoElasticWiechert;
    // E, nu, m, n, nMaxwell, minTau, timeToDays
    const double valid[7] = { 2e5, 0.2, 0.5, 0.1, 10., 1e-4, 1. };

    throwExceptionOnFailure( throwsInvalidArgument( [&] { LinearViscoElasticWiechert material( valid, 6, 1 ); } ),
                             "Fewer than seven properties must be rejected." );
    throwExceptionOnFailure( throwsInvalidArgument( [&] { LinearViscoElasticWiechert material( nullptr, 7, 1 ); } ),
                             "A missing property array must be rejected." );

    const auto withProperty = [&]( int index, double value ) {
      std::array< double, 7 > properties;
      std::copy( valid, valid + 7, properties.begin() );
      properties[index] = value;
      return properties;
    };

    const auto noMaxwellElement = withProperty( 4, 0. );
    throwExceptionOnFailure( throwsInvalidArgument(
                               [&] { LinearViscoElasticWiechert material( noMaxwellElement.data(), 7, 1 ); } ),
                             "Fewer than one Maxwell element must be rejected." );

    const auto nonPositiveMinTau = withProperty( 5, 0. );
    throwExceptionOnFailure( throwsInvalidArgument(
                               [&] { LinearViscoElasticWiechert material( nonPositiveMinTau.data(), 7, 1 ); } ),
                             "A non-positive minimum relaxation time must be rejected." );

    const auto negativeM = withProperty( 2, -1. );
    throwExceptionOnFailure( throwsInvalidArgument(
                               [&] { LinearViscoElasticWiechert material( negativeM.data(), 7, 1 ); } ),
                             "A negative m must be rejected." );

    const auto nonPositiveN = withProperty( 3, 0. );
    throwExceptionOnFailure( throwsInvalidArgument(
                               [&] { LinearViscoElasticWiechert material( nonPositiveN.data(), 7, 1 ); } ),
                             "A non-positive n must be rejected." );
  }

  void testZeroIncrementReturnsInstantaneousElasticTangent()
  {
    auto materialProperties = getMaterialPropertiesLinearViscoElasticWiechert();
    Marmot::Materials::LinearViscoElasticWiechert material( materialProperties.data(), materialProperties.size(), 1 );

    Eigen::VectorXd stateVars( material.getNumberOfRequiredStateVars() );
    material.initializeYourself( stateVars.data(), stateVars.size() );

    MarmotMaterialHypoElastic::state3D state{ Marmot::Vector6d::Zero(), 0.0, 0.0, stateVars.data() };
    Marmot::Matrix6d                   tangent = Marmot::Matrix6d::Zero();

    // neither strain nor time advance: the tangent is the instantaneous isotropic stiffness
    material.computeStress( state, tangent, Marmot::Vector6d::Zero(), { 0.0, 0.0 } );

    const Marmot::Matrix6d
      expected = Marmot::ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( materialProperties( 0 ),
                                                                                     materialProperties( 1 ) );
    throwExceptionOnFailure( checkIfEqual< double >( tangent, expected, 1e-8 * expected.maxCoeff() ),
                             "Zero-increment tangent is not the instantaneous elastic stiffness." );
    throwExceptionOnFailure( state.stress.isZero(), "A zero increment must not change the stress." );
  }

  void testDensityMustBeProvided()
  {
    const double                                  properties[7] = { 2e5, 0.2, 0.5, 0.1, 10., 1e-4, 1. };
    Marmot::Materials::LinearViscoElasticWiechert material( properties, 7, 1 );

    bool rejected = false;
    try {
      material.getDensity( nullptr );
    }
    catch ( const std::runtime_error& ) {
      rejected = true;
    }
    throwExceptionOnFailure( rejected, "A missing density property must be reported." );
  }

} // namespace

int main()
{
  std::vector< std::function< void() > > tests = {
    testLinearViscoElasticWiechertReferenceStress,
    testDensity,
    testConstructorRejectsInvalidProperties,
    testZeroIncrementReturnsInstantaneousElasticTangent,
    testDensityMustBeProvided,
  };
  executeTestsAndCollectExceptions( tests );
  return 0;
}
