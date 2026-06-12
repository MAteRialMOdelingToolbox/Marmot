#include "Marmot/KelvinChainInterface.h"
#include "Marmot/LinearViscoelasticPowerLaw.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotInterfaceMaterialHelperFunctions.h"
#include "Marmot/MarmotInterfaceMaterialHypoElastic.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotVoigt.h"

#include <Eigen/Dense>

#include <functional>
#include <memory>
#include <stdexcept>
#include <vector>

using namespace Marmot::Testing;

namespace {

  struct InterfaceResponse {
    Eigen::Vector3d force = Eigen::Vector3d::Zero();
    Eigen::Matrix< double, 3, 3, Eigen::RowMajor >
      surfaceStress = Eigen::Matrix< double, 3, 3, Eigen::RowMajor >::Zero();
  };

  std::unique_ptr< MarmotInterfaceMaterialHypoElastic > createInterfaceMaterial( const double* properties,
                                                                                 int           nProperties )
  {
    auto material = std::unique_ptr< MarmotInterfaceMaterialHypoElastic >(
      MarmotLibrary::MarmotInterfaceMaterialHypoElasticFactory::createMaterial( "KELVINCHAININTERFACE",
                                                                                properties,
                                                                                nProperties,
                                                                                1 ) );
    if ( !material ) {
      throw std::runtime_error( "KelvinChainInterface registration failed." );
    }
    return material;
  }

  void testAgainstBulkKelvinChain()
  {
    // Interface properties: [E, nu, h, m, n, nKelvin, minTau, timeToDays]
    const double interfaceProperties[8] = { 2e5, 0.2, 0.01, 0.5, 0.1, 10., 0.0001, 1. };
    // Wrapped LinearViscoelasticPowerLaw properties, without h.
    const double bulkProperties[7] = { 2e5, 0.2, 0.5, 0.1, 10., 0.0001, 1. };
    const double h                 = interfaceProperties[2];

    auto                                          interfaceMaterial = createInterfaceMaterial( interfaceProperties, 8 );
    Marmot::Materials::LinearViscoelasticPowerLaw bulkMaterial( bulkProperties, 7, 1 );

    Eigen::VectorXd interfaceStateVars( interfaceMaterial->getNumberOfRequiredStateVars() );
    Eigen::VectorXd bulkStateVars( bulkMaterial.getNumberOfRequiredStateVars() );
    interfaceMaterial->initializeYourself( interfaceStateVars.data(), interfaceStateVars.size() );
    bulkMaterial.initializeYourself( bulkStateVars.data(), bulkStateVars.size() );

    InterfaceResponse interfaceResponse;
    Marmot::Vector6d  bulkStress = Marmot::Vector6d::Zero();

    const double normal[3] = { 0., 0., 1. };

    struct Increment {
      double dT;
      double jumpY;
      double surfaceShear;
    };
    const std::vector< Increment > increments = {
      { 0.01, 1e-4, 2e-4 },
      { 10.0, 0.0, 0.0 },
    };

    double timeOld = 0.0;
    for ( const auto& increment : increments ) {
      const double dU[6]              = { 0., increment.jumpY, 0., 0., 0., 0. };
      const double dSurfaceStrain[18] = { 0.,
                                          increment.surfaceShear,
                                          0.,
                                          increment.surfaceShear,
                                          0.,
                                          0.,
                                          0.,
                                          0.,
                                          0.,
                                          0.,
                                          increment.surfaceShear,
                                          0.,
                                          increment.surfaceShear,
                                          0.,
                                          0.,
                                          0.,
                                          0.,
                                          0. };

      double Q[9]  = { 0. };
      double Z[81] = { 0. };
      double H[27] = { 0. };
      double Y[81] = { 0. };

      MarmotInterfaceMaterialHypoElastic::State         state{ interfaceResponse.force.data(),
                                                       interfaceResponse.surfaceStress.data(),
                                                       interfaceStateVars.data() };
      MarmotInterfaceMaterialHypoElastic::Tangents      tangents{ Q, Z, H, Y };
      MarmotInterfaceMaterialHypoElastic::Deformation   deformation{ dU, dSurfaceStrain, normal };
      MarmotInterfaceMaterialHypoElastic::TimeIncrement timeIncrement{ timeOld, increment.dT };
      interfaceMaterial->computeStress( state, tangents, deformation, timeIncrement );

      Marmot::Vector6d bulkStrainIncrement = Marmot::Vector6d::Zero();
      bulkStrainIncrement[3]               = 2. * increment.surfaceShear;
      bulkStrainIncrement[5]               = increment.jumpY / h;

      Marmot::Matrix6d                    bulkTangent = Marmot::Matrix6d::Zero();
      MarmotMaterialHypoElastic::state3D  bulkState{ bulkStress, 0.0, 0.0, bulkStateVars.data() };
      MarmotMaterialHypoElastic::timeInfo bulkTime{ timeOld, increment.dT };
      bulkMaterial.computeStress( bulkState, bulkTangent, bulkStrainIncrement, bulkTime );
      bulkStress = bulkState.stress;

      const Eigen::Matrix3d expectedStress = Marmot::ContinuumMechanics::VoigtNotation::voigtToStress( bulkStress );
      const Eigen::Vector3d expectedForce  = expectedStress * Eigen::Vector3d::UnitZ();

      throwExceptionOnFailure( checkIfEqual< double >( interfaceResponse.force, expectedForce, 1e-10 ),
                               "Interface force does not match wrapped Kelvin-chain stress." );
      throwExceptionOnFailure( checkIfEqual< double >( interfaceResponse.surfaceStress, h * expectedStress, 1e-10 ),
                               "Interface surface stress does not match wrapped Kelvin-chain stress." );
      throwExceptionOnFailure( checkIfEqual< double >( interfaceStateVars, bulkStateVars, 1e-10 ),
                               "Interface and wrapped Kelvin-chain state variables differ." );

      const Marmot::FastorStandardTensors::Tensor3d normalTensor( normal );
      auto [expectedZ, expectedQ, expectedH, expectedY] = Marmot::Materials::InterfaceMaterialHelperFunctions::
        calculateInterfaceMaterialParameters( normalTensor, bulkTangent );

      const Eigen::Map< const Eigen::VectorXd > actualQ( Q, 9 );
      const Eigen::Map< const Eigen::VectorXd > actualZ( Z, 81 );
      const Eigen::Map< const Eigen::VectorXd > actualH( H, 27 );
      const Eigen::Map< const Eigen::VectorXd > actualY( Y, 81 );
      const Eigen::Map< const Eigen::VectorXd > expectedQVector( expectedQ.data(), 9 );
      const Eigen::Map< const Eigen::VectorXd > expectedZVector( expectedZ.data(), 81 );
      const Eigen::Map< const Eigen::VectorXd > expectedHVector( expectedH.data(), 27 );
      const Eigen::Map< const Eigen::VectorXd > expectedYVector( expectedY.data(), 81 );

      throwExceptionOnFailure( checkIfEqual< double >( actualQ, ( 1. / h ) * expectedQVector, 1e-10 ),
                               "Interface Q tangent does not match the wrapped Kelvin-chain tangent." );
      throwExceptionOnFailure( checkIfEqual< double >( actualZ, h * expectedZVector, 1e-10 ),
                               "Interface Z tangent does not match the wrapped Kelvin-chain tangent." );
      throwExceptionOnFailure( checkIfEqual< double >( actualH, expectedHVector, 1e-10 ),
                               "Interface H tangent does not match the wrapped Kelvin-chain tangent." );
      throwExceptionOnFailure( checkIfEqual< double >( actualY, h * expectedYVector, 1e-10 ),
                               "Interface Y tangent does not match the wrapped Kelvin-chain tangent." );

      timeOld += increment.dT;
    }
  }

} // namespace

int main()
{
  executeTestsAndCollectExceptions( { testAgainstBulkKelvinChain } );
  return 0;
}
