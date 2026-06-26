#include "Marmot/MarmotInterfaceMaterialHelperFunctions.h"
#include "Marmot/MarmotInterfaceMaterialHypoElastic.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotVoigt.h"

#include <Eigen/Dense>

#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Marmot::Testing;

namespace {

  struct InterfaceResponse {
    Eigen::Vector3d force = Eigen::Vector3d::Zero();
    Eigen::Matrix< double, 3, 3, Eigen::RowMajor >
                    surfaceStress = Eigen::Matrix< double, 3, 3, Eigen::RowMajor >::Zero();
    Eigen::VectorXd Q             = Eigen::VectorXd::Zero( 9 );
    Eigen::VectorXd Z             = Eigen::VectorXd::Zero( 81 );
    Eigen::VectorXd H             = Eigen::VectorXd::Zero( 27 );
    Eigen::VectorXd Y             = Eigen::VectorXd::Zero( 81 );
  };

  std::unique_ptr< MarmotInterfaceMaterialHypoElastic > createInterfaceMaterial( const std::string& materialName,
                                                                                 const double*      properties,
                                                                                 int                nProperties )
  {
    auto material = std::unique_ptr< MarmotInterfaceMaterialHypoElastic >(
      MarmotLibrary::MarmotInterfaceMaterialHypoElasticFactory::createMaterial( materialName,
                                                                                properties,
                                                                                nProperties,
                                                                                1 ) );
    if ( !material ) {
      throw std::runtime_error( "MarmotInterfaceMaterialHypoElastic creation failed." );
    }
    return material;
  }

  std::unique_ptr< MarmotMaterialHypoElastic > createBulkMaterial( const std::string& materialName,
                                                                   const double*      properties,
                                                                   int                nProperties )
  {
    return std::unique_ptr< MarmotMaterialHypoElastic >(
      MarmotLibrary::MarmotMaterialHypoElasticFactory::createMaterial( materialName, properties, nProperties, 1 ) );
  }

  void computeInterfaceStress( MarmotInterfaceMaterialHypoElastic& interfaceMaterial,
                               InterfaceResponse&                  response,
                               double*                             stateVars,
                               const double*                       dU,
                               const double*                       dSurfaceStrain,
                               const double*                       normal,
                               double                              timeOld,
                               double                              dT )
  {
    double Q[9]  = { 0. };
    double Z[81] = { 0. };
    double H[27] = { 0. };
    double Y[81] = { 0. };

    MarmotInterfaceMaterialHypoElastic::State state{ response.force.data(), response.surfaceStress.data(), stateVars };
    MarmotInterfaceMaterialHypoElastic::Tangents      tangents{ Q, Z, H, Y };
    MarmotInterfaceMaterialHypoElastic::Deformation   deformation{ dU, dSurfaceStrain, normal };
    MarmotInterfaceMaterialHypoElastic::TimeIncrement timeIncrement{ timeOld, dT };

    interfaceMaterial.computeStress( state, tangents, deformation, timeIncrement );

    response.Q = Eigen::Map< const Eigen::VectorXd >( Q, 9 );
    response.Z = Eigen::Map< const Eigen::VectorXd >( Z, 81 );
    response.H = Eigen::Map< const Eigen::VectorXd >( H, 27 );
    response.Y = Eigen::Map< const Eigen::VectorXd >( Y, 81 );
  }

  void testGenericInterfaceAgainstBulkMaterial( const std::string& materialName,
                                                const double*      interfaceProperties,
                                                int                nInterfaceProperties,
                                                const double*      bulkProperties,
                                                int                nBulkProperties )
  {
    constexpr double h         = 0.01;
    const double     normal[3] = { 0., 0., 1. };

    auto interfaceMaterial = createInterfaceMaterial( materialName, interfaceProperties, nInterfaceProperties );
    auto bulkMaterial      = createBulkMaterial( materialName, bulkProperties, nBulkProperties );

    Eigen::VectorXd interfaceStateVars( interfaceMaterial->getNumberOfRequiredStateVars() );
    Eigen::VectorXd bulkStateVars( bulkMaterial->getNumberOfRequiredStateVars() );
    interfaceMaterial->initializeYourself( interfaceStateVars.data(), interfaceStateVars.size() );
    bulkMaterial->initializeYourself( bulkStateVars.data(), bulkStateVars.size() );

    InterfaceResponse interfaceResponse;
    Marmot::Vector6d  bulkStress = Marmot::Vector6d::Zero();

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

      computeInterfaceStress( *interfaceMaterial,
                              interfaceResponse,
                              interfaceStateVars.data(),
                              dU,
                              dSurfaceStrain,
                              normal,
                              timeOld,
                              increment.dT );

      Marmot::Vector6d bulkStrainIncrement = Marmot::Vector6d::Zero();
      bulkStrainIncrement[3]               = 2. * increment.surfaceShear;
      bulkStrainIncrement[5]               = increment.jumpY / h;

      Marmot::Matrix6d                   bulkTangent = Marmot::Matrix6d::Zero();
      MarmotMaterialHypoElastic::state3D bulkState{ bulkStress, 0.0, 0.0, bulkStateVars.data() };
      bulkMaterial->computeStress( bulkState, bulkTangent, bulkStrainIncrement, { timeOld, increment.dT } );
      bulkStress = bulkState.stress;

      const Eigen::Matrix3d expectedStress = Marmot::ContinuumMechanics::VoigtNotation::voigtToStress( bulkStress );
      const Eigen::Vector3d expectedForce  = expectedStress * Eigen::Vector3d::UnitZ();

      throwExceptionOnFailure( checkIfEqual< double >( interfaceResponse.force, expectedForce, 1e-10 ),
                               materialName + ": interface force does not match bulk stress." );
      throwExceptionOnFailure( checkIfEqual< double >( interfaceResponse.surfaceStress, h * expectedStress, 1e-10 ),
                               materialName + ": interface surface stress does not match bulk stress." );
      throwExceptionOnFailure( checkIfEqual< double >( interfaceStateVars, bulkStateVars, 1e-10 ),
                               materialName + ": interface state variables do not match bulk state variables." );

      const Marmot::FastorStandardTensors::Tensor3d normalTensor( normal );
      auto [expectedZ, expectedQ, expectedH, expectedY] = Marmot::Materials::InterfaceMaterialHelperFunctions::
        calculateInterfaceMaterialParameters( normalTensor, bulkTangent );

      throwExceptionOnFailure( checkIfEqual< double >( interfaceResponse.Q,
                                                       ( 1. / h ) *
                                                         Eigen::Map< const Eigen::VectorXd >( expectedQ.data(), 9 ),
                                                       1e-10 ),
                               materialName + ": Q tangent mismatch." );
      throwExceptionOnFailure( checkIfEqual< double >( interfaceResponse.Z,
                                                       h * Eigen::Map< const Eigen::VectorXd >( expectedZ.data(), 81 ),
                                                       1e-10 ),
                               materialName + ": Z tangent mismatch." );
      throwExceptionOnFailure( checkIfEqual< double >( interfaceResponse.H,
                                                       Eigen::Map< const Eigen::VectorXd >( expectedH.data(), 27 ),
                                                       1e-10 ),
                               materialName + ": H tangent mismatch." );
      throwExceptionOnFailure( checkIfEqual< double >( interfaceResponse.Y,
                                                       h * Eigen::Map< const Eigen::VectorXd >( expectedY.data(), 81 ),
                                                       1e-10 ),
                               materialName + ": Y tangent mismatch." );

      timeOld += increment.dT;
    }
  }

  void testGenericVonMisesInterface()
  {
    const double interfaceProperties[8] = { 1e5, 0.3, 0.01, 100., 10., 0., 1., 2400. };
    const double bulkProperties[7]      = { 1e5, 0.3, 100., 10., 0., 1., 2400. };
    testGenericInterfaceAgainstBulkMaterial( "VONMISES", interfaceProperties, 8, bulkProperties, 7 );
  }

  void testGenericKelvinChainInterface()
  {
    const double interfaceProperties[8] = { 2e5, 0.2, 0.01, 0.5, 0.1, 10., 0.0001, 1. };
    const double bulkProperties[7]      = { 2e5, 0.2, 0.5, 0.1, 10., 0.0001, 1. };
    testGenericInterfaceAgainstBulkMaterial( "LINEARVISCOELASTICPOWERLAW", interfaceProperties, 8, bulkProperties, 7 );
  }

  void testGenericWiechertInterface()
  {
    const double interfaceProperties[9] = { 1e8, 0.3, 0.01, 2e7, 0.25, 6., 1e-4, 1., 2400. };
    const double bulkProperties[8]      = { 1e8, 0.3, 2e7, 0.25, 6., 1e-4, 1., 2400. };
    testGenericInterfaceAgainstBulkMaterial( "LINEARVISCOELASTICWIECHERT", interfaceProperties, 9, bulkProperties, 8 );

    auto interfaceMaterial = createInterfaceMaterial( "LINEARVISCOELASTICWIECHERT", interfaceProperties, 9 );
    throwExceptionOnFailure( checkIfEqual( interfaceMaterial->getDensity(), interfaceProperties[8] ),
                             "Generic Wiechert interface density delegation failed." );
  }

  void testLegacyInterfaceNamesResolveToGenericAdapter()
  {
    const double vonMisesProperties[8] = { 1e5, 0.3, 0.01, 100., 10., 0., 1., 2400. };
    const double kelvinProperties[8]   = { 2e5, 0.2, 0.01, 0.5, 0.1, 10., 0.0001, 1. };
    const double wiechertProperties[9] = { 1e8, 0.3, 0.01, 2e7, 0.25, 6., 1e-4, 1., 2400. };

    auto vonMisesInterface = createInterfaceMaterial( "VONMISESINTERFACE", vonMisesProperties, 8 );
    auto kelvinInterface   = createInterfaceMaterial( "KELVINCHAININTERFACE", kelvinProperties, 8 );
    auto wiechertInterface = createInterfaceMaterial( "WIECHERTINTERFACE", wiechertProperties, 9 );

    throwExceptionOnFailure( vonMisesInterface->getNumberOfRequiredStateVars() > 0,
                             "VONMISESINTERFACE legacy alias did not create a usable material." );
    throwExceptionOnFailure( kelvinInterface->getNumberOfRequiredStateVars() > 0,
                             "KELVINCHAININTERFACE legacy alias did not create a usable material." );
    throwExceptionOnFailure( wiechertInterface->getNumberOfRequiredStateVars() > 0,
                             "WIECHERTINTERFACE legacy alias did not create a usable material." );
  }

} // namespace

int main()
{
  std::vector< std::function< void() > > tests = {
    testGenericVonMisesInterface,
    testGenericKelvinChainInterface,
    testGenericWiechertInterface,
    testLegacyInterfaceNamesResolveToGenericAdapter,
  };
  executeTestsAndCollectExceptions( tests );
  return 0;
}
