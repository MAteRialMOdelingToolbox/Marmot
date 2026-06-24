#include "Marmot/MarmotInterfaceMaterialHypoElastic.h"
#include "Marmot/MarmotTesting.h"

#include <Eigen/Dense>

#include <functional>
#include <memory>
#include <stdexcept>
#include <vector>

using namespace Marmot::Testing;

namespace {

  std::unique_ptr< MarmotInterfaceMaterialHypoElastic > createMaterial( const double* props, int nProps )
  {
    const int elLabel = 1;
    auto      mat     = std::unique_ptr< MarmotInterfaceMaterialHypoElastic >(
      MarmotLibrary::MarmotInterfaceMaterialHypoElasticFactory::createMaterial( "WIECHERTINTERFACE",
                                                                                props,
                                                                                nProps,
                                                                                elLabel ) );
    if ( !mat )
      throw std::runtime_error( "WiechertInterfaceMaterial creation failed." );
    return mat;
  }

  void computeStress( MarmotInterfaceMaterialHypoElastic& mat,
                      double*                             stateVars,
                      double*                             force,
                      double*                             surfaceStress,
                      double*                             Q_ij,
                      double*                             Z_ijkl,
                      double*                             H_ijk,
                      double*                             Y_ijkl,
                      const double*                       dU,
                      const double*                       dSurfaceStrain,
                      const double*                       normal,
                      const double                        timeOld,
                      const double                        dT )
  {
    MarmotInterfaceMaterialHypoElastic::State         state{ force, surfaceStress, stateVars };
    MarmotInterfaceMaterialHypoElastic::Tangents      tangents{ Q_ij, Z_ijkl, H_ijk, Y_ijkl };
    MarmotInterfaceMaterialHypoElastic::Deformation   deformation{ dU, dSurfaceStrain, normal };
    MarmotInterfaceMaterialHypoElastic::TimeIncrement timeIncrement{ timeOld, dT };
    mat.computeStress( state, tangents, deformation, timeIncrement );
  }

  void runSingleIncrement( const double* props,
                           const int     nProps,
                           const double* dU,
                           const double* dSurfaceStrain,
                           const double* normal,
                           const double  dT,
                           double*       force,
                           double*       surfaceStress )
  {
    auto mat = createMaterial( props, nProps );

    Eigen::VectorXd stateVars( mat->getNumberOfRequiredStateVars() );
    mat->initializeYourself( stateVars.data(), stateVars.size() );

    double       Q_ij[9]    = { 0. };
    double       Z_ijkl[81] = { 0. };
    double       H_ijk[27]  = { 0. };
    double       Y_ijkl[81] = { 0. };
    const double timeOld    = 0.0;

    computeStress( *mat,
                   stateVars.data(),
                   force,
                   surfaceStress,
                   Q_ij,
                   Z_ijkl,
                   H_ijk,
                   Y_ijkl,
                   dU,
                   dSurfaceStrain,
                   normal,
                   timeOld,
                   dT );
  }

  void testDisplacementJump()
  {
    const double props[8]           = { 1e8, 0.3, 0.01, 2e7, 0.25, 6., 1e-4, 1. };
    const double dU[6]              = { 0., 1e-4, 0., 0., 0., 0. };
    const double dSurfaceStrain[18] = { 0. };
    const double normal[3]          = { 0., 0., 1. };

    double force[3]         = { 0. };
    double surfaceStress[9] = { 0. };

    runSingleIncrement( props, 8, dU, dSurfaceStrain, normal, 1e-6, force, surfaceStress );

    const double forceTarget[3]         = { 0., 1.31415499891368486e6, 0. };
    const double surfaceStressTarget[9] = { 0., 0., 0., 0., 0., 1.31415499891368490e4, 0., 1.31415499891368490e4, 0. };

    Eigen::Map< const Eigen::Vector3d > forceVec( force );
    Eigen::Map< const Eigen::Vector3d > forceTgt( forceTarget );
    Eigen::Map< const Eigen::VectorXd > surfaceStressVec( surfaceStress, 9 );
    Eigen::Map< const Eigen::VectorXd > surfaceStressTgt( surfaceStressTarget, 9 );

    throwExceptionOnFailure( checkIfEqual< double >( forceVec, forceTgt, 1e-8 ),
                             "force mismatch in " + std::string( __PRETTY_FUNCTION__ ) );
    throwExceptionOnFailure( checkIfEqual< double >( surfaceStressVec, surfaceStressTgt, 1e-8 ),
                             "surface stress mismatch in " + std::string( __PRETTY_FUNCTION__ ) );
  }

  void testSurfaceStrain()
  {
    const double props[8] = { 1e8, 0.3, 0.01, 2e7, 0.25, 6., 1e-4, 1. };
    const double dU[6]    = { 0. };
    const double dSurfaceStrain[18] =
      { 0., 2e-4, 0., 2e-4, 0., 0., 0., 0., 0., 0., 2e-4, 0., 2e-4, 0., 0., 0., 0., 0. };
    const double normal[3] = { 0., 0., 1. };

    double force[3]         = { 0. };
    double surfaceStress[9] = { 0. };

    runSingleIncrement( props, 8, dU, dSurfaceStrain, normal, 1e-6, force, surfaceStress );

    const double forceTarget[3]         = { 0., 0., 0. };
    const double surfaceStressTarget[9] = { 0., 5.25661999565474048e2, 0., 5.25661999565474048e2, 0., 0., 0., 0., 0. };

    Eigen::Map< const Eigen::Vector3d > forceVec( force );
    Eigen::Map< const Eigen::Vector3d > forceTgt( forceTarget );
    Eigen::Map< const Eigen::VectorXd > surfaceStressVec( surfaceStress, 9 );
    Eigen::Map< const Eigen::VectorXd > surfaceStressTgt( surfaceStressTarget, 9 );

    throwExceptionOnFailure( checkIfEqual< double >( forceVec, forceTgt, 1e-8 ),
                             "force mismatch in " + std::string( __PRETTY_FUNCTION__ ) );
    throwExceptionOnFailure( checkIfEqual< double >( surfaceStressVec, surfaceStressTgt, 1e-8 ),
                             "surface stress mismatch in " + std::string( __PRETTY_FUNCTION__ ) );
  }

  void testDensityDelegation()
  {
    const double props[9] = { 1e8, 0.3, 0.01, 2e7, 0.25, 6., 1e-4, 1., 2400. };
    auto         mat      = createMaterial( props, 9 );
    throwExceptionOnFailure( checkIfEqual( mat->getDensity(), props[8] ),
                             "WiechertInterfaceMaterial density delegation failed." );
  }

} // namespace

int main()
{
  std::vector< std::function< void() > > tests = {
    testDisplacementJump,
    testSurfaceStrain,
    testDensityDelegation,
  };
  executeTestsAndCollectExceptions( tests );
  return 0;
}
