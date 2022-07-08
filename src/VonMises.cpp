#include "Marmot/HaighWestergaard.h"
#include "Marmot/MarmotConstants.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotTensor.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"
#include "Marmot/VonMises.h"
#include <iostream>
#include <map>

namespace Marmot::Materials {

  using namespace Eigen;
  using namespace Marmot;

  void VonMisesModel::assignStateVars( double* stateVars, int nStateVars )
  {
    if ( nStateVars < getNumberOfRequiredStateVars() )
      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": Not sufficient stateVars!" );

    managedStateVars = std::make_unique< VonMisesModelStateVarManager >( stateVars );
  }

  StateView VonMisesModel::getStateView( const std::string& stateName )
  {
    return managedStateVars->getStateView( stateName );
  }

  void VonMisesModel::computeStress( double*       stress,
                                     double*       dStress_dStrain,
                                     const double* dStrain,
                                     const double* timeOld,
                                     const double  dT,
                                     double&       pNewDT )

  {
    mVector6d  S( stress );
    mMatrix6d  dS_dE( dStress_dStrain );
    const auto dE = Map< const Vector6d >( dStrain );

    int i = 0;

    // elasticity parameters
    const double& E  = this->materialProperties[i++];
    const double& nu = this->materialProperties[i++];
    // plasticity parameters
    const double& yieldStress      = this->materialProperties[i++];
    const double& HLin             = this->materialProperties[i++];
    const double& deltaYieldStress = this->materialProperties[i++];
    const double& delta            = this->materialProperties[i++];

    const auto Cel = ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( E, nu );

    double& kappa = managedStateVars->kappa;

    if ( dE.isZero( 1e-14 ) ) {
      dS_dE = Cel;
      return;
    }

    // isotropic hardening law
    auto fy = [&]( double kappa_ ) {
      return yieldStress + HLin * kappa_ + deltaYieldStress * ( 1. - std::exp( -delta * kappa_ ) );
    };

    // yield function
    auto f = [&]( double rho_, double kappa_ ) { return rho_ - Constants::sqrt2_3 * fy( kappa_ ); };

    Vector6d trialStress = S + Cel * dE;

    std::cout << trialStress << std::endl;

    const auto hwTrial = ContinuumMechanics::HaighWestergaard::haighWestergaard( trialStress );

    if ( f( hwTrial.rho, kappa ) >= 0. ) {
      // plastic stepi
      const double G = E / ( 2. * ( 1. - nu ) );

      auto g = [&]( double deltaKappa ) {
        return hwTrial.rho - Constants::sqrt6 * G * deltaKappa - Constants::sqrt2_3 * fy( kappa + deltaKappa );
      };

      double       dKappa  = 0.0;
      const double tol     = 1e-12;
      int          counter = 0;
      double       dg_ddKappa;

      while ( std::abs( g( dKappa ) ) > tol ) {
        dg_ddKappa = -Constants::sqrt6 * G -
                     Constants::sqrt2_3 * ( HLin + deltaYieldStress * delta * std::exp( -delta * ( kappa + dKappa ) ) );

        // std::cout << "dg_dkappa =" << dg_ddKappa << std::endl;
        // std::cout << "g = " << g( dKappa ) << std::endl;
        // std::cout << dKappa << std::endl;
        // std::cout << counter << std::endl;

        dKappa -= g( dKappa ) / dg_ddKappa;
        counter += 1;

        if ( counter == 5 ) {
          pNewDT = 0.25;
          return;
        }
      }

      double dLambda = Constants::sqrt3_2 * dKappa;
      Vector6d n = ContinuumMechanics::VoigtNotation::IDev * trialStress / hwTrial.rho;
      Vector6d dEp = dLambda * n;
      std::cout << "dEp = " << dEp << std::endl;
      std::cout << " n = " << n << std::endl;
      S            = trialStress - Cel * dEp;

      const auto hw = ContinuumMechanics::HaighWestergaard::haighWestergaard( Vector6d( S ) );
      std::cout << hw.rho << std::endl;

      std::cout << "f =" << f( hw.rho, kappa + dKappa ) << std::endl;
      dS_dE = Cel;
      kappa += dKappa;
    }
    else {
      // elastic step
      S     = trialStress;
      dS_dE = Cel;
    }

    std::cout << S << std::endl;
  }

} // namespace Marmot::Materials
