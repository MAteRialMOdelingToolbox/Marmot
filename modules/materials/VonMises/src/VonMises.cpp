#include "Marmot/VonMises.h"
#include "Marmot/HaighWestergaard.h"
#include "Marmot/MarmotConstants.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotTensor.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"
#include "Marmot/VonMisesConstants.h"
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
    return MarmotMaterialHypoElastic::assignStateVars( stateVars, nStateVars );
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
    // elasticity parameters
    const double& E  = this->materialProperties[0];
    const double& nu = this->materialProperties[1];
    // plasticity parameters
    const double& yieldStress      = this->materialProperties[2];
    const double& HLin             = this->materialProperties[3];
    const double& deltaYieldStress = this->materialProperties[4];
    const double& delta            = this->materialProperties[5];

    // map to stress, strain and tangent
    mVector6d  S( stress );
    mMatrix6d  dS_dE( dStress_dStrain );
    const auto dE = Map< const Vector6d >( dStrain );

    // compute elastic stiffness
    const auto Cel = ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( E, nu );

    // handle zero strain increment
    if ( dE.isZero( 1e-14 ) ) {
      dS_dE = Cel;
      return;
    }

    // get current hardening variable
    double& kappa = managedStateVars->kappa;

    // isotropic hardening law
    auto fy = [&]( double kappa_ ) {
      return yieldStress + HLin * kappa_ + deltaYieldStress * ( 1. - std::exp( -delta * kappa_ ) );
    };

    // derivative of fy wrt dKappa
    auto dfy_ddKappa = [&]( double kappa_ ) { return HLin + deltaYieldStress * delta * std::exp( -delta * kappa_ ); };

    // yield function
    auto f = [&]( double rho_, double kappa_ ) { return rho_ - Constants::sqrt2_3 * fy( kappa_ ); };

    // compute elastic predictor
    const Vector6d trialStress = S + Cel * dE;

    using namespace ContinuumMechanics::VoigtNotation;
    const double rhoTrial = std::sqrt( 2. * Invariants::J2( trialStress ) );

    if ( f( rhoTrial, kappa ) >= 0.0 ) {
      // plastic step
      const double G = E / ( 2. * ( 1. + nu ) );

      auto g = [&]( double deltaKappa ) {
        return rhoTrial - Constants::sqrt6 * G * deltaKappa - Constants::sqrt2_3 * fy( kappa + deltaKappa );
      };

      // variables for return mapping
      int    counter    = 0;
      double dKappa     = 0;
      double dLambda    = 0;
      double dg_ddKappa = 0;

      // compute return mapping direction
      Vector6d n = ContinuumMechanics::VoigtNotation::IDev * trialStress / rhoTrial;

      while ( std::abs( g( dKappa ) ) > VonMisesConstants::innerNewtonTol ) {

        if ( counter == VonMisesConstants::nMaxInnerNewtonCycles ) {
          pNewDT = 0.5;
          return;
        }
        // compute derivative of g wrt kappa
        dg_ddKappa = -Constants::sqrt6 * G - Constants::sqrt2_3 * dfy_ddKappa( kappa + dKappa );

        // update dKappa and iteration counter
        dKappa -= g( dKappa ) / dg_ddKappa;
        counter += 1;
      }

      dLambda = Constants::sqrt3_2 * dKappa;

      // update material state
      S     = trialStress - 2. * G * dLambda * n;
      kappa = kappa + dKappa;

      // compute consistent tangent in Voigt Notation
      Matrix6d IDevHalfShear = ContinuumMechanics::VoigtNotation::IDev;
      IDevHalfShear.block< 6, 3 >( 0, 3 ) *= 0.5;

      dS_dE = Cel -
              2. * G * ( 1. / ( 1. + dfy_ddKappa( kappa ) / ( 3. * G ) ) - 2. * G * dLambda / rhoTrial ) *
                ( n * n.transpose() ) -
              4. * G * G * dLambda / rhoTrial * IDevHalfShear;
    }
    else {
      // elastic step
      S     = trialStress;
      dS_dE = Cel;
    }
  }

} // namespace Marmot::Materials
