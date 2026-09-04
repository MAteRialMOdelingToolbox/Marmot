#include "Marmot/VonMises.h"
#include "Marmot/MarmotConstants.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotTypedefs.h"

namespace Marmot::Materials {

  using namespace Eigen;
  using namespace Marmot;

  double VonMisesModel::getDensity( const double* stateVars ) const
  {
    if ( this->nMaterialProperties < 7 ) {
      throw std::runtime_error(
        std::string( MakeString() << __PRETTY_FUNCTION__ << ": No density given! nMaterialProperties < 7." ) );
    }
    return this->materialProperties[6];
  }

  VonMisesModel::VonMisesModel( const double* materialProperties,
                                const int     nMaterialProperties,
                                const int     materialLabel )
    : MarmotMaterialHypoElastic( materialProperties, nMaterialProperties, materialLabel )
  {
    stateLayout.add( "kappa", 1 );
    stateLayout.finalize();
  }

  void VonMisesModel::computeStress( state3D&        state,
                                     Matrix6d&       dStress_dStrain,
                                     const Vector6d& dStrain,
                                     const timeInfo& timeInfo ) const

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
    mVector6d  S( state.stress.data() );
    mMatrix6d  dS_dE( dStress_dStrain.data() );
    const auto dE = Map< const Vector6d >( dStrain.data() );

    // compute elastic stiffness
    const auto Cel = ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( E, nu );

    // handle zero strain increment
    if ( dE.isZero( 1e-14 ) ) {
      dS_dE = Cel;
      return;
    }

    // get current hardening variable
    double& kappa = stateLayout.getAs< double& >( state.stateVars, "kappa" );

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

    using namespace ContinuumMechanics::Voigt;
    using namespace ContinuumMechanics::Invariants;
    const double rhoTrial = std::sqrt( 2. * J2( trialStress ) );

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
      Vector6d n = ContinuumMechanics::Voigt::IDev * trialStress / rhoTrial;

      while ( std::abs( g( dKappa ) ) > innerNewtonTol ) {

        if ( counter == nMaxInnerNewtonCycles ) {
          throw StressUpdateFailed( "return mapping failed to converge in VonMisesModel::computeStress" );
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
      Matrix6d IDevHalfShear = ContinuumMechanics::Voigt::IDev;
      IDevHalfShear.block< 6, 3 >( 0, 3 ) *= 0.5;

      dS_dE = Cel -
              2. * G * ( 1. / ( 1. + dfy_ddKappa( kappa ) / ( 3. * G ) ) - 2. * G * dLambda / rhoTrial ) *
                ( n * n.transpose() ) -
              4. * G * G * dLambda / rhoTrial * IDevHalfShear;
      state.elasticEnergyDensity += 0.5 * dE.dot( Cel * dE ) - 2. * G * dLambda * n.dot( dE );
      state.dissipation += 2. * G * dLambda * n.dot( dE );
    }
    else {
      // elastic step
      S     = trialStress;
      dS_dE = Cel;
      state.elasticEnergyDensity += 0.5 * dE.dot( Cel * dE );
    }
  }

  void VonMisesModel::computeStressExplicit( state3D& state, const Vector6d& dStrain, const timeInfo& timeInfo ) const

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
    mVector6d  S( state.stress.data() );
    const auto dE = dStrain;

    // compute elastic stiffness
    const auto Cel = ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( E, nu );

    // handle zero strain increment
    if ( dE.isZero( 1e-14 ) ) {
      return;
    }

    // get current hardening variable
    double& kappa = stateLayout.getAs< double& >( state.stateVars, "kappa" );

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

    using namespace ContinuumMechanics::Voigt;
    using namespace ContinuumMechanics::Invariants;
    const double rhoTrial = std::sqrt( 2. * J2( trialStress ) );

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
      Vector6d n = ContinuumMechanics::Voigt::IDev * trialStress / rhoTrial;

      double g_val = g( dKappa );

      // always perform 5 iterations for explicit scheme
      // this allows to avoid if-conditions inside the time stepping loop
      // and improves performance
      while ( counter < 5 ) {

        // compute derivative of g wrt kappa
        dg_ddKappa = -Constants::sqrt6 * G - Constants::sqrt2_3 * dfy_ddKappa( kappa + dKappa );

        // update dKappa and iteration counter
        dKappa -= g_val / dg_ddKappa;
        g_val = g( dKappa );
        counter += 1;
      }

      if ( std::abs( g_val ) > innerNewtonTol ) {
        throw Marmot::StressUpdateFailed( "return mapping failed to converge in VonMisesModel::computeStressExplicit" );
      }

      dLambda = Constants::sqrt3_2 * dKappa;

      // update material state
      S     = trialStress - 2. * G * dLambda * n;
      kappa = kappa + dKappa;
      state.elasticEnergyDensity += 0.5 * dE.dot( Cel * dE ) - 2. * G * dLambda * n.dot( dE );
      state.dissipation += 2. * G * dLambda * n.dot( dE );
    }
    else {
      // elastic step
      S = trialStress;
      state.elasticEnergyDensity += 0.5 * dE.dot( Cel * dE );
    }
  }

} // namespace Marmot::Materials
