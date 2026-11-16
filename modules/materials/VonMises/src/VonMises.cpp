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
      throw std::runtime_error( std::string(
        Marmot::Utilities::MakeString() << __PRETTY_FUNCTION__ << ": No density given! nMaterialProperties < 7." ) );
    }
    return this->materialProperties[6];
  }

  VonMisesModel::VonMisesModel( const double* materialProperties,
                                const int     nMaterialProperties,
                                const int     materialLabel )
    : Marmot::Materials::MaterialHypoElastic( materialProperties, nMaterialProperties, materialLabel )
  {
    stateLayout.add( "kappa", 1 );
    stateLayout.finalize();
  }

  void VonMisesModel::computeStress( state3D&                      state,
                                     Marmot::Math::Matrix6d&       dStress_dStrain,
                                     const Marmot::Math::Vector6d& dStrain,
                                     const timeInfo&               timeInfo ) const
  {
    computeStressImplicit( state, dStress_dStrain, dStrain, timeInfo, nullptr, nullptr );
  }

  void VonMisesModel::computeStressWithEnergyTangents( state3D&                      state,
                                                       Marmot::Math::Matrix6d&       dStress_dStrain,
                                                       const Marmot::Math::Vector6d& dStrain,
                                                       const timeInfo&               timeInfo,
                                                       Marmot::Math::Vector6d&       dElasticEnergyDensity_dStrain,
                                                       Marmot::Math::Vector6d&       dDissipation_dStrain ) const
  {
    computeStressImplicit( state,
                           dStress_dStrain,
                           dStrain,
                           timeInfo,
                           &dElasticEnergyDensity_dStrain,
                           &dDissipation_dStrain );
  }

  void VonMisesModel::computeStressImplicit( state3D&                      state,
                                             Marmot::Math::Matrix6d&       dStress_dStrain,
                                             const Marmot::Math::Vector6d& dStrain,
                                             const timeInfo&               timeInfo,
                                             Marmot::Math::Vector6d*       dPsi_dE,
                                             Marmot::Math::Vector6d*       dD_dE ) const
  {
    const bool energyTangents = dPsi_dE && dD_dE;

    // elasticity parameters
    const double& E  = this->materialProperties[0];
    const double& nu = this->materialProperties[1];
    // plasticity parameters
    const double& yieldStress      = this->materialProperties[2];
    const double& HLin             = this->materialProperties[3];
    const double& deltaYieldStress = this->materialProperties[4];
    const double& delta            = this->materialProperties[5];

    // map to stress, strain and tangent
    Marmot::Math::mVector6d S( state.stress.data() );
    Marmot::Math::mMatrix6d dS_dE( dStress_dStrain.data() );
    const auto              dE = Map< const Marmot::Math::Vector6d >( dStrain.data() );

    // compute elastic stiffness and compliance
    const auto                   Cel    = ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( E, nu );
    const Marmot::Math::Matrix6d CelInv = ContinuumMechanics::Elasticity::Isotropic::complianceTensor( E, nu );

    // handle zero strain increment
    if ( dE.isZero( 1e-14 ) ) {
      dS_dE = Cel;
      if ( energyTangents ) {
        *dPsi_dE = S;
        dD_dE->setZero();
      }
      return;
    }

    const double psiOld = 0.5 * S.dot( CelInv * S );

    // get current hardening variable
    double& kappa = stateLayout.getAs< double& >( state.stateVars, "kappa" );

    // isotropic hardening law, yield stress floored at zero
    auto fyUnfloored = [&]( double kappa_ ) {
      return yieldStress + HLin * kappa_ + deltaYieldStress * ( 1. - std::exp( -delta * kappa_ ) );
    };
    auto fy = [&]( double kappa_ ) { return std::max( fyUnfloored( kappa_ ), 0. ); };

    // derivative of fy wrt dKappa
    auto dfy_ddKappa = [&]( double kappa_ ) {
      return fyUnfloored( kappa_ ) > 0. ? HLin + deltaYieldStress * delta * std::exp( -delta * kappa_ ) : 0.;
    };

    // yield function
    auto f = [&]( double rho_, double kappa_ ) { return rho_ - Marmot::Math::Constants::sqrt2_3 * fy( kappa_ ); };

    // compute elastic predictor
    const Marmot::Math::Vector6d trialStress = S + Cel * dE;

    using namespace ContinuumMechanics::Voigt;
    const double rhoTrial = std::sqrt( 2. * ContinuumMechanics::Invariants::J2( trialStress ) );

    if ( f( rhoTrial, kappa ) >= 0.0 && rhoTrial > 0.0 ) {
      // plastic step
      const double G = E / ( 2. * ( 1. + nu ) );

      auto g = [&]( double deltaKappa ) {
        return rhoTrial - Marmot::Math::Constants::sqrt6 * G * deltaKappa -
               Marmot::Math::Constants::sqrt2_3 * fy( kappa + deltaKappa );
      };

      // variables for return mapping
      int    counter    = 0;
      double dKappa     = 0;
      double dLambda    = 0;
      double dg_ddKappa = 0;

      // compute return mapping direction
      Marmot::Math::Vector6d n = ContinuumMechanics::Voigt::IDev * trialStress / rhoTrial;

      while ( std::abs( g( dKappa ) ) > Constants::innerNewtonTol ) {

        if ( counter == Constants::nMaxInnerNewtonCycles ) {
          throw Marmot::Utilities::StressUpdateFailed(
            "return mapping failed to converge in VonMisesModel::computeStress" );
        }
        // compute derivative of g wrt kappa
        dg_ddKappa = -Marmot::Math::Constants::sqrt6 * G -
                     Marmot::Math::Constants::sqrt2_3 * dfy_ddKappa( kappa + dKappa );

        // update dKappa and iteration counter
        dKappa -= g( dKappa ) / dg_ddKappa;
        counter += 1;
      }

      dLambda = Marmot::Math::Constants::sqrt3_2 * dKappa;

      // update material state
      S     = trialStress - 2. * G * dLambda * n;
      kappa = kappa + dKappa;

      // compute consistent tangent in Voigt Notation
      Marmot::Math::Matrix6d IDevHalfShear = ContinuumMechanics::Voigt::IDev;
      IDevHalfShear.block< 6, 3 >( 0, 3 ) *= 0.5;

      dS_dE = Cel -
              2. * G * ( 1. / ( 1. + dfy_ddKappa( kappa ) / ( 3. * G ) ) - 2. * G * dLambda / rhoTrial ) *
                ( n * n.transpose() ) -
              4. * G * G * dLambda / rhoTrial * IDevHalfShear;
      state.elasticEnergyDensity += 0.5 * S.dot( CelInv * S ) - psiOld;
      // plastic work sigma : dEpsP = rho dLambda = fy dKappa, with kappa already updated
      state.dissipation += fy( kappa ) * dKappa;

      // energy tangents; dKappa follows from g = 0 with d rhoTrial / d dE = 2 G n
      if ( energyTangents ) {
        const double dDKappa_dRhoTrial = 1. / ( Marmot::Math::Constants::sqrt6 * G +
                                                Marmot::Math::Constants::sqrt2_3 * dfy_ddKappa( kappa ) );
        *dPsi_dE                       = dS_dE.transpose() * ( CelInv * S );
        *dD_dE = ( fy( kappa ) + dfy_ddKappa( kappa ) * dKappa ) * dDKappa_dRhoTrial * 2. * G * n;
      }
    }
    else {
      // elastic step
      S     = trialStress;
      dS_dE = Cel;
      state.elasticEnergyDensity += 0.5 * S.dot( CelInv * S ) - psiOld;
      if ( energyTangents ) {
        *dPsi_dE = S;
        dD_dE->setZero();
      }
    }
  }

  void VonMisesModel::computeStressExplicit( state3D&                      state,
                                             const Marmot::Math::Vector6d& dStrain,
                                             const timeInfo&               timeInfo ) const

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
    Marmot::Math::mVector6d S( state.stress.data() );
    const auto              dE = dStrain;

    // compute elastic stiffness and compliance
    const auto                   Cel    = ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( E, nu );
    const Marmot::Math::Matrix6d CelInv = ContinuumMechanics::Elasticity::Isotropic::complianceTensor( E, nu );

    // handle zero strain increment
    if ( dE.isZero( 1e-14 ) ) {
      return;
    }

    const double psiOld = 0.5 * S.dot( CelInv * S );

    // get current hardening variable
    double& kappa = stateLayout.getAs< double& >( state.stateVars, "kappa" );

    // isotropic hardening law, yield stress floored at zero
    auto fyUnfloored = [&]( double kappa_ ) {
      return yieldStress + HLin * kappa_ + deltaYieldStress * ( 1. - std::exp( -delta * kappa_ ) );
    };
    auto fy = [&]( double kappa_ ) { return std::max( fyUnfloored( kappa_ ), 0. ); };

    // derivative of fy wrt dKappa
    auto dfy_ddKappa = [&]( double kappa_ ) {
      return fyUnfloored( kappa_ ) > 0. ? HLin + deltaYieldStress * delta * std::exp( -delta * kappa_ ) : 0.;
    };

    // yield function
    auto f = [&]( double rho_, double kappa_ ) { return rho_ - Marmot::Math::Constants::sqrt2_3 * fy( kappa_ ); };

    // compute elastic predictor
    const Marmot::Math::Vector6d trialStress = S + Cel * dE;

    using namespace ContinuumMechanics::Voigt;
    const double rhoTrial = std::sqrt( 2. * ContinuumMechanics::Invariants::J2( trialStress ) );

    if ( f( rhoTrial, kappa ) >= 0.0 && rhoTrial > 0.0 ) {
      // plastic step
      const double G = E / ( 2. * ( 1. + nu ) );

      auto g = [&]( double deltaKappa ) {
        return rhoTrial - Marmot::Math::Constants::sqrt6 * G * deltaKappa -
               Marmot::Math::Constants::sqrt2_3 * fy( kappa + deltaKappa );
      };

      // variables for return mapping
      int    counter    = 0;
      double dKappa     = 0;
      double dLambda    = 0;
      double dg_ddKappa = 0;

      // compute return mapping direction
      Marmot::Math::Vector6d n = ContinuumMechanics::Voigt::IDev * trialStress / rhoTrial;

      double g_val = g( dKappa );

      // always perform 5 iterations for explicit scheme
      // this allows to avoid if-conditions inside the time stepping loop
      // and improves performance
      while ( counter < 5 ) {

        // compute derivative of g wrt kappa
        dg_ddKappa = -Marmot::Math::Constants::sqrt6 * G -
                     Marmot::Math::Constants::sqrt2_3 * dfy_ddKappa( kappa + dKappa );

        // update dKappa and iteration counter
        dKappa -= g_val / dg_ddKappa;
        g_val = g( dKappa );
        counter += 1;
      }

      if ( std::abs( g_val ) > Constants::innerNewtonTol ) {
        throw Marmot::Utilities::StressUpdateFailed(
          "return mapping failed to converge in VonMisesModel::computeStressExplicit" );
      }

      dLambda = Marmot::Math::Constants::sqrt3_2 * dKappa;

      // update material state
      S     = trialStress - 2. * G * dLambda * n;
      kappa = kappa + dKappa;
      state.elasticEnergyDensity += 0.5 * S.dot( CelInv * S ) - psiOld;
      // plastic work sigma : dEpsP = rho dLambda = fy dKappa, with kappa already updated
      state.dissipation += fy( kappa ) * dKappa;
    }
    else {
      // elastic step
      S = trialStress;
      state.elasticEnergyDensity += 0.5 * S.dot( CelInv * S ) - psiOld;
    }
  }

} // namespace Marmot::Materials
