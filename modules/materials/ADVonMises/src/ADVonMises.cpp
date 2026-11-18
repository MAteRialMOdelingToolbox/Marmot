#include "Marmot/ADVonMises.h"
#include "Marmot/MarmotAutomaticDifferentiationTypes.h"
#include "Marmot/MarmotConstants.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"

namespace Marmot::Materials {

  using namespace Marmot;
  using namespace Marmot::Math::Differentiation::Automatic;
  using namespace Eigen;
  using namespace ContinuumMechanics::Elasticity;

  double ADVonMises::getDensity( const double* stateVars ) const
  {
    if ( nMaterialProperties < 7 ) {
      throw std::runtime_error( Marmot::Utilities::MakeString()
                                << __PRETTY_FUNCTION__ << ": Density not provided in material properties array!" );
    }
    else {
      return materialProperties[6];
    }
  }

  ADVonMises::ADVonMises( const double* materialProperties, int nMaterialProperties, int materialNumber )
    : Marmot::Materials::MaterialHypoElasticAD::MaterialHypoElasticAD( materialProperties,
                                                                       nMaterialProperties,
                                                                       materialNumber ),
      E( materialProperties[0] ),
      nu( materialProperties[1] ),
      yieldStress( materialProperties[2] ),
      HLin( materialProperties[3] ),
      deltaYieldStress( materialProperties[4] ),
      delta( materialProperties[5] ),
      G( E / ( 2. * ( 1. + nu ) ) )
  {
    stateLayout.add( "kappa", 1 );
    stateLayout.finalize();
  }

  void ADVonMises::computeStressAD( state3DAD& state, const Vector6AD& dStrain, const timeInfo& timeInfo ) const
  {
    Eigen::Map< Vector6AD >       S( state.stress.data() );
    Eigen::Map< const Vector6AD > dE( dStrain.data() );

    // compute elastic stiffness
    const Marmot::Math::Matrix6d Cel = ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( E, nu );

    // get current hardening variable
    double& kappa = stateLayout.getAs< double& >( state.stateVars, "kappa" );

    // compute elastic predictor
    const Vector6AD trialStress = S + Cel * dE;

    using namespace ContinuumMechanics::Voigt;
    const ADScalar rhoTrial = sqrt( 2. * ContinuumMechanics::Invariants::J2( trialStress ) );

    if ( Math::Tensors::makeReal( f( rhoTrial, kappa ) ) >= 0.0 ) {

      // variables for return mapping
      size_t   counter = 0;
      ADScalar dKappa( 0.0 );
      ADScalar dLambda( 0.0 );
      double   dg_ddKappa( 0.0 );

      while ( abs( g( value( rhoTrial ), kappa, value( dKappa ) ) ) > Constants::innerNewtonTol * value( rhoTrial ) ) {

        if ( counter == Constants::nMaxInnerNewtonCycles ) {
          throw Marmot::Utilities::StressUpdateFailed(
            Marmot::Utilities::MakeString()
            << __PRETTY_FUNCTION__ << ": Return mapping did not converge within maximum number of iterations!" );
        }

        // compute derivative of g wrt kappa at constant rhoTrial
        dual dKappa_( value( dKappa ) );
        seed< 1 >( dKappa_, 1.0 );
        dg_ddKappa = derivative< 1 >( g( dual( value( rhoTrial ) ), kappa, dKappa_ ) );

        // update dKappa and iteration counter
        dKappa -= g( rhoTrial, kappa, dKappa ) / dg_ddKappa;
        counter += 1;
      }

      // compute plastic corrector
      dLambda = Marmot::Math::Constants::sqrt3_2 * dKappa;

      // compute return mapping direction
      const Vector6AD n = ContinuumMechanics::Voigt::IDev * trialStress / rhoTrial;

      // update stress and hardening variable
      S          = trialStress - 2. * G * dLambda * n;
      auto dEpsP = dLambda * n;
      auto dEpsE = dE - dEpsP;
      auto dS    = Cel * dEpsE;
      kappa += value( dKappa );
      auto stressAvg = S - 0.5 * Cel * dEpsE;
      state.elasticEnergyDensity += Math::Tensors::makeReal( dEpsE.dot( stressAvg ) );
      state.dissipation += Math::Tensors::makeReal( dEpsP.dot( dS ) );
    }
    else {
      // elastic step
      S              = trialStress;
      auto stressAvg = S - 0.5 * Cel * dE;
      state.elasticEnergyDensity += Math::Tensors::makeReal( stressAvg.dot( dE ) );
      state.dissipation += 0.0;
    }
  }
} // namespace Marmot::Materials
