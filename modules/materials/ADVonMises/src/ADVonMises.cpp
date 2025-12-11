#include "Marmot/ADVonMises.h"
#include "Marmot/ADVonMisesConstants.h"
#include "Marmot/MarmotConstants.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"
#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>

namespace Marmot::Materials {

  using namespace Marmot;
  using namespace autodiff;
  using namespace Eigen;
  using namespace ContinuumMechanics::Elasticity;

  ADVonMises::ADVonMises( const double* materialProperties, int nMaterialProperties, int materialNumber )
    : MarmotMaterialHypoElasticAD::MarmotMaterialHypoElasticAD( materialProperties,
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
    assert( nMaterialProperties == 6 );
    initializeStateLayout();
  }

  void ADVonMises::computeStressAD( state3DAD& state, const autodiff::dual* dStrain, const timeInfo& timeInfo ) const
  {
    mVector6dual            S( state.stress );
    const mVector6dualConst dE( dStrain );

    // compute elastic stiffness
    const Matrix6d Cel = ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( E, nu );

    // get current hardening variable
    double& kappa = stateLayout.getAs< double& >( state.stateVars, "kappa" );

    // compute elastic predictor
    const Vector6dual trialStress = S + Cel * dE;

    using namespace ContinuumMechanics::VoigtNotation;
    const dual rhoTrial = sqrt( 2. * Invariants::J2( trialStress ) );

    if ( Math::makeReal( f( rhoTrial, kappa ) ) >= 0.0 ) {

      // variables for return mapping
      size_t counter = 0;
      dual   dKappa( 0.0 );
      dual   dLambda( 0.0 );
      double dg_ddKappa( 0.0 );

      while ( abs( g( (double)rhoTrial, kappa, (double)dKappa ) ) > ADVonMisesConstants::innerNewtonTol ) {

        if ( counter == ADVonMisesConstants::nMaxInnerNewtonCycles ) {
          throw std::runtime_error( MakeString()
                                    << __PRETTY_FUNCTION__
                                    << ": Return mapping did not converge within maximum number of iterations!" );
        }

        // compute derivative of g wrt kappa at constant rhoTrial
        dual dKappa_( dKappa );
        seed< 1 >( dKappa_, 1.0 );
        dg_ddKappa = derivative< 1 >( g( dual( rhoTrial.val ), kappa, dKappa_ ) );

        // update dKappa and iteration counter
        dKappa -= g( rhoTrial, kappa, dKappa ) / dg_ddKappa;
        counter += 1;
      }

      // compute plastic corrector
      dLambda = Constants::sqrt3_2 * dKappa;

      // compute return mapping direction
      const Vector6dual n = ContinuumMechanics::VoigtNotation::IDev * trialStress / rhoTrial;

      // update stress and hardening variable
      S = trialStress - 2. * G * dLambda * n;
      kappa += dKappa.val;
    }
    else {
      // elastic step
      S = trialStress;
    }
  }
} // namespace Marmot::Materials
