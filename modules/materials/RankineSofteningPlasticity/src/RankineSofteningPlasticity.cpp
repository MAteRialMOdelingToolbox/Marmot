#include "Marmot/RankineSofteningPlasticity.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotTensor.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"
#include "Marmot/RankineSoftening.h"
#include "Marmot/RankineSofteningPlasticityConstants.h"
#include <Eigen/Eigenvalues>

namespace Marmot::Materials {

  using namespace Eigen;
  using namespace Marmot;

  double RankineSofteningPlasticity::getDensity( const double* stateVars ) const
  {
    if ( this->nMaterialProperties < 7 ) {
      throw std::runtime_error( std::string(
        Utilities::MakeString() << __PRETTY_FUNCTION__ << ": No density given! nMaterialProperties < 7." ) );
    }
    return this->materialProperties[6];
  }

  RankineSofteningPlasticity::RankineSofteningPlasticity( const double* materialProperties,
                                                          const int     nMaterialProperties,
                                                          const int     materialLabel )
    : MaterialHypoElastic( materialProperties, nMaterialProperties, materialLabel )
  {
    if ( this->nMaterialProperties < 6 ) {
      throw std::runtime_error( Utilities::MakeString()
                                << __PRETTY_FUNCTION__ << ": Insufficient material properties." );
    }
    stateLayout.add( "kappa", 1 );
    stateLayout.add( "acYield", 1 );
    stateLayout.finalize();
  }

  void RankineSofteningPlasticity::computeStress( state3D&              state,
                                                  Math::Matrix6d&       dStress_dStrain,
                                                  const Math::Vector6d& dStrain,
                                                  const timeInfo&       timeInfo ) const
  {
    updateStress( state, dStress_dStrain, dStrain, timeInfo );

    // a diverging iterate (e.g., of an arc-length iteration) may overflow the results (e.g., the elastic energy):
    // request a cutback instead of passing non-finite values to the solver
    if ( !( state.stress.allFinite() && dStress_dStrain.allFinite() && std::isfinite( state.elasticEnergyDensity ) &&
            std::isfinite( state.dissipation ) ) )
      throw Utilities::StressUpdateFailed(
        Utilities::MakeString() << __PRETTY_FUNCTION__ << ": non-finite results, e.g., from a diverging iterate" );
  }

  void RankineSofteningPlasticity::updateStress( state3D&              state,
                                                 Math::Matrix6d&       dStress_dStrain,
                                                 const Math::Vector6d& dStrain,
                                                 const timeInfo&       timeInfo ) const
  {
    using namespace RankineSofteningPlasticityConstants;
    using namespace Marmot::ContinuumMechanics::Voigt;
    using namespace Math::Tensors;
    using RankineSoftening::YieldSurfFlagArr;

    // elasticity parameters
    const double& E  = this->materialProperties[0];
    const double& nu = this->materialProperties[1];
    // isotropic hardening (softening for negative HLin, deltaYieldStress) law, identical to VonMisesModel
    const RankineSoftening::Hardening hardening{ this->materialProperties[2],
                                                 this->materialProperties[3],
                                                 this->materialProperties[4],
                                                 this->materialProperties[5] };

    // map to stress, strain and tangent
    Math::mVector6d S( state.stress.data() );
    Math::mMatrix6d dS_dE( dStress_dStrain.data() );
    const auto      dE = Map< const Math::Vector6d >( dStrain.data() );

    // compute elastic stiffness
    const Math::Matrix6d Cel = ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( E, nu );

    // get current hardening variable and yield flag
    double& kappa   = stateLayout.getAs< double& >( state.stateVars, "kappa" );
    double& acYield = stateLayout.getAs< double& >( state.stateVars, "acYield" );

    // handle zero strain increment
    if ( dE.isZero( 1e-14 ) ) {
      dS_dE   = Cel;
      acYield = 0;
      return;
    }

    const Math::Matrix6d CelInv = ContinuumMechanics::Elasticity::Isotropic::complianceTensor( E, nu );
    const double         psiOld = 0.5 * S.dot( CelInv * S );

    // compute elastic predictor
    const Math::Vector6d trialStress = S + Cel * dE;

    // spectral decomposition; the eigenvectors are preserved by the return mapping due to isotropy
    const SelfAdjointEigenSolver< Math::Matrix3d > eigenSolver( voigtToStress( trialStress ) );
    const Math::Vector3d&                          sigTrial = eigenSolver.eigenvalues();
    const Math::Matrix3d&                          V        = eigenSolver.eigenvectors();

    // initial guess for the active set: all violated yield surfaces
    const YieldSurfFlagArr activeSurfaces = ( sigTrial.transpose().array() - hardening.fy( kappa ) >= 0.0 );

    if ( !activeSurfaces.any() ) {
      // elastic step
      S       = trialStress;
      dS_dE   = Cel;
      acYield = 0;
      state.elasticEnergyDensity += 0.5 * S.dot( CelInv * S ) - psiOld;
      return;
    }

    // plastic step
    const double G      = E / ( 2. * ( 1. + nu ) );
    const double lambda = E * nu / ( ( 1. + nu ) * ( 1. - 2. * nu ) );

    // elastic stiffness in principal space
    const Math::Matrix3d Cp          = lambda * Math::Matrix3d::Ones() + 2. * G * Math::Matrix3d::Identity();
    const double         stressScale = std::max( 1.0, sigTrial.cwiseAbs().maxCoeff() );

    // results of the return mapping: principal stresses, increment of kappa, d sigma_A / d sigma^trial_B
    Math::Vector3d sig;
    double         dKappa = 0;
    Math::Matrix3d D;

    // backward Euler return mapping in principal space for a given set of active yield surfaces
    auto returnMappingAttempt = [&]( const YieldSurfFlagArr& active ) -> bool {
      // the system is solved for all three surfaces; inactive surfaces keep dLambda = 0 by identity rows
      const Math::Matrix3d activeMask = active.cast< double >().matrix().asDiagonal();
      const Math::Matrix3d inactiveId = Math::Matrix3d::Identity() - activeMask;

      Math::Vector3d dLambda = Math::Vector3d::Zero();

      auto residual = [&]() -> Math::Vector3d {
        return activeMask *
               ( sigTrial - Cp * dLambda - Math::Vector3d::Constant( hardening.fy( kappa + dLambda.sum() ) ) );
      };

      // negative Jacobian of the residual wrt dLambda
      auto jacobian = [&]() -> Math::Matrix3d {
        return activeMask * ( Cp + hardening.dfy_dKappa( kappa + dLambda.sum() ) * Math::Matrix3d::Ones() ) *
                 activeMask +
               inactiveId;
      };

      Math::Vector3d R       = residual();
      int            counter = 0;

      while ( R.norm() > innerNewtonTol * stressScale ) {
        if ( counter == nMaxInnerNewtonCycles )
          return false;

        const auto lu = jacobian().fullPivLu();
        if ( !lu.isInvertible() )
          return false;

        dLambda += lu.solve( R );
        R = residual();
        counter += 1;
      }

      // check consistency parameters and admissibility of the returned stress
      if ( ( dLambda.array() < 0.0 ).any() )
        return false;

      dKappa = dLambda.sum();
      sig    = sigTrial - Cp * dLambda;

      if ( ( sig.array() - hardening.fy( kappa + dKappa ) > yieldSurfaceResultTol * stressScale ).any() )
        return false;

      const auto lu = jacobian().fullPivLu();
      if ( !lu.isInvertible() )
        return false;

      // algorithmic tangent in principal space
      D = Math::Matrix3d::Identity() - Cp * activeMask * lu.solve( activeMask );

      return true;
    };

    RankineSoftening::searchActiveSurfaces( activeSurfaces, returnMappingAttempt, __PRETTY_FUNCTION__ );

    // update material state
    kappa += dKappa;
    acYield = 1;
    S       = stressToVoigt< double >( V * sig.asDiagonal() * V.transpose() );

    // consistent tangent from the spectral representation of the isotropic stress function
    const auto           basis = RankineSoftening::spectralBasis( V );
    const Math::Vector3d c     = RankineSoftening::spinCoefficients( sigTrial, sig, D, eigenValueTol * stressScale );
    const Math::Matrix9d dSig_dSigTrial = basis.E * D * basis.E.transpose() +
                                          basis.N * c.asDiagonal() * basis.N.transpose();

    // conversions between stress in Voigt notation and the flattened stress tensor
    static const auto [voigtToFlat, flatToVoigt] = [] {
      Matrix< double, 9, 6 > toFlat;
      Matrix< double, 6, 9 > toVoigt;
      for ( int i = 0; i < 6; i++ ) {
        const Math::Matrix3d s = voigtToStress( Math::Vector6d( Math::Vector6d::Unit( i ) ) );
        toFlat.col( i )        = flatten( s );
      }
      for ( int i = 0; i < 9; i++ ) {
        const Math::Vector9d f = Math::Vector9d::Unit( i );
        toVoigt.col( i )       = stressToVoigt< double >( as< 3, 3 >( f ) );
      }
      return std::pair{ toFlat, toVoigt };
    }();

    dS_dE = flatToVoigt * dSig_dSigTrial * voigtToFlat * Cel;

    state.elasticEnergyDensity += 0.5 * S.dot( CelInv * S ) - psiOld;
    state.dissipation += hardening.fy( kappa ) * dKappa;
  }

} // namespace Marmot::Materials
