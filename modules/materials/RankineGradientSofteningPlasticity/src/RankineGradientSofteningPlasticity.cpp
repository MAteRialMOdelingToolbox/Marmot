#include "Marmot/RankineGradientSofteningPlasticity.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotTensor.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"
#include "Marmot/RankineSofteningPlasticityConstants.h"
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>

namespace Marmot::Materials {

  using namespace Eigen;
  using namespace Marmot;

  namespace {
    // E, nu, fy0, HLin, deltaYieldStress, delta, nonlocal length, over-nonlocal parameter m
    constexpr int nRequiredProperties  = 8;
    constexpr int idxNonlocalLength    = 6;
    constexpr int idxM                 = 7;
    constexpr int idxDensity           = 8;
    constexpr int idxNonlocalViscosity = 9;
  } // namespace

  RankineGradientSofteningPlasticity::RankineGradientSofteningPlasticity( const double* materialProperties,
                                                                          int           nMaterialProperties,
                                                                          int           materialLabel )
    : Base( materialProperties, nMaterialProperties, materialLabel ), hardening( [&] {
        if ( nMaterialProperties < nRequiredProperties )
          throw std::runtime_error( Utilities::MakeString()
                                    << __PRETTY_FUNCTION__ << ": Insufficient material properties." );
        if ( !( materialProperties[idxNonlocalLength] >= 0. ) )
          throw std::invalid_argument( Utilities::MakeString()
                                       << __PRETTY_FUNCTION__ << ": nonlocal length must be >= 0." );
        return RankineSoftening::Hardening{ materialProperties[2],
                                            materialProperties[3],
                                            materialProperties[4],
                                            materialProperties[5] };
      }() )
  {
    stateLayout.add( "kappa", 1 );
    stateLayout.add( "acYield", 1 );
    stateLayout.add( "kappaNonlocal", 1 );
    stateLayout.add( "kappaWeighted", 1 );
    stateLayout.finalize();
  }

  double RankineGradientSofteningPlasticity::getDensity( const double* ) const
  {
    if ( nMaterialProperties <= idxDensity )
      throw std::runtime_error( Utilities::MakeString()
                                << __PRETTY_FUNCTION__ << ": No density given! nMaterialProperties < 9." );
    return materialProperties[idxDensity];
  }

  std::vector< double > RankineGradientSofteningPlasticity::getNonlocalViscosity( const double* ) const
  {
    return { nMaterialProperties > idxNonlocalViscosity ? materialProperties[idxNonlocalViscosity] : 0.0 };
  }

  void RankineGradientSofteningPlasticity::computeStress( response& res, tangents& tan, const increment& inc ) const
  {
    // a NaN nonlocal field would not enter the checked results, since it floors the yield stress at zero
    if ( !inc.K.allFinite() )
      throw Utilities::StressUpdateFailed( Utilities::MakeString()
                                           << __PRETTY_FUNCTION__ << ": non-finite nonlocal field" );

    updateStress( res, tan, inc );

    // a diverging iterate may overflow the results (e.g., the elastic energy): request a cutback instead of passing
    // non-finite values to the solver
    if ( !( res.stress.allFinite() && res.KLocal.allFinite() && std::isfinite( res.elasticEnergyDensity ) &&
            std::isfinite( res.dissipation ) && tan.dStressddStrain.allFinite() && tan.dStressddK.allFinite() &&
            tan.dKLocalddStrain.allFinite() && tan.dKLocalddK.allFinite() ) )
      throw Utilities::StressUpdateFailed(
        Utilities::MakeString() << __PRETTY_FUNCTION__ << ": non-finite results, e.g., from a diverging iterate" );
  }

  void RankineGradientSofteningPlasticity::updateStress( response& res, tangents& tan, const increment& inc ) const
  {
    using namespace RankineSofteningPlasticityConstants;
    using namespace Marmot::ContinuumMechanics::Voigt;
    using namespace Math::Tensors;
    using RankineSoftening::YieldSurfFlagArr;

    const double& E              = materialProperties[0];
    const double& nu             = materialProperties[1];
    const double& nonlocalLength = materialProperties[idxNonlocalLength];
    const double& m              = materialProperties[idxM];

    Math::mVector6d S( res.stress.data() );
    const auto      dE = Map< const Math::Vector6d >( inc.dStrain.data() );
    const double    K  = inc.K( 0 );

    double& kappa         = stateLayout.getAs< double& >( res.stateVars, "kappa" );
    double& acYield       = stateLayout.getAs< double& >( res.stateVars, "acYield" );
    double& kappaNonlocal = stateLayout.getAs< double& >( res.stateVars, "kappaNonlocal" );
    double& kappaWeighted = stateLayout.getAs< double& >( res.stateVars, "kappaWeighted" );

    res.c( 0 ) = nonlocalLength * nonlocalLength;
    tan.dcddK.setZero();
    tan.d2cddK2.setZero();

    // yield stress and its derivative wrt the over-nonlocal softening variable, for a local softening variable
    auto weighted = [&]( double kappa_ ) { return m * K + ( 1. - m ) * kappa_; };
    auto fy       = [&]( double kappa_ ) { return hardening.fy( weighted( kappa_ ) ); };
    auto dfy      = [&]( double kappa_ ) { return hardening.dfy_dKappa( weighted( kappa_ ) ); };

    const Math::Matrix6d Cel    = ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( E, nu );
    const Math::Matrix6d CelInv = ContinuumMechanics::Elasticity::Isotropic::complianceTensor( E, nu );
    const double         psiOld = 0.5 * S.dot( CelInv * S );

    // elastic predictor; also without a strain increment, since the yield stress depends on the nonlocal variable
    const Math::Vector6d trialStress = S + Cel * dE;

    // spectral decomposition; the eigenvectors are preserved by the return mapping due to isotropy
    const SelfAdjointEigenSolver< Math::Matrix3d > eigenSolver( voigtToStress( trialStress ) );
    const Math::Vector3d&                          sigTrial = eigenSolver.eigenvalues();
    const Math::Matrix3d&                          V        = eigenSolver.eigenvectors();

    // initial guess for the active set: all violated yield surfaces
    const YieldSurfFlagArr activeSurfaces = ( sigTrial.transpose().array() - fy( kappa ) >= 0.0 );

    kappaNonlocal = K;

    if ( !activeSurfaces.any() ) {
      // elastic step
      S                   = trialStress;
      acYield             = 0;
      kappaWeighted       = weighted( kappa );
      tan.dStressddStrain = Cel;
      tan.dStressddK.setZero();
      res.KLocal( 0 ) = kappa;
      tan.dKLocalddStrain.setZero();
      tan.dKLocalddK.setZero();
      res.elasticEnergyDensity += 0.5 * S.dot( CelInv * S ) - psiOld;
      return;
    }

    // plastic step
    const double G      = E / ( 2. * ( 1. + nu ) );
    const double lambda = E * nu / ( ( 1. + nu ) * ( 1. - 2. * nu ) );

    // elastic stiffness in principal space
    const Math::Matrix3d Cp          = lambda * Math::Matrix3d::Ones() + 2. * G * Math::Matrix3d::Identity();
    const double         stressScale = std::max( 1.0, sigTrial.cwiseAbs().maxCoeff() );

    // results of the return mapping: principal stresses, increment of kappa, d sigma_A / d sigma^trial_B, and the
    // derivatives of the principal stresses and of kappa wrt the nonlocal variable
    Math::Vector3d sig;
    double         dKappa = 0;
    Math::Matrix3d D;
    Math::Vector3d dSig_dK;
    double         dKappa_dK = 0;

    // backward Euler return mapping in principal space for a given set of active yield surfaces
    auto returnMappingAttempt = [&]( const YieldSurfFlagArr& active ) -> bool {
      // the system is solved for all three surfaces; inactive surfaces keep dLambda = 0 by identity rows
      const Math::Matrix3d activeMask = active.cast< double >().matrix().asDiagonal();
      const Math::Matrix3d inactiveId = Math::Matrix3d::Identity() - activeMask;

      Math::Vector3d dLambda = Math::Vector3d::Zero();

      auto residual = [&]() -> Math::Vector3d {
        return activeMask * ( sigTrial - Cp * dLambda - Math::Vector3d::Constant( fy( kappa + dLambda.sum() ) ) );
      };

      // negative Jacobian of the residual wrt dLambda; the local softening variable enters the yield stress by 1 - m
      auto jacobian = [&]() -> Math::Matrix3d {
        return activeMask * ( Cp + ( 1. - m ) * dfy( kappa + dLambda.sum() ) * Math::Matrix3d::Ones() ) * activeMask +
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

      if ( ( sig.array() - fy( kappa + dKappa ) > yieldSurfaceResultTol * stressScale ).any() )
        return false;

      const auto lu = jacobian().fullPivLu();
      if ( !lu.isInvertible() )
        return false;

      // algorithmic tangent in principal space
      D = Math::Matrix3d::Identity() - Cp * activeMask * lu.solve( activeMask );

      // the nonlocal variable shifts the yield stress of the active surfaces by m fy'
      const Math::Vector3d dLambda_dK = -m * dfy( kappa + dKappa ) * lu.solve( activeMask * Math::Vector3d::Ones() );
      dSig_dK                         = -Cp * dLambda_dK;
      dKappa_dK                       = dLambda_dK.sum();

      return true;
    };

    RankineSoftening::searchActiveSurfaces( activeSurfaces, returnMappingAttempt, __PRETTY_FUNCTION__ );

    // update material state
    kappa += dKappa;
    acYield       = 1;
    kappaWeighted = weighted( kappa );
    S             = stressToVoigt< double >( V * sig.asDiagonal() * V.transpose() );

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

    tan.dStressddStrain = flatToVoigt * dSig_dSigTrial * voigtToFlat * Cel;
    tan.dStressddK      = flatToVoigt * ( basis.E * dSig_dK );

    // dKappa = tr( dE - Cel^-1 dStress ), since the flow directions have unit trace; I is the Voigt identity
    res.KLocal( 0 )              = kappa;
    tan.dKLocalddStrain.row( 0 ) = ( I - tan.dStressddStrain.transpose() * ( CelInv * I ) ).transpose();
    tan.dKLocalddK( 0, 0 )       = dKappa_dK;

    res.elasticEnergyDensity += 0.5 * S.dot( CelInv * S ) - psiOld;
    res.dissipation += fy( kappa ) * dKappa;
  }

} // namespace Marmot::Materials
