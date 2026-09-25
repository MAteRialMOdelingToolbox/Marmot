#include "Marmot/GradientEnhancedFiniteStrainDruckerPragerPlasticity.h"
#include "Marmot/GradientEnhancedFiniteStrainDruckerPragerConstants.h"
#include <Eigen/Dense>
#include <cmath>

namespace Marmot::Materials {

  using namespace Eigen;
  namespace Constants = GradientEnhancedFiniteStrainDruckerPragerConstants;

  using Plasticity = GradientEnhancedFiniteStrainDruckerPragerPlasticity;

  Plasticity::GradientEnhancedFiniteStrainDruckerPragerPlasticity( const MaterialParameters& materialParameters,
                                                                   const ModelParameters&    modelParameters )
    : materialParameters( materialParameters ), modelParameters( modelParameters )
  {
  }

  Plasticity::StressState Plasticity::computeStressState( const MaterialState& state ) const
  {
    const auto& [K, G]   = materialParameters;
    const Vector3d& eps  = state.elasticLogStrainPrincipal;
    const double    th   = eps.sum();
    const Vector3d  E    = ( 2.0 * ( eps.array() - th / 3.0 ) ).exp();
    const double    Emid = E.sum() / 3.0;

    StressState s;
    s.mandelPrincipal = Vector3d::Constant( 0.5 * K * std::sinh( 2.0 * th ) ) + G * ( E.array() - Emid ).matrix();
    return s;
  }

  Matrix3d Plasticity::dStress_dElasticLogStrain( const MaterialState& state ) const
  {
    const auto& [K, G]   = materialParameters;
    const Vector3d& eps  = state.elasticLogStrainPrincipal;
    const double    th   = eps.sum();
    const Vector3d  E    = ( 2.0 * ( eps.array() - th / 3.0 ) ).exp();
    const double    Emid = E.sum() / 3.0;

    Matrix3d D = Matrix3d::Constant( K * std::cosh( 2.0 * th ) );
    for ( int a = 0; a < 3; a++ )
      for ( int c = 0; c < 3; c++ )
        D( a, c ) += G * ( 2.0 * E( a ) * ( ( a == c ) - 1. / 3 ) - 2. / 3 * ( E( c ) - Emid ) );
    return D;
  }

  double Plasticity::yieldFunction( const StressState& stress, double alphaP ) const
  {
    const Vector3d& S   = stress.mandelPrincipal;
    const double    p   = S.sum() / 3.0;
    const Vector3d  dev = S.array() - p;
    const double    q   = std::sqrt( 0.5 * dev.squaredNorm() );
    return q + modelParameters.eta * p - modelParameters.xi * cohesion( alphaP );
  }

  bool Plasticity::checkIfYielding( const StressState& stress, const MaterialState& state ) const
  {
    return yieldFunction( stress, state.alphaP ) > 0.0;
  }

  Plasticity::ReturnMapResult Plasticity::performReturnMapping( const MaterialState& trialState ) const
  {
    ReturnMapResult result{ computeStressState( trialState ), trialState, Vector3d::Zero(), ReturnMode::Elastic };

    if ( !checkIfYielding( result.newStressState, trialState ) )
      return result;

    if ( returnToCone( trialState, result ) )
      return result;

    if ( returnToApex( trialState, result ) )
      return result;

    throw ReturnMappingFailedException();
  }

  bool Plasticity::returnToCone( const MaterialState& trialState, ReturnMapResult& result ) const
  {
    const auto& [eta, xi, etaBar, c0, H] = modelParameters;
    const double G                       = materialParameters.G; // scales the yield residual to a strain

    const Vector3d& epsTrial  = trialState.elasticLogStrainPrincipal;
    const Vector3d  SdevTrial = computeStressState( trialState ).mandelPrincipal.array() -
                               computeStressState( trialState ).mandelPrincipal.mean();

    const Matrix3d Pdev = Matrix3d::Identity() - Matrix3d::Constant( 1. / 3 );

    MaterialState state   = trialState;
    double        dLambda = 0.0;

    for ( int i = 0; i < Constants::nMaxInnerNewtonCycles; i++ ) {

      const Vector3d S   = computeStressState( state ).mandelPrincipal;
      const double   p   = S.mean();
      const Vector3d dev = S.array() - p;
      const double   q   = std::sqrt( 0.5 * dev.squaredNorm() );

      if ( q <= Constants::apexTol * xi * c0 )
        return false; // the cone vertex is reached: the apex return has to take over

      const Vector3d dq_dS = dev / ( 2.0 * q );
      const Vector3d n     = dq_dS + Vector3d::Constant( etaBar / 3.0 );
      const Vector3d df_dS = dq_dS + Vector3d::Constant( eta / 3.0 );

      Vector4d R;
      R.head< 3 >() = state.elasticLogStrainPrincipal - epsTrial + dLambda * n;
      R( 3 )        = ( q + eta * p - xi * cohesion( trialState.alphaP + xi * dLambda ) ) / G;

      if ( R.norm() < Constants::innerNewtonTol ) {

        // a valid cone return keeps the deviatoric direction and a non-negative multiplier
        if ( dLambda < 0.0 || dev.dot( SdevTrial ) <= 0.0 )
          return false;

        result.newMaterialState               = { state.elasticLogStrainPrincipal, trialState.alphaP + xi * dLambda };
        result.newStressState                 = { S };
        result.deltaPlasticLogStrainPrincipal = epsTrial - state.elasticLogStrainPrincipal;
        result.mode                           = ReturnMode::Cone;
        return true;
      }

      const Matrix3d D     = dStress_dElasticLogStrain( state );
      const Matrix3d dn_dS = Pdev / ( 2.0 * q ) - dev * dev.transpose() / ( 4.0 * q * q * q );

      Matrix4d dR;
      dR.topLeftCorner< 3, 3 >()    = Matrix3d::Identity() + dLambda * dn_dS * D;
      dR.topRightCorner< 3, 1 >()   = n;
      dR.bottomLeftCorner< 1, 3 >() = df_dS.transpose() * D / G;
      dR( 3, 3 )                    = -xi * xi * H / G;

      const Vector4d dX = dR.partialPivLu().solve( -R );
      state.elasticLogStrainPrincipal += dX.head< 3 >();
      dLambda += dX( 3 );

      if ( !dX.allFinite() )
        return false;
    }

    return false;
  }

  bool Plasticity::returnToApex( const MaterialState& trialState, ReturnMapResult& result ) const
  {
    const auto& [eta, xi, etaBar, c0, H] = modelParameters;
    const auto& [K, G]                   = materialParameters;

    if ( eta <= 0.0 || etaBar <= 0.0 )
      return false; // there is no apex to return to without friction and dilatancy

    const double thetaTrial = trialState.elasticLogStrainPrincipal.sum();

    double dEpsVol = 0.0;
    for ( int i = 0; i < Constants::nMaxInnerNewtonCycles; i++ ) {

      const double theta = thetaTrial - dEpsVol;
      const double p     = 0.5 * K * std::sinh( 2.0 * theta );
      const double alpha = trialState.alphaP + xi / etaBar * dEpsVol;
      const double r     = ( eta * p - xi * cohesion( alpha ) ) / G;

      if ( std::abs( r ) < Constants::innerNewtonTol ) {
        if ( dEpsVol < 0.0 )
          return false;

        result.newMaterialState               = { Vector3d::Constant( theta / 3.0 ), alpha };
        result.newStressState                 = computeStressState( result.newMaterialState );
        result.deltaPlasticLogStrainPrincipal = trialState.elasticLogStrainPrincipal -
                                                result.newMaterialState.elasticLogStrainPrincipal;
        result.mode = ReturnMode::Apex;
        return true;
      }

      const double dr_dEpsVol = ( -eta * K * std::cosh( 2.0 * theta ) - xi * xi / etaBar * H ) / G;
      dEpsVol -= r / dr_dEpsVol;

      if ( !std::isfinite( dEpsVol ) )
        return false;
    }

    return false;
  }

} // namespace Marmot::Materials
