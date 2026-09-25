#include "Marmot/GradientEnhancedFiniteStrainDruckerPrager.h"
#include "Marmot/GradientEnhancedFiniteStrainDruckerPragerConstants.h"
#include "Marmot/GradientEnhancedFiniteStrainDruckerPragerDamage.h"
#include "Marmot/GradientEnhancedFiniteStrainDruckerPragerPlasticity.h"
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotUtils.h"
#include <Eigen/Dense>
#include <cmath>
#include <cstring>
#include <stdexcept>
#include <vector>

namespace Marmot::Materials {

  using Tensor33d     = Fastor::Tensor< double, 3, 3 >;
  namespace Constants = GradientEnhancedFiniteStrainDruckerPragerConstants;

  namespace {

    Eigen::Matrix3d toEigen( const Tensor33d& T )
    {
      Eigen::Matrix3d M;
      for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ )
          M( i, j ) = T( i, j );
      return M;
    }

    Tensor33d toFastor( const Eigen::Matrix3d& M )
    {
      Tensor33d T;
      for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ )
          T( i, j ) = M( i, j );
      return T;
    }

  } // namespace

  GradientEnhancedFiniteStrainDruckerPrager::GradientEnhancedFiniteStrainDruckerPrager(
    const double* materialProperties,
    int           nMaterialProperties,
    int           materialNumber )
    : MarmotMaterialGradientEnhancedFiniteStrain( materialProperties, nMaterialProperties, materialNumber ),
      K( materialProperties[0] ),
      G( materialProperties[1] ),
      c0( materialProperties[2] ),
      frictionAngle( materialProperties[3] ),
      dilatancyAngle( materialProperties[4] ),
      H( materialProperties[5] ),
      As( materialProperties[6] ),
      softeningModulus( materialProperties[7] ),
      maxDamage( materialProperties[8] ),
      nonLocalRadius( materialProperties[9] ),
      weightingParameter( materialProperties[10] )
  {
    if ( nMaterialProperties < 11 )
      throw std::invalid_argument( MakeString()
                                   << __PRETTY_FUNCTION__ << ": expected at least 11 material properties, got "
                                   << nMaterialProperties );
    if ( c0 <= 0.0 || softeningModulus <= 0.0 )
      throw std::invalid_argument( MakeString()
                                   << __PRETTY_FUNCTION__ << ": cohesion and softening modulus must be positive" );
    if ( dilatancyAngle > frictionAngle || frictionAngle < 0.0 || dilatancyAngle < 0.0 || frictionAngle >= 90.0 )
      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__
                                                << ": expected 0 <= dilatancy angle <= friction angle < 90 deg" );

    stateLayout.add( "Fp", 9 );     // plastic deformation gradient
    stateLayout.add( "alphaP", 1 ); // plastic hardening variable
    stateLayout.add( "alphaD", 1 ); // local damage variable, the source L of the nonlocal balance
    stateLayout.add( "kappa", 1 );  // damage history
    stateLayout.add( "omega", 1 );  // scalar damage
    stateLayout.finalize();
  }

  std::pair< double, double > GradientEnhancedFiniteStrainDruckerPrager::outerConeParameters( double angle )
  {
    const double s = std::sin( Math::degToRad( angle ) );
    const double c = std::cos( Math::degToRad( angle ) );
    const double d = std::sqrt( 3.0 ) * ( 3.0 - s );
    return { 6.0 * s / d, 6.0 * c / d };
  }

  double GradientEnhancedFiniteStrainDruckerPrager::getDensity( const double* stateVars ) const
  {
    if ( nMaterialProperties <= 11 )
      throw std::runtime_error( MakeString() << __PRETTY_FUNCTION__ << ": density not provided (property 11)" );
    return materialProperties[11];
  }

  void GradientEnhancedFiniteStrainDruckerPrager::initializeYourself( double* stateVars, int nStateVars )
  {
    for ( int i = 0; i < nStateVars; ++i )
      stateVars[i] = 0.0;

    const Tensor33d I = toFastor( Eigen::Matrix3d::Identity() );
    std::memcpy( stateLayout.getPtr( stateVars, "Fp" ), I.data(), 9 * sizeof( double ) );
  }

  std::tuple< Tensor33d, double, double > GradientEnhancedFiniteStrainDruckerPrager::stressUpdate( const Tensor33d& F_,
                                                                                                   double nonLocalField,
                                                                                                   double* sv ) const
  {
    using Plasticity = GradientEnhancedFiniteStrainDruckerPragerPlasticity;

    const auto [eta, xi] = outerConeParameters( frictionAngle );
    const double etaBar  = outerConeParameters( dilatancyAngle ).first;

    const Plasticity                                      plasticity( { K, G }, { eta, xi, etaBar, c0, H } );
    const GradientEnhancedFiniteStrainDruckerPragerDamage damageLaw(
      { As, softeningModulus, weightingParameter, maxDamage } );

    // read the plastic deformation gradient through the very same Tensor33d( ptr ) construction that writes it,
    // so that the two can never disagree about the storage order of the flat state array
    const Eigen::Matrix3d Fp     = toEigen( Tensor33d( stateLayout.getPtr( sv, "Fp" ) ) );
    double&               alphaP = stateLayout.getAs< double& >( sv, "alphaP" );
    double&               alphaD = stateLayout.getAs< double& >( sv, "alphaD" );
    double&               kappa  = stateLayout.getAs< double& >( sv, "kappa" );
    double&               omega  = stateLayout.getAs< double& >( sv, "omega" );

    const Eigen::Matrix3d F = toEigen( F_ );

    // trial elastic state, decomposed spectrally through the elastic left Cauchy-Green tensor
    const Eigen::Matrix3d                            FeTrial = F * Fp.inverse();
    Eigen::SelfAdjointEigenSolver< Eigen::Matrix3d > es( FeTrial * FeTrial.transpose() );
    const Eigen::Vector3d                            b = es.eigenvalues();
    const Eigen::Matrix3d                            n = es.eigenvectors();
    if ( b.minCoeff() <= 0.0 || !b.allFinite() )
      throw StressUpdateFailed( MakeString() << __PRETTY_FUNCTION__ << ": non-positive elastic stretch" );

    const Plasticity::MaterialState trial{ 0.5 * b.array().log().matrix(), alphaP };

    Plasticity::ReturnMapResult result;
    try {
      result = plasticity.performReturnMapping( trial );
    }
    catch ( const Plasticity::ReturnMappingFailedException& ) {
      throw StressUpdateFailed( MakeString() << __PRETTY_FUNCTION__ << ": return mapping not successful" );
    }

    // exponential-map update of the plastic deformation gradient, elastic rotation frozen:
    // Fp_new = FeTrial^-1 exp( dEp ) F, i.e. F Fp_new^-1 = V_e,new R_e
    const Eigen::Vector3d& dEp    = result.deltaPlasticLogStrainPrincipal;
    const Eigen::Matrix3d  expDEp = n * dEp.array().exp().matrix().asDiagonal() * n.transpose();
    const Eigen::Matrix3d  FpNew  = FeTrial.inverse() * expDEp * F;

    const Tensor33d FpNew_ = toFastor( FpNew );
    std::memcpy( stateLayout.getPtr( sv, "Fp" ), FpNew_.data(), 9 * sizeof( double ) );
    alphaP = result.newMaterialState.alphaP;

    // implicit-gradient damage
    alphaD += damageLaw.deltaAlphaLocal( dEp );
    const auto damage = damageLaw.computeDamage( alphaD, nonLocalField, kappa );
    kappa             = damage.kappa;
    omega             = damage.omega;

    // Kirchhoff stress: the principal Mandel stresses on the principal directions of the elastic stretch
    const Eigen::Vector3d& S      = result.newStressState.mandelPrincipal;
    const Eigen::Matrix3d  tauEff = n * S.asDiagonal() * n.transpose();

    // elastic energy density of PenceGouPotentialB in the principal elastic log strains
    const Eigen::Vector3d& epsE   = result.newMaterialState.elasticLogStrainPrincipal;
    const double           theta  = epsE.sum();
    const double           J2     = std::exp( 2.0 * theta );
    const double           psiEff = K / 8. * ( J2 + 1. / J2 - 2. ) +
                          G / 2. * ( ( 2.0 * ( epsE.array() - theta / 3.0 ) ).exp().sum() - 3. );

    return { toFastor( ( 1.0 - omega ) * tauEff ), alphaD, ( 1.0 - omega ) * psiEff };
  }

  void GradientEnhancedFiniteStrainDruckerPrager::computeStress( ConstitutiveResponse< 3 >& response,
                                                                 AlgorithmicModuli< 3 >&    tangents,
                                                                 const Deformation< 3 >&    deformation,
                                                                 const TimeIncrement&       timeIncrement ) const
  {
    const int                   nsv = getNumberOfRequiredStateVars();
    const std::vector< double > oldState( response.stateVars, response.stateVars + nsv );

    const Tensor33d F( deformation.F );
    const double    N = deformation.N;

    const auto [tau, L, psi] = stressUpdate( F, N, response.stateVars );

    response.tau                  = tau;
    response.L                    = L;
    response.nonLocalRadius       = nonLocalRadius;
    response.elasticEnergyDensity = psi;
    response.dissipation          = 0.0;

    // algorithmic tangents by forward finite differences of the full state update; every probe starts from the
    // state at the beginning of the increment
    const double eps = Constants::tangentPerturbation;
    for ( int k = 0; k < 3; ++k )
      for ( int l = 0; l < 3; ++l ) {
        Tensor33d Fpert = F;
        Fpert( k, l ) += eps;
        std::vector< double > sv             = oldState;
        const auto [tauPert, LPert, psiPert] = stressUpdate( Fpert, N, sv.data() );
        for ( int i = 0; i < 3; ++i )
          for ( int j = 0; j < 3; ++j )
            tangents.dTau_dF( i, j, k, l ) = ( tauPert( i, j ) - tau( i, j ) ) / eps;
        tangents.dL_dF( k, l ) = ( LPert - L ) / eps;
      }

    {
      const double          epsN           = Constants::tangentPerturbationNonLocal;
      std::vector< double > sv             = oldState;
      const auto [tauPert, LPert, psiPert] = stressUpdate( F, N + epsN, sv.data() );
      for ( int i = 0; i < 3; ++i )
        for ( int j = 0; j < 3; ++j )
          tangents.dTau_dN( i, j ) = ( tauPert( i, j ) - tau( i, j ) ) / epsN;
      tangents.dL_dN = ( LPert - L ) / epsN;
    }
  }

} // namespace Marmot::Materials
