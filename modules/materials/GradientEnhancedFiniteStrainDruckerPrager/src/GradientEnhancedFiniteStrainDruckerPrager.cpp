#include "Marmot/GradientEnhancedFiniteStrainDruckerPrager.h"
#include "Marmot/MarmotEigenSystems.h"
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotNumericalDifferentiation.h"
#include "Marmot/MarmotStressMeasures.h"
#include "Marmot/MarmotTensorExponential.h"
#include "Marmot/MarmotUtils.h"
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstring>
#include <functional>
#include <stdexcept>

namespace Marmot::Materials {

  using namespace Fastor;
  using namespace FastorStandardTensors;
  namespace Differentiation = NumericalAlgorithms::Differentiation;

  namespace {

    /// absolute tolerance of the return-map residual (the yield function is scaled by the cohesion)
    constexpr double innerNewtonTol = 1e-12;
    /// maximum number of iterations of the local Newton
    constexpr int nMaxInnerNewtonCycles = 50;
    /// maximum number of step halvings of the line search
    constexpr int nMaxHalvings = 10;
    /// a cone solution whose sqrt(J2) is below this fraction of the cohesive strength xi c0 counts as the apex
    constexpr double apexTol = 1e-8;
    /// the apex is admissible if its plastic increment violates the subdifferential condition by less than this
    /// fraction of xi c0, measured in stress (G times the strain violation). It overlaps apexTol, the threshold of
    /// the cone, so that every trial state beyond the cone has a return (no gap between the cone and the apex)
    constexpr double apexAdmissibilityTol = 1e-6;
    /// tolerance of the principal-space return (strains and the scaled yield function; its round-off floor is ~1e-12)
    constexpr double principalNewtonTol = 1e-10;
    /// a principal-space solution is accepted for the full return to the cone if it solves it to this tolerance
    constexpr double principalAcceptanceTol = 1e-9;
    /// relative tolerance of the coaxiality of a cone solution with the trial state
    constexpr double coaxialityTol = 1e-6;
    /// number of material properties without the (optional) density
    constexpr int nRequiredProperties = 10;

    using namespace FastorIndices;

    /// the local Newton and its implicit function theorem work on Eigen vectors (as Marmot's complex step does);
    /// these views map Fastor's row-major tensors to them
    using Matrix9dRowMajor = Eigen::Matrix< double, 9, 9, Eigen::RowMajor >;
    using Vector9d         = Eigen::Matrix< double, 9, 1 >;

    /**
     * @brief Derivative of the trial elastic deformation gradient.
     * @param[in] FpInv Inverse of the plastic deformation gradient at the beginning of the increment.
     * @return @f$ \partial( F_{iK} F^{p,-1}_{KJ} ) / \partial F_{kL} = \delta_{ik} F^{p,-1}_{LJ} @f$.
     */
    Tensor3333d dFeTrial_dF( const Tensor33d& FpInv )
    {
      return einsum< IK, JL, to_IJKL >( Spatial3D::I, Fastor::transpose( FpInv ) );
    }

    /**
     * @brief Derivative of the logarithmic volume ratio.
     * @param[in] F Deformation gradient.
     * @return @f$ \partial \ln\det\boldsymbol{F} / \partial\boldsymbol{F} = \boldsymbol{F}^{-T} @f$.
     */
    Tensor33d dLogDet_dF( const Tensor33d& F )
    {
      return Fastor::transpose( Fastor::inverse( F ) );
    }

    /**
     * @brief Newton's method with the Jacobian by the complex step, and a backtracking line search.
     * @details A step is halved while the residual cannot be evaluated, does not decrease, or leads to an
     * inadmissible iterate.
     * @param[in] residual The residual function @f$ \boldsymbol{R}(\boldsymbol{X}) @f$.
     * @param[in,out] X Initial guess; the solution on success.
     * @param[out] R Residual at @p X on success.
     * @param[out] dR_dX Jacobian at @p X on success.
     * @param[in] admissible Whether an iterate is admissible.
     * @param[in] tolerance Tolerance of the residual norm.
     * @return Whether the iteration converged.
     */
    bool newtonWithBacktracking( const Differentiation::Complex::vector_to_vector_function_type& residual,
                                 Eigen::VectorXd&                                                X,
                                 Eigen::VectorXd&                                                R,
                                 Eigen::MatrixXd&                                                dR_dX,
                                 const std::function< bool( const Eigen::VectorXd& ) >&          admissible,
                                 double tolerance = innerNewtonTol )
    {
      try {
        std::tie( R, dR_dX ) = Differentiation::Complex::forwardDifference( residual, X );
      }
      catch ( const ContinuumMechanics::TensorUtility::TensorExponential::ExponentialMapFailed& ) {
        return false;
      }
      for ( int counter = 0; counter <= nMaxInnerNewtonCycles; counter++ ) {
        if ( !R.allFinite() || !dR_dX.allFinite() )
          return false;
        if ( R.norm() < tolerance )
          return true;

        const Eigen::VectorXd dX   = -dR_dX.colPivHouseholderQr().solve( R );
        double                step = 1.0;
        for ( int halving = 0; halving <= nMaxHalvings; halving++, step *= 0.5 ) {
          const Eigen::VectorXd Xtrial = X + step * dX;
          if ( !admissible( Xtrial ) ) {
            if ( halving == nMaxHalvings )
              return false;
            continue;
          }
          try {
            const auto [Rtrial, Jtrial] = Differentiation::Complex::forwardDifference( residual, Xtrial );
            if ( Rtrial.allFinite() && ( Rtrial.norm() < R.norm() || halving == nMaxHalvings ) ) {
              X     = Xtrial;
              R     = Rtrial;
              dR_dX = Jtrial;
              break;
            }
          }
          catch ( const ContinuumMechanics::TensorUtility::TensorExponential::ExponentialMapFailed& ) {
            // the iterate is too far off for the exponential map: a shorter step
            if ( halving == nMaxHalvings )
              return false;
          }
        }
      }
      return false;
    }

    /**
     * @brief The rotation of the polar decomposition.
     * @param[in] F A tensor with positive determinant.
     * @return @f$ \boldsymbol{R} @f$ of @f$ \boldsymbol{F} = \boldsymbol{R}\boldsymbol{U} @f$.
     */
    Tensor33d rotationOf( const Tensor33d& F )
    {
      const Eigen::JacobiSVD< Eigen::Matrix3d > svd( Eigen::Map< const Eigen::Matrix< double, 3, 3, Eigen::RowMajor > >(
                                                       F.data() ),
                                                     Eigen::ComputeFullU | Eigen::ComputeFullV );
      Tensor33d                                 R;
      Eigen::Map< Eigen::Matrix< double, 3, 3, Eigen::RowMajor > >( R.data() ) = svd.matrixU() *
                                                                                 svd.matrixV().transpose();
      return R;
    }

    /**
     * @brief Principal values of a symmetric tensor.
     * @param[in] symmetric Symmetric tensor.
     * @return Its eigenvalues.
     */
    Fastor::Tensor< double, 3 > principalValues( const Tensor33d& symmetric )
    {
      return Math::computeEigenSystemJacobi( symmetric ).first;
    }

  } // namespace

  GradientEnhancedFiniteStrainDruckerPrager::GradientEnhancedFiniteStrainDruckerPrager(
    const double* materialProperties,
    int           nMaterialProperties,
    int           materialNumber )
    : MarmotMaterialGradientEnhancedFiniteStrain( materialProperties, nMaterialProperties, materialNumber ),
      K( checkedMaterialProperty( materialProperties, nMaterialProperties, 0 ) ),
      G( checkedMaterialProperty( materialProperties, nMaterialProperties, 1 ) ),
      c0( checkedMaterialProperty( materialProperties, nMaterialProperties, 2 ) ),
      frictionAngle( checkedMaterialProperty( materialProperties, nMaterialProperties, 3 ) ),
      dilatancyAngle( checkedMaterialProperty( materialProperties, nMaterialProperties, 4 ) ),
      H( checkedMaterialProperty( materialProperties, nMaterialProperties, 5 ) ),
      softeningModulus( checkedMaterialProperty( materialProperties, nMaterialProperties, 6 ) ),
      maxDamage( checkedMaterialProperty( materialProperties, nMaterialProperties, 7 ) ),
      nonLocalRadius( checkedMaterialProperty( materialProperties, nMaterialProperties, 8 ) ),
      weightingParameter( checkedMaterialProperty( materialProperties, nMaterialProperties, 9 ) ),
      eta( outerConeParameters( frictionAngle ).first ),
      xi( outerConeParameters( frictionAngle ).second ),
      etaBar( outerConeParameters( dilatancyAngle ).first )
  {
    if ( K <= 0.0 || G <= 0.0 )
      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": bulk and shear modulus must be positive" );
    if ( maxDamage < 0.0 || maxDamage >= 1.0 )
      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": expected 0 <= maximum damage < 1" );
    if ( c0 <= 0.0 || softeningModulus <= 0.0 )
      throw std::invalid_argument( MakeString()
                                   << __PRETTY_FUNCTION__ << ": cohesion and softening modulus must be positive" );
    if ( dilatancyAngle > frictionAngle || frictionAngle < 0.0 || dilatancyAngle < 0.0 || frictionAngle >= 90.0 )
      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__
                                                << ": expected 0 <= dilatancy angle <= friction angle < 90 deg" );
    // a softening cohesion would let the cone shrink to its apex; the softening is the damage's
    if ( H < 0.0 )
      throw std::invalid_argument( MakeString()
                                   << __PRETTY_FUNCTION__ << ": the hardening modulus must not be negative" );
    if ( nonLocalRadius <= 0.0 )
      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": the nonlocal radius must be positive" );
    // m > 1 is the over-nonlocal formulation
    if ( weightingParameter < 0.0 )
      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__
                                                << ": the nonlocal weighting parameter must not be negative" );

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
    if ( nMaterialProperties <= nRequiredProperties )
      throw std::runtime_error( MakeString() << __PRETTY_FUNCTION__ << ": density not provided (property 10)" );
    return materialProperties[nRequiredProperties];
  }

  void GradientEnhancedFiniteStrainDruckerPrager::initializeYourself( double* stateVars, int nStateVars )
  {
    for ( int i = 0; i < nStateVars; ++i )
      stateVars[i] = 0.0;

    std::memcpy( stateLayout.getPtr( stateVars, "Fp" ), Spatial3D::I.data(), 9 * sizeof( double ) );
  }

  GradientEnhancedFiniteStrainDruckerPrager::ReturnMapping GradientEnhancedFiniteStrainDruckerPrager::returnMapping(
    const Tensor33d&    F,
    const TensorMap33d& FpOld,
    double              alphaPOld ) const
  {
    // computed once, for the elastic trial and both returns
    const Tensor33d FpOldInv = Fastor::inverse( Tensor33d( FpOld ) );
    const Tensor33d FeTrial  = F % FpOldInv;

    const double JeTrial = Fastor::determinant( FeTrial );
    if ( !( JeTrial > 0.0 ) || !std::isfinite( JeTrial ) )
      throw StressUpdateFailed( MakeString() << __PRETTY_FUNCTION__ << ": non-positive elastic volume ratio" );

    if ( yieldFunction( mandelStress( FeTrial ), alphaPOld ) <= 0.0 ) {
      ReturnMapping elastic;
      elastic.Fe              = FeTrial;
      elastic.FpNew           = FpOld;
      elastic.alphaP          = alphaPOld;
      elastic.dEpVol          = 0.0;
      elastic.plasticWork     = 0.0;
      elastic.dFe_dF          = dFeTrial_dF( FpOldInv );
      elastic.dDeltaAlphaD_dF = Tensor33d( 0.0 );
      return elastic;
    }

    // the return in the principal elastic log strains decides between the cone and the apex by the sign of the
    // (signed) deviatoric radius of its solution, which is smooth across the vertex: no gap between the two
    Tensor33d    FePrincipal;
    double       dLambdaPrincipal = 0.0, rho = 0.0;
    const bool   principal = principalReturn( FeTrial, alphaPOld, FePrincipal, dLambdaPrincipal, rho );
    const double rhoApex   = apexTol * xi * c0 / ( std::sqrt( 2.0 ) * G ); // sqrt(J2) ~ sqrt(2) G rho

    bool converged = false;
    if ( principal && rho > rhoApex ) {
      const auto cone = returnToCone( FeTrial, FpOld, FpOldInv, alphaPOld, converged, &FePrincipal, dLambdaPrincipal );
      if ( converged )
        return cone;
    }

    if ( eta > 0.0 && etaBar > 0.0 ) {
      bool       admissible = false;
      const auto apex       = returnToApex( F, FeTrial, FpOldInv, alphaPOld, admissible );
      if ( admissible )
        return apex;
    }

    // fallbacks: the full return to the cone from the linearized guesses, and, as the last resort, a cone solution
    // with a deviator below the apex threshold (from which the apex cannot be reached)
    const auto cone = returnToCone( FeTrial, FpOld, FpOldInv, alphaPOld, converged );
    if ( converged )
      return cone;
    const auto coneNearVertex = returnToCone( FeTrial,
                                              FpOld,
                                              FpOldInv,
                                              alphaPOld,
                                              converged,
                                              principal && rho > 0.0 ? &FePrincipal : nullptr,
                                              dLambdaPrincipal,
                                              true );
    if ( converged )
      return coneNearVertex;

    // e.g. a tension beyond the apex without dilatancy, or an iterate of the global scheme far off
    throw StressUpdateFailed( MakeString() << __PRETTY_FUNCTION__ << ": no admissible return, neither to the cone "
                                           << "nor to the apex" );
  }

  bool GradientEnhancedFiniteStrainDruckerPrager::principalReturn( const Tensor33d& FeTrial,
                                                                   double           alphaPOld,
                                                                   Tensor33d&       Fe,
                                                                   double&          dLambda,
                                                                   double&          rho ) const
  {
    using namespace Eigen;
    using complexDouble = std::complex< double >;

    // principal directions and elastic log strains of the trial state, in the reference frame of Ce
    const Tensor33d                          CeTrial = Fastor::transpose( FeTrial ) % FeTrial;
    const SelfAdjointEigenSolver< Matrix3d > eig(
      Eigen::Map< const Eigen::Matrix< double, 3, 3, Eigen::RowMajor > >( CeTrial.data() ) );
    if ( eig.info() != Success || !( eig.eigenvalues().minCoeff() > 0.0 ) )
      return false;
    const Vector3d epsTrial = 0.5 * eig.eigenvalues().array().log();
    const Matrix3d Q        = eig.eigenvectors();

    const Vector3d eTrial = epsTrial.array() - epsTrial.mean();
    const Vector3d b1     = deviatoricBasis( 0 );
    const Vector3d b2     = deviatoricBasis( 1 );
    if ( !( eTrial.norm() > 0.0 ) )
      return false; // a hydrostatic trial state has no direction in the deviatoric plane: only the apex

    // initial guess: the radial return of the linearized (Hencky) problem, also beyond the vertex (rho < 0)
    Tensor33d FeTrialPrincipal( 0.0 );
    for ( int a = 0; a < 3; a++ )
      FeTrialPrincipal( a, a ) = std::exp( epsTrial( a ) );
    const Tensor33d MTrial     = mandelStress( FeTrialPrincipal );
    const double    qTrial     = std::sqrt( 0.5 * Fastor::inner( deviatoric( MTrial ), deviatoric( MTrial ) ) );
    const double    dLambdaLin = std::max( 0.0,
                                        yieldFunction( MTrial, alphaPOld ) / ( G + K * eta * etaBar + xi * xi * H ) );

    VectorXd X( 4 );
    X( 0 ) = eTrial.norm() * ( 1.0 - G * dLambdaLin / qTrial );
    X( 1 ) = std::atan2( eTrial.dot( b2 ), eTrial.dot( b1 ) );
    X( 2 ) = 3. * epsTrial.mean() - etaBar * dLambdaLin;
    X( 3 ) = dLambdaLin;

    auto residual = [&]( const VectorXcd& X_ ) -> VectorXcd {
      return principalResidual< complexDouble >( X_, epsTrial, alphaPOld );
    };

    VectorXd R;
    MatrixXd dR_dX;
    if ( !newtonWithBacktracking(
           residual, X, R, dR_dX, []( const VectorXd& X_ ) { return X_( 3 ) >= 0.0; }, principalNewtonTol ) )
      return false;

    rho     = X( 0 );
    dLambda = X( 3 );

    const Vector3d eps = rho * ( std::cos( X( 1 ) ) * b1 + std::sin( X( 1 ) ) * b2 ) + X( 2 ) / 3. * Vector3d::Ones();
    Vector3d       stretchRatio;
    for ( int a = 0; a < 3; a++ )
      stretchRatio( a ) = std::exp( eps( a ) - epsTrial( a ) );
    Tensor33d dFeRatio;
    Eigen::Map< Eigen::Matrix< double, 3, 3, Eigen::RowMajor > >( dFeRatio.data() ) = Q * stretchRatio.asDiagonal() *
                                                                                      Q.transpose();
    Fe = FeTrial % dFeRatio;
    return true;
  }

  GradientEnhancedFiniteStrainDruckerPrager::ReturnMapping GradientEnhancedFiniteStrainDruckerPrager::returnToCone(
    const Tensor33d&    FeTrial,
    const TensorMap33d& FpOld,
    const Tensor33d&    FpOldInv,
    double              alphaPOld,
    bool&               converged,
    const Tensor33d*    FeGuess,
    double              dLambdaGuess,
    bool                acceptVanishingDeviator ) const
  {
    using namespace Eigen;
    using complexDouble = std::complex< double >;

    auto residual = [&]( const VectorXcd& X_ ) -> VectorXcd {
      return coneResidual< complexDouble >( X_, FeTrial, alphaPOld );
    };

    // initial guesses: the return of the linearized (Hencky) problem, exact for small strains, and fractions of it,
    // which keep the iterates away from the singular vertex of the cone
    const Tensor33d MTrial        = mandelStress( FeTrial );
    const Tensor33d devTrial      = deviatoric( MTrial );
    const double    dLambdaLinear = std::max( 0.0,
                                           yieldFunction( MTrial, alphaPOld ) /
                                             ( G + K * eta * etaBar + xi * xi * H ) );

    // an iterate whose deviator has turned against the trial deviator has crossed the vertex of the cone
    const auto admissible = [&]( const VectorXd& X_ ) {
      const Tensor33d dev_ = deviatoric( mandelStress( Tensor33d( X_.head( 9 ).eval().data() ) ) );
      return Fastor::inner( dev_, devTrial ) > 0.0;
    };

    VectorXd X( 11 ), R;
    MatrixXd dR_dX;
    converged = false;

    if ( FeGuess ) {
      // the solution of the principal-space return: the full problem converges at once. Should the Newton not
      // reach its tolerance (the vertex region is ill-conditioned for the full problem), the principal solution is
      // taken if it solves the full problem to a slightly looser tolerance
      X.head( 9 )           = Map< const Matrix< double, 9, 1 > >( FeGuess->data() );
      X( 9 )                = alphaPOld + xi * dLambdaGuess;
      X( 10 )               = dLambdaGuess;
      const VectorXd XGuess = X;
      converged             = admissible( X ) && newtonWithBacktracking( residual, X, R, dR_dX, admissible );
      if ( !converged ) {
        X = XGuess;
        try {
          std::tie( R, dR_dX ) = Differentiation::Complex::forwardDifference( residual, X );
          converged            = R.allFinite() && dR_dX.allFinite() && R.norm() < principalAcceptanceTol;
        }
        catch ( const ContinuumMechanics::TensorUtility::TensorExponential::ExponentialMapFailed& ) {
          converged = false;
        }
      }
    }
    else
      for ( const double fraction : { 1.0, 0.5, 0.25, 0.9, 0.0 } ) {
        const double dLambda0 = fraction * dLambdaLinear;
        Tensor33d    dFp0;
        try {
          dFp0 = ContinuumMechanics::FiniteStrain::Plasticity::FlowIntegration::exponentialMapScalingAndSquaring(
            Tensor33d( dLambda0 * flowDirection( MTrial ) ) );
        }
        catch ( const ContinuumMechanics::TensorUtility::TensorExponential::ExponentialMapFailed& ) {
          continue; // e.g. a vanishing trial deviator: no flow direction for this guess
        }
        const Tensor33d Fe0 = FeTrial % Fastor::inverse( dFp0 );
        X.head( 9 )         = Map< const Matrix< double, 9, 1 > >( Fe0.data() );
        X( 9 )              = alphaPOld + xi * dLambda0;
        X( 10 )             = dLambda0;
        if ( admissible( X ) && newtonWithBacktracking( residual, X, R, dR_dX, admissible ) ) {
          converged = true;
          break;
        }
      }
    if ( !converged )
      return {};

    ReturnMapping r;
    r.plastic               = true;
    r.Fe                    = Tensor33d( X.head( 9 ).eval().data() );
    r.alphaP                = X( 9 );
    const double    dLambda = X( 10 );
    const Tensor33d M       = mandelStress( r.Fe );

    // a solution on the cone needs a positive multiplier and a deviatoric stress. For an isotropic model it is
    // moreover coaxial with the trial state (M and M_trial commute), with a deviator that points the same way as the
    // trial deviator (the finite-strain analogue of sqrt(J2_trial) - G dLambda >= 0); near the apex, Newton may
    // also converge to a spurious, non-coaxial solution. In all these cases the apex takes over.
    const Tensor33d dev        = deviatoric( M );
    const double    normDev    = std::sqrt( Fastor::inner( dev, dev ) );
    const double    normTr     = std::sqrt( Fastor::inner( devTrial, devTrial ) );
    const Tensor33d commutator = dev % devTrial - devTrial % dev;
    // the coaxiality is measured relative to the deviators, with a floor for the round-off of tiny deviators
    const double normM           = std::sqrt( Fastor::inner( M, M ) );
    const double commutatorFloor = 1e-12 * normM * normM;
    const bool   vanishingDev    = std::sqrt( 0.5 ) * normDev <= apexTol * xi * c0;
    if ( dLambda < 0.0 || ( vanishingDev && !acceptVanishingDeviator ) || Fastor::inner( dev, devTrial ) <= 0.0 ||
         std::sqrt( Fastor::inner( commutator, commutator ) ) > coaxialityTol * normDev * normTr + commutatorFloor ) {
      converged = false;
      return {};
    }

    r.FpNew = Fastor::inverse( r.Fe ) % FeTrial % FpOld;

    const Tensor33d dEp = dLambda * flowDirection( M );
    r.dEpVol            = etaBar * dLambda; // tr dg/dM = etaBar
    r.plasticWork       = Fastor::inner( M, dEp );

    // implicit function theorem: R( X, FeTrial ) = 0 with dR/dFeTrial = [-I; 0; 0]
    const Tensor3333d dFeTr_dF                 = dFeTrial_dF( FpOldInv );
    MatrixXd          dR_dF                    = MatrixXd::Zero( 11, 9 );
    dR_dF.topRows( 9 )                         = -Map< const Matrix9dRowMajor >( dFeTr_dF.data() );
    const MatrixXd dX_dF                       = -dR_dX.colPivHouseholderQr().solve( dR_dF );
    Map< Matrix9dRowMajor >( r.dFe_dF.data() ) = dX_dF.topRows( 9 );

    // the local damage increment etaBar dLambda (dLambda >= 0)
    Map< Vector9d >( r.dDeltaAlphaD_dF.data() ) = etaBar * dX_dF.row( 10 ).transpose();

    return r;
  }

  GradientEnhancedFiniteStrainDruckerPrager::ReturnMapping GradientEnhancedFiniteStrainDruckerPrager::returnToApex(
    const Tensor33d& F,
    const Tensor33d& FeTrial,
    const Tensor33d& FpOldInv,
    double           alphaPOld,
    bool&            admissible ) const
  {
    using namespace Eigen;
    using complexDouble = std::complex< double >;

    admissible = false;
    if ( eta <= 0.0 || etaBar <= 0.0 )
      return {}; // there is no apex to return to without friction and dilatancy

    const double thetaTrial = std::log( Fastor::determinant( FeTrial ) );

    auto residual = [&]( const VectorXcd& X_ ) -> VectorXcd {
      return apexResidual< complexDouble >( X_, thetaTrial, alphaPOld );
    };

    VectorXd X( 2 );
    X << thetaTrial, alphaPOld;
    VectorXd   R;
    MatrixXd   dR_dX;
    const bool converged = newtonWithBacktracking( residual, X, R, dR_dX, []( const VectorXd& ) { return true; } );
    if ( !converged || X( 1 ) < alphaPOld )
      return {};

    const double theta = X( 0 );

    // by isotropy, Fp is determined up to a rotation: keep the elastic rotation of the trial state (as the cone return
    // does), Fe = Je^(1/3) R_trial, so that Fp = Fe^-1 F is invariant under a superposed rotation
    const Tensor33d Rtrial = rotationOf( FeTrial );
    ReturnMapping   r;
    r.plastic = true;
    r.Fe      = std::exp( theta / 3. ) * Rtrial;
    r.FpNew   = std::exp( -theta / 3. ) * Tensor33d( Fastor::transpose( Rtrial ) % F );
    r.alphaP  = X( 1 );

    // the plastic log strain increment: the trial elastic log strain minus the new, spherical one
    const auto dEpPrincipal = [&]( const Tensor33d& F_, double theta_ ) {
      const Tensor33d                   FeTrial_ = F_ % FpOldInv;
      const Fastor::Tensor< double, 3 > b = principalValues( Tensor33d( Fastor::transpose( FeTrial_ ) % FeTrial_ ) );
      Fastor::Tensor< double, 3 >       dEp;
      for ( int a = 0; a < 3; a++ )
        dEp( a ) = 0.5 * std::log( b( a ) ) - theta_ / 3.;
      return dEp;
    };
    r.dEpVol = thetaTrial - theta;

    // the apex is admissible only if the plastic increment lies in the subdifferential of g at the vertex:
    // sqrt(2) |dev dEp| <= dEp_v / etaBar (otherwise the state belongs to the cone). The violation is measured in
    // stress, G ( sqrt(2) |dev dEp| - dEp_v / etaBar ), the small-strain sqrt(J2) of the cone solution, against the
    // same scale as the cone's threshold (a relative measure fails for small increments)
    {
      const Fastor::Tensor< double, 3 > dEp    = dEpPrincipal( F, theta );
      const double                      dEpVol = dEp( 0 ) + dEp( 1 ) + dEp( 2 );
      double                            dev2   = 0.0;
      for ( int a = 0; a < 3; a++ )
        dev2 += std::pow( dEp( a ) - dEpVol / 3., 2 );
      const double violation = std::sqrt( 2.0 * dev2 ) - dEpVol / etaBar;
      if ( G * violation > apexAdmissibilityTol * xi * c0 )
        return {};
    }
    admissible     = true;
    const double p = Fastor::trace( mandelStress( r.Fe ) ) / 3.;
    r.plasticWork  = p * ( thetaTrial - theta );

    // implicit function theorem: R( X, thetaTrial ) = 0 with dR/dthetaTrial = [0; -xi/etaBar]
    const Vector2d  dR_dThetaTrial( 0.0, -xi / etaBar );
    const Vector2d  dX_dThetaTrial = -dR_dX.colPivHouseholderQr().solve( dR_dThetaTrial );
    const Tensor33d dTheta_dF      = dX_dThetaTrial( 0 ) * dLogDet_dF( F );

    // Fe = exp( theta / 3 ) R_trial; the rotation does not change the spherical stress, d tau / dFe : dR = 0
    r.dFe_dF = std::exp( theta / 3. ) / 3. * Tensor3333d( outer( Rtrial, dTheta_dF ) );

    // the local damage increment thetaTrial - theta (>= 0 at an admissible apex)
    r.dDeltaAlphaD_dF = dLogDet_dF( F ) - dTheta_dF;

    return r;
  }

  void GradientEnhancedFiniteStrainDruckerPrager::computeStress( ConstitutiveResponse< 3 >& response,
                                                                 AlgorithmicModuli< 3 >&    tangents,
                                                                 const Deformation< 3 >&    deformation,
                                                                 const TimeIncrement&       timeIncrement ) const
  {
    using namespace Eigen;

    double*      sv     = response.stateVars;
    TensorMap33d Fp     = stateLayout.getAs< TensorMap33d >( sv, "Fp" );
    double&      alphaP = stateLayout.getAs< double& >( sv, "alphaP" );
    double&      alphaD = stateLayout.getAs< double& >( sv, "alphaD" );
    double&      kappa  = stateLayout.getAs< double& >( sv, "kappa" );
    double&      omega  = stateLayout.getAs< double& >( sv, "omega" );

    const Tensor33d F( deformation.F );
    const double    N = deformation.N;

    const ReturnMapping r = returnMapping( F, Fp, alphaP );

    std::memcpy( Fp.data(), r.FpNew.data(), 9 * sizeof( double ) );
    alphaP = r.alphaP;

    // the local damage variable grows by the dilatant part of the volumetric plastic strain increment; its
    // derivative is cut off together with it, so that L and dL/dF stay consistent at the kink
    const bool dilatant = r.dEpVol > 0.0;
    if ( dilatant )
      alphaD += r.dEpVol;

    // implicit-gradient damage, irreversible through the history maximum kappa
    const double m             = weightingParameter;
    const double omegaOld      = omega;
    const double kappaOld      = kappa;
    const double alphaWeighted = m * N + ( 1.0 - m ) * alphaD;
    kappa                      = std::max( kappaOld, alphaWeighted );
    const double omegaFree     = kappa > 0.0 ? 1.0 - std::exp( -kappa / softeningModulus ) : 0.0;
    omega                      = std::min( omegaFree, maxDamage );
    const bool   loading       = alphaWeighted > kappaOld && kappa > 0.0 && omegaFree < maxDamage;
    const double dOmega_dKappa = loading ? std::exp( -kappa / softeningModulus ) / softeningModulus : 0.0;

    // effective Kirchhoff stress and its sensitivity to Fe, as in FiniteStrainJ2Plasticity
    using namespace ContinuumMechanics;
    const auto [Ce, dCe_dFe]             = DeformationMeasures::FirstOrderDerived::rightCauchyGreen( r.Fe );
    const auto [psiEff, dPsi_dCe, d2Psi] = EnergyDensityFunctions::SecondOrderDerived::PenceGouPotentialB( Ce, K, G );
    const Tensor33d PK2                  = 2. * dPsi_dCe;
    const auto [tauEff, dTau_dPK2, dTau_dFePartial] = StressMeasures::FirstOrderDerived::KirchhoffStressFromPK2( PK2,
                                                                                                                 r.Fe );
    const Tensor3333d dPK2_dFe                      = einsum< ijKL, KLMN >( Tensor3333d( 2. * d2Psi ), dCe_dFe );
    const Tensor3333d dTauEff_dFe                   = einsum< ijKL, KLMN >( dTau_dPK2, dPK2_dFe ) + dTau_dFePartial;

    response.tau                  = ( 1.0 - omega ) * tauEff;
    response.L                    = alphaD;
    response.nonLocalRadius       = nonLocalRadius;
    response.elasticEnergyDensity = ( 1.0 - omega ) * psiEff;
    // cumulative: the host carries the dissipation of the previous increments in
    response.dissipation += ( 1.0 - omega ) * r.plasticWork + psiEff * ( omega - omegaOld );

    // tangents: tau = (1 - omega) tauEff( Fe( F ) ), omega( kappa ), kappa = m N + (1 - m) alphaD( F )
    const Tensor33d dL_dF     = dilatant ? r.dDeltaAlphaD_dF : Tensor33d( 0.0 );
    const Tensor33d dOmega_dF = dOmega_dKappa * ( 1.0 - m ) * dL_dF;
    const double    dOmega_dN = dOmega_dKappa * m;

    tangents.dTau_dF = ( 1.0 - omega ) * Tensor3333d( einsum< ijKL, KLMN >( dTauEff_dFe, r.dFe_dF ) ) -
                       Tensor3333d( outer( tauEff, dOmega_dF ) );
    tangents.dTau_dN = -dOmega_dN * tauEff;
    tangents.dL_dF   = dL_dF;
    tangents.dL_dN   = 0.0;
  }

} // namespace Marmot::Materials
