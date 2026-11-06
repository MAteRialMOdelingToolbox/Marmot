#include "Marmot/MarmotMaterialFiniteStrain.h"
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotTypedefs.h"

void applyEigenDeformation_( Fastor::Tensor< double, 3, 3 >& F, double F0_XX, double F0_YY, double F0_ZZ )
{
  F( 0, 0 ) *= F0_XX;
  F( 1, 1 ) *= F0_YY;
  F( 2, 2 ) *= F0_ZZ;
}

void applyEigenDeformationToTangent_( Fastor::Tensor< double, 3, 3 >& dT0_dF, double F0_XX, double F0_YY, double F0_ZZ )
{
  dT0_dF( 0, 0 ) *= F0_XX;
  dT0_dF( 1, 1 ) *= F0_YY;
  dT0_dF( 2, 2 ) *= F0_ZZ;
}

void applyEigenDeformationToTangent_( Fastor::Tensor< double, 3, 3, 3, 3 >& dT2_dF,
                                      double                                F0_XX,
                                      double                                F0_YY,
                                      double                                F0_ZZ )
{
  for ( int i = 0; i < 3; i++ ) {
    for ( int j = 0; j < 3; j++ ) {
      dT2_dF( i, j, 0, 0 ) *= F0_XX;
      dT2_dF( i, j, 1, 1 ) *= F0_YY;
      dT2_dF( i, j, 2, 2 ) *= F0_ZZ;
    }
  }
}
void Marmot::Materials::MaterialFiniteStrain::computeStress(
  ConstitutiveResponse< 3 >&                  response,
  AlgorithmicModuli< 3 >&                     tangents,
  const Deformation< 3 >&                     deformation,
  const TimeIncrement&                        timeIncrement,
  const std::tuple< double, double, double >& eigenDeformation ) const
{
  // TODO: Think about if we simply modify the input increment.

  const auto& [F0_XX, F0_YY, F0_ZZ]    = eigenDeformation;
  auto deformationWithEigenDeformation = deformation;
  applyEigenDeformation_( deformationWithEigenDeformation.F, F0_XX, F0_YY, F0_ZZ );

  computeStress( response, tangents, deformationWithEigenDeformation, timeIncrement );

  applyEigenDeformationToTangent_( tangents.dTau_dF, F0_XX, F0_YY, F0_ZZ );

  return;
}

void Marmot::Materials::MaterialFiniteStrain::computePlaneStrain( ConstitutiveResponse< 3 >& response,
                                                                  AlgorithmicModuli< 3 >&    algorithmicModuli,
                                                                  const Deformation< 3 >&    deformation,
                                                                  const TimeIncrement&       timeIncrement ) const
{
  return computeStress( response, algorithmicModuli, deformation, timeIncrement );
}

void Marmot::Materials::MaterialFiniteStrain::computePlaneStrain(
  ConstitutiveResponse< 3 >&                  response,
  AlgorithmicModuli< 3 >&                     algorithmicModuli,
  const Deformation< 3 >&                     deformation,
  const TimeIncrement&                        timeIncrement,
  const std::tuple< double, double, double >& eigenDeformation ) const
{
  return computeStress( response, algorithmicModuli, deformation, timeIncrement, eigenDeformation );
}

void Marmot::Materials::MaterialFiniteStrain::computePlaneStress( ConstitutiveResponse< 2 >& response,
                                                                  AlgorithmicModuli< 2 >&    algorithmicModuli,
                                                                  const Deformation< 2 >&    deformation,
                                                                  const TimeIncrement&       timeIncrement ) const
{
  throw std::invalid_argument( Marmot::Utilities::MakeString() << __PRETTY_FUNCTION__ << "Not yet implemented." );
}

/** The basic implementation here assumes non-chirial, isotropic elastic behavior.
 */
void Marmot::Materials::MaterialFiniteStrain::computeUniaxialStress( ConstitutiveResponse< 3 >& response,
                                                                     double&                    dTau11_dF11,
                                                                     Deformation< 3 >&          deformation,
                                                                     const TimeIncrement&       timeIncrement ) const
{
  using Marmot::Utilities::MakeString;
  using Marmot::Utilities::StressUpdateFailed;

  constexpr int    maxIterations = 25;
  constexpr double tolerance     = 1e-12; // on the Newton correction of the lateral stretches, i.e., a strain

  Eigen::Map< Eigen::VectorXd > stateVars( response.stateVars, stateLayout.totalSize() );
  const Eigen::VectorXd         stateVarsOld = stateVars;
  const auto                    responseOld  = response;

  AlgorithmicModuli< 3 > tangents;
  Eigen::Matrix2d        dTauLateral_dFLateral;

  bool finalEvaluation = false;
  for ( int iteration = 0;; iteration++ ) {

    // every evaluation starts from the beginning of the increment, a path dependent material must not see the
    // trial states of this iteration
    stateVars = stateVarsOld;
    response  = responseOld;

    computeStress( response, tangents, deformation, timeIncrement );

    const auto& C = tangents.dTau_dF;
    dTauLateral_dFLateral << C( 1, 1, 1, 1 ), C( 1, 1, 2, 2 ), C( 2, 2, 1, 1 ), C( 2, 2, 2, 2 );

    if ( finalEvaluation )
      break;

    const Eigen::Vector2d residual( response.tau( 1, 1 ), response.tau( 2, 2 ) );

    const auto            solver     = dTauLateral_dFLateral.fullPivLu();
    const Eigen::Vector2d correction = -solver.solve( residual );

    if ( !solver.isInvertible() || !correction.allFinite() )
      throw StressUpdateFailed( MakeString() << __PRETTY_FUNCTION__
                                             << ": singular lateral tangent in the uniaxial stress iteration" );

    if ( correction.norm() <= tolerance ) {
      // apply the last (tiny) correction and evaluate once more: with the quadratic convergence of Newton's method,
      // the remaining lateral stress is then at round-off level; returning the state before the correction would
      // leave lateral stresses of (stiffness x tolerance), i.e., noise in the axial stress
      deformation.F( 1, 1 ) += correction( 0 );
      deformation.F( 2, 2 ) += correction( 1 );
      finalEvaluation = true;
      continue;
    }

    if ( iteration >= maxIterations ) {
      Marmot::Utilities::Journal::warningToMSG( "UniaxialStress (finite strain) requires cutback" );
      throw StressUpdateFailed( MakeString() << __PRETTY_FUNCTION__ << ": uniaxial stress iteration did not converge" );
    }

    // a lateral stretch must stay positive: damp the correction if it would invert the lateral fibers
    double damping = 1.0;
    while ( deformation.F( 1, 1 ) + damping * correction( 0 ) <= 0.0 ||
            deformation.F( 2, 2 ) + damping * correction( 1 ) <= 0.0 )
      damping *= 0.5;

    deformation.F( 1, 1 ) += damping * correction( 0 );
    deformation.F( 2, 2 ) += damping * correction( 1 );
  }

  // condensation of the lateral stretches: dTau11/dF11 under tau22 = tau33 = 0
  const auto&           C = tangents.dTau_dF;
  const Eigen::Vector2d dTau11_dFLateral( C( 0, 0, 1, 1 ), C( 0, 0, 2, 2 ) );
  const Eigen::Vector2d dTauLateral_dF11( C( 1, 1, 0, 0 ), C( 2, 2, 0, 0 ) );

  dTau11_dF11 = C( 0, 0, 0, 0 ) - dTau11_dFLateral.dot( dTauLateral_dFLateral.fullPivLu().solve( dTauLateral_dF11 ) );
}

std::tuple< double, double, double > Marmot::Materials::MaterialFiniteStrain::findEigenDeformationForEigenStress(
  const std::tuple< double, double, double >& initialGuess,
  const std::tuple< double, double, double >& eigenStress,
  double*                                     stateVars ) const
{
  using namespace Marmot;

  Deformation< 3 >          deformation = { 0.0 };
  ConstitutiveResponse< 3 > response;
  response.stateVars = stateVars;
  AlgorithmicModuli< 3 > tangents;

  const double        time[] = { 0.0, 0.0 };
  const TimeIncrement timeIncrement{ time[0], 0.0 };

  Eigen::Map< Eigen::VectorXd > theStateVars( stateVars, stateLayout.totalSize() );

  auto evaluateStress = [&]( const Marmot::Math::Vector3d& F0 ) {
    deformation.F( 0, 0 ) = F0( 0 );
    deformation.F( 1, 1 ) = F0( 1 );
    deformation.F( 2, 2 ) = F0( 2 );

    computeStress( response, tangents, deformation, timeIncrement );

    Marmot::Math::Vector3d tau = { response.tau( 0, 0 ), response.tau( 1, 1 ), response.tau( 2, 2 ) };
    Marmot::Math::Matrix3d dTau_dF;
    dTau_dF << tangents.dTau_dF( 0, 0, 0, 0 ), tangents.dTau_dF( 0, 0, 1, 1 ), tangents.dTau_dF( 0, 0, 2, 2 ),
      tangents.dTau_dF( 1, 1, 0, 0 ), tangents.dTau_dF( 1, 1, 1, 1 ), tangents.dTau_dF( 1, 1, 2, 2 ),
      tangents.dTau_dF( 2, 2, 0, 0 ), tangents.dTau_dF( 2, 2, 1, 1 ), tangents.dTau_dF( 2, 2, 2, 2 );

    return std::tuple< Eigen::Vector3d, Eigen::Matrix3d >( tau, dTau_dF );
  };

  const auto& [F0_XX, F0_YY, F0_ZZ] = initialGuess;
  const auto& [S_XX, S_YY, S_ZZ]    = eigenStress;
  Eigen::Vector3d def               = { F0_XX, F0_YY, F0_ZZ };
  Eigen::Vector3d eigenNormalStress = { S_XX, S_YY, S_ZZ };
  Eigen::Vector3d R, dF;

  Eigen::VectorXd materialStateVarsBackup = theStateVars;

  int itCounter = 0;
  while ( true ) {
    auto [normalStress, dNormalStress_dF] = evaluateStress( def );
    theStateVars                          = materialStateVarsBackup;

    R = normalStress - eigenNormalStress;

    if ( R.norm() / std::min( normalStress.norm(), 1.0 ) <= 1e-10 )
      break;

    if ( itCounter > 5 )
      throw std::invalid_argument( Marmot::Utilities::MakeString()
                                   << __PRETTY_FUNCTION__ << ": failed to find eigen deformation within 5 iterations" );

    itCounter++;

    dF = -dNormalStress_dF.inverse() * R;

    def += dF;
  }

  return { def( 0 ), def( 1 ), def( 2 ) };
}
