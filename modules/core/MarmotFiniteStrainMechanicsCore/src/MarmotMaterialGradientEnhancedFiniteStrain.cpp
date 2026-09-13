#include "Marmot/MarmotMaterialGradientEnhancedFiniteStrain.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotTypedefs.h"

namespace {

  /// @brief Scale the diagonal of a deformation gradient by the eigen deformation.
  void applyEigenDeformation( Fastor::Tensor< double, 3, 3 >& F, double F0_XX, double F0_YY, double F0_ZZ )
  {
    F( 0, 0 ) *= F0_XX;
    F( 1, 1 ) *= F0_YY;
    F( 2, 2 ) *= F0_ZZ;
  }

  /// @brief Scale a first-order derivative with respect to the deformation gradient by the eigen deformation.
  void applyEigenDeformationToTangent( Fastor::Tensor< double, 3, 3 >& dT0_dF,
                                       double                          F0_XX,
                                       double                          F0_YY,
                                       double                          F0_ZZ )
  {
    dT0_dF( 0, 0 ) *= F0_XX;
    dT0_dF( 1, 1 ) *= F0_YY;
    dT0_dF( 2, 2 ) *= F0_ZZ;
  }

  /// @brief Scale a second-order derivative with respect to the deformation gradient by the eigen deformation.
  void applyEigenDeformationToTangent( Fastor::Tensor< double, 3, 3, 3, 3 >& dT2_dF,
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

} // namespace

void MarmotMaterialGradientEnhancedFiniteStrain::computeStress(
  ConstitutiveResponse< 3 >&                  response,
  AlgorithmicModuli< 3 >&                     tangents,
  const Deformation< 3 >&                     deformation,
  const TimeIncrement&                        timeIncrement,
  const std::tuple< double, double, double >& eigenDeformation ) const
{
  const auto& [F0_XX, F0_YY, F0_ZZ]    = eigenDeformation;
  auto deformationWithEigenDeformation = deformation;
  applyEigenDeformation( deformationWithEigenDeformation.F, F0_XX, F0_YY, F0_ZZ );

  computeStress( response, tangents, deformationWithEigenDeformation, timeIncrement );

  // Both tangents differentiate with respect to F and therefore pick up the eigen deformation by the
  // chain rule; dTau_dN does not, since the nonlocal field is not deformed.
  applyEigenDeformationToTangent( tangents.dTau_dF, F0_XX, F0_YY, F0_ZZ );
  applyEigenDeformationToTangent( tangents.dL_dF, F0_XX, F0_YY, F0_ZZ );
}

void MarmotMaterialGradientEnhancedFiniteStrain::computePlaneStrain( ConstitutiveResponse< 3 >& response,
                                                                     AlgorithmicModuli< 3 >&    algorithmicModuli,
                                                                     const Deformation< 3 >&    deformation,
                                                                     const TimeIncrement&       timeIncrement ) const
{
  return computeStress( response, algorithmicModuli, deformation, timeIncrement );
}

void MarmotMaterialGradientEnhancedFiniteStrain::computePlaneStrain(
  ConstitutiveResponse< 3 >&                  response,
  AlgorithmicModuli< 3 >&                     algorithmicModuli,
  const Deformation< 3 >&                     deformation,
  const TimeIncrement&                        timeIncrement,
  const std::tuple< double, double, double >& eigenDeformation ) const
{
  return computeStress( response, algorithmicModuli, deformation, timeIncrement, eigenDeformation );
}

/** The basic implementation here assumes non-chiral, isotropic elastic behavior, and an undamaged
 * initial state: the nonlocal field is held at zero throughout the iteration.
 */
std::tuple< double, double, double > MarmotMaterialGradientEnhancedFiniteStrain::findEigenDeformationForEigenStress(
  const std::tuple< double, double, double >& initialGuess,
  const std::tuple< double, double, double >& eigenStressComponents,
  double*                                     stateVars ) const
{
  using namespace Marmot;

  Deformation< 3 >          deformation = { Fastor::Tensor< double, 3, 3 >( 0.0 ), 0.0 };
  ConstitutiveResponse< 3 > response;
  response.stateVars = stateVars;
  AlgorithmicModuli< 3 > tangents;

  const TimeIncrement timeIncrement{ 0.0, 0.0 };

  Eigen::Map< Eigen::VectorXd > theStateVars( stateVars, stateLayout.totalSize() );

  auto evaluateStress = [&]( const Vector3d& F0 ) {
    deformation.F( 0, 0 ) = F0( 0 );
    deformation.F( 1, 1 ) = F0( 1 );
    deformation.F( 2, 2 ) = F0( 2 );

    computeStress( response, tangents, deformation, timeIncrement );

    Vector3d tau = { response.tau( 0, 0 ), response.tau( 1, 1 ), response.tau( 2, 2 ) };
    Matrix3d dTau_dF;
    dTau_dF << tangents.dTau_dF( 0, 0, 0, 0 ), tangents.dTau_dF( 0, 0, 1, 1 ), tangents.dTau_dF( 0, 0, 2, 2 ),
      tangents.dTau_dF( 1, 1, 0, 0 ), tangents.dTau_dF( 1, 1, 1, 1 ), tangents.dTau_dF( 1, 1, 2, 2 ),
      tangents.dTau_dF( 2, 2, 0, 0 ), tangents.dTau_dF( 2, 2, 1, 1 ), tangents.dTau_dF( 2, 2, 2, 2 );

    return std::tuple< Eigen::Vector3d, Eigen::Matrix3d >( tau, dTau_dF );
  };

  const auto& [F0_XX, F0_YY, F0_ZZ] = initialGuess;
  const auto& [S_XX, S_YY, S_ZZ]    = eigenStressComponents;
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
      throw std::invalid_argument( MakeString()
                                   << __PRETTY_FUNCTION__ << ": failed to find eigen deformation within 5 iterations" );

    itCounter++;

    dF = -dNormalStress_dF.inverse() * R;

    def += dF;
  }

  return { def( 0 ), def( 1 ), def( 2 ) };
}
