#include "Marmot/MarmotMaterialFiniteStrain.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotTypedefs.h"

void applyEigenDeformation_( Fastor::Tensor< double, 3, 3 >& F, double F0_XX, double F0_YY, double F0_ZZ )
{
  F( 0, 0 ) *= F0_XX;
  F( 1, 1 ) *= F0_YY;
  F( 2, 2 ) = F0_ZZ;
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
void MarmotMaterialFiniteStrain::computeStress( ConstitutiveResponse< 3 >&                  response,
                                                AlgorithmicModuli< 3 >&                     tangents,
                                                const Deformation< 3 >&                     deformation,
                                                const TimeIncrement&                        timeIncrement,
                                                double&                                     pNewDT,
                                                const std::tuple< double, double, double >& eigenDeformation )
{
  // TODO: Think about if we simply modify the input increment.

  const auto& [F0_XX, F0_YY, F0_ZZ]    = eigenDeformation;
  auto deformationWithEigenDeformation = deformation;
  applyEigenDeformation_( deformationWithEigenDeformation.F, F0_XX, F0_YY, F0_ZZ );

  computeStress( response, tangents, deformationWithEigenDeformation, timeIncrement, pNewDT );

  applyEigenDeformationToTangent_( tangents.dS_dF, F0_XX, F0_YY, F0_ZZ );

  return;
}

void MarmotMaterialFiniteStrain::computePlaneStrain( ConstitutiveResponse< 3 >& response,
                                                     AlgorithmicModuli< 3 >&    algorithmicModuli,
                                                     const Deformation< 3 >&    deformation,
                                                     const TimeIncrement&       timeIncrement,
                                                     double&                    pNewDT )
{
  return computeStress( response, algorithmicModuli, deformation, timeIncrement, pNewDT );
}

void MarmotMaterialFiniteStrain::computePlaneStrain( ConstitutiveResponse< 3 >&                  response,
                                                     AlgorithmicModuli< 3 >&                     algorithmicModuli,
                                                     const Deformation< 3 >&                     deformation,
                                                     const TimeIncrement&                        timeIncrement,
                                                     double&                                     pNewDT,
                                                     const std::tuple< double, double, double >& eigenDeformation )
{
  return computeStress( response, algorithmicModuli, deformation, timeIncrement, pNewDT, eigenDeformation );
}

void MarmotMaterialFiniteStrain::computePlaneStress( ConstitutiveResponse< 2 >& response,
                                                     AlgorithmicModuli< 2 >&    algorithmicModuli,
                                                     const Deformation< 2 >&    deformation,
                                                     const TimeIncrement&       timeIncrement,
                                                     double&                    pNewDT )
{
  throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << "Not yet implemented." );
}

/** The basic implementation here assumes non-chirial, isotropic elastic behavior.
 */
std::tuple< double, double, double > MarmotMaterialFiniteStrain::findEigenDeformationForEigenStress(
  const std::tuple< double, double, double >& initialGuess,
  const std::tuple< double, double, double >& eigenStressComponents )
{
  using namespace Marmot;

  Deformation< 3 >          deformation = { 0.0 };
  ConstitutiveResponse< 3 > response    = { 0.0, 0.0, 0.0 };
  AlgorithmicModuli< 3 >    tangents;

  const double        time[] = { 0.0, 0.0 };
  const TimeIncrement timeIncrement{ time[0], 0.0 };
  double              pNewDT = 1e36;

  Eigen::Map< Eigen::VectorXd > theStateVars( stateVars, nStateVars );

  auto evaluateStress = [&]( const Vector3d& F0 ) {
    deformation.F( 0, 0 ) = F0( 0 );
    deformation.F( 1, 1 ) = F0( 1 );
    deformation.F( 2, 2 ) = F0( 2 );

    computeStress( response, tangents, deformation, timeIncrement, pNewDT );

    Vector3d S = { response.S( 0, 0 ), response.S( 1, 1 ), response.S( 2, 2 ) };
    Matrix3d dS_dF;
    dS_dF << tangents.dS_dF( 0, 0, 0, 0 ), tangents.dS_dF( 0, 0, 1, 1 ), tangents.dS_dF( 0, 0, 2, 2 ),
      tangents.dS_dF( 1, 1, 0, 0 ), tangents.dS_dF( 1, 1, 1, 1 ), tangents.dS_dF( 1, 1, 2, 2 ),
      tangents.dS_dF( 2, 2, 0, 0 ), tangents.dS_dF( 2, 2, 1, 1 ), tangents.dS_dF( 2, 2, 2, 2 );

    return std::tuple< Eigen::Vector3d, Eigen::Matrix3d >( S, dS_dF );
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
