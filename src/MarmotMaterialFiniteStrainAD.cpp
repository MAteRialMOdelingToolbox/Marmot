#include "Marmot/MarmotMaterialFiniteStrainAD.h"
#include "Marmot/MarmotAutomaticDifferentiationForFastor.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include "Marmot/MarmotMicromorphicTensorBasics.h"
#include "Marmot/MarmotTypedefs.h"
#include "autodiff/forward/dual.hpp"

void applyEigenDeformation_( Fastor::Tensor< autodiff::dual, 3, 3 >& F, double F0_XX, double F0_YY, double F0_ZZ )
{
  F( 0, 0 ) *= F0_XX;
  F( 1, 1 ) *= F0_YY;
  F( 2, 2 ) *= F0_ZZ;
}

void MarmotMaterialFiniteStrainAD::computeStress( ConstitutiveResponse< 3 >&                  response,
                                                  const Deformation< 3 >&                     deformation,
                                                  const TimeIncrement&                        timeIncrement,
                                                  const std::tuple< double, double, double >& eigenDeformation )
{
  // TODO: Think about if we simply modify the input increment.

  const auto& [F0_XX, F0_YY, F0_ZZ]    = eigenDeformation;
  auto deformationWithEigenDeformation = deformation;
  applyEigenDeformation_( deformationWithEigenDeformation.F, F0_XX, F0_YY, F0_ZZ );

  computeStress( response, deformationWithEigenDeformation, timeIncrement );

  // applyEigenDeformationToTangent_( tangents.dS_dF, F0_XX, F0_YY, F0_ZZ );

  return;
}

void MarmotMaterialFiniteStrainAD::computePlaneStrain( ConstitutiveResponse< 3 >& response,
                                                       const Deformation< 3 >&    deformation,
                                                       const TimeIncrement&       timeIncrement )
{
  return computeStress( response, deformation, timeIncrement );
}

void MarmotMaterialFiniteStrainAD::computePlaneStrain( ConstitutiveResponse< 3 >&                  response,
                                                       const Deformation< 3 >&                     deformation,
                                                       const TimeIncrement&                        timeIncrement,
                                                       const std::tuple< double, double, double >& eigenDeformation )
{
  return computeStress( response, deformation, timeIncrement, eigenDeformation );
}

void MarmotMaterialFiniteStrainAD::computePlaneStress( ConstitutiveResponse< 2 >& response,
                                                       const Deformation< 2 >&    deformation,
                                                       const TimeIncrement&       timeIncrement )
{
  throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << "Not yet implemented." );
}

/** The basic implementation here assumes non-chirial, isotropic elastic behavior.
 */
std::tuple< double, double, double > MarmotMaterialFiniteStrainAD::findEigenDeformationForEigenStress(
  const std::tuple< double, double, double >& initialGuess,
  const std::tuple< double, double, double >& eigenStressComponents )
{
  using namespace Marmot;
  using namespace FastorStandardTensors;

  Tensor33d deformation( 0.0 );
  Tensor33d stress( 0.0 );

  const double        time[] = { 0.0, 0.0 };
  const TimeIncrement timeIncrement{ time[0], 0.0 };

  Eigen::Map< Eigen::VectorXd > theStateVars( stateVars, nStateVars );

  auto evaluateStress = [&]( const Vector3d& F0 ) {
    deformation( 0, 0 ) = F0( 0 );
    deformation( 1, 1 ) = F0( 1 );
    deformation( 2, 2 ) = F0( 2 );

    /* computeStress( response, deformation, timeIncrement ); */

    // compute tangent using autodiff
    using func_type = std::function< Tensor33t< autodiff::dual >( const Tensor33t< autodiff::dual >& ) >;
    func_type func  = [&]( const Tensor33t< autodiff::dual >& F ) {
      Deformation< 3 >          def  = { F };
      ConstitutiveResponse< 3 > resp = { 0.0, 0.0, 0.0 };
      computeStress( resp, def, timeIncrement );
      return resp.tau;
    };

    Tensor3333d dS_dF_;

    std::tie( stress, dS_dF_ ) = Marmot::AutomaticDifferentiation::dF_dT( func, deformation );

    Vector3d S = { stress( 0, 0 ), stress( 1, 1 ), stress( 2, 2 ) };
    Matrix3d dS_dF;
    dS_dF << dS_dF_( 0, 0, 0, 0 ), dS_dF_( 0, 0, 1, 1 ), dS_dF_( 0, 0, 2, 2 ), dS_dF_( 1, 1, 0, 0 ),
      dS_dF_( 1, 1, 1, 1 ), dS_dF_( 1, 1, 2, 2 ), dS_dF_( 2, 2, 0, 0 ), dS_dF_( 2, 2, 1, 1 ), dS_dF_( 2, 2, 2, 2 );

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
