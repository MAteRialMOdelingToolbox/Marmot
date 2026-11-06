#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotLowerDimensionalStress.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotVoigt.h"
#include <algorithm>
#include <array>

using namespace Eigen;

void Marmot::Materials::MaterialHypoElastic::setCharacteristicElementLength( double length )
{
  characteristicElementLength = length;
}

void Marmot::Materials::MaterialHypoElastic::computePlaneStress( state2D&                      state2D_,
                                                                 Marmot::Math::Matrix3d&       dStress_dStrain2D_,
                                                                 const Marmot::Math::Vector3d& dStrain2D_,
                                                                 const timeInfo&               timeInfo ) const
{
  using namespace Marmot;
  using namespace ContinuumMechanics::Voigt;

  Map< const Matrix< double, 3, 1 > > dStrain2D( dStrain2D_.data() );
  Map< Matrix< double, 3, 1 > >       stress2D( state2D_.stress.data() );
  Map< Matrix< double, 3, 3 > >       dStress_dStrain2D( dStress_dStrain2D_.data() );
  Map< VectorXd >                     stateVars( state2D_.stateVars, this->getNumberOfRequiredStateVars() );

  Marmot::Math::Matrix6d dStress_dStrain3D;

  VectorXd               stateVarsOld  = stateVars;
  Marmot::Math::Vector6d dStrain3DTemp = Marmot::ContinuumMechanics::Voigt::make3DVoigt< VoigtSize::TwoD >( dStrain2D );

  // assumption of isochoric deformation for initial guess
  dStrain3DTemp( 2 ) = ( -dStrain2D( 0 ) - dStrain2D( 1 ) );

  state3D state( Marmot::ContinuumMechanics::Voigt::make3DVoigt< VoigtSize::TwoD >( stress2D ),
                 state2D_.elasticEnergyDensity,
                 state2D_.dissipation,
                 stateVars.data() );

  int planeStressCount = 1;
  while ( true ) {

    stateVars = stateVarsOld;

    state.stress               = Marmot::ContinuumMechanics::Voigt::make3DVoigt< VoigtSize::TwoD >( stress2D );
    state.elasticEnergyDensity = state2D_.elasticEnergyDensity;
    state.dissipation          = state2D_.dissipation;

    computeStress( state, dStress_dStrain3D, dStrain3DTemp, timeInfo );

    double residual = state.stress.array().abs()[2];

    if ( residual < 1.e-10 || ( planeStressCount > 7 && residual < 1e-8 ) ) {
      break;
    }

    double tangentCompliance = 1. / dStress_dStrain3D( 2, 2 );
    if ( Math::isNaN( tangentCompliance ) || std::abs( tangentCompliance ) > 1e10 )
      tangentCompliance = 1e10;

    dStrain3DTemp( 2 ) -= tangentCompliance * state.stress( 2 );

    planeStressCount += 1;
    if ( planeStressCount > 13 ) {
      Marmot::Utilities::Journal::warningToMSG( "PlaneStressWrapper requires cutback" );
      throw Marmot::Utilities::StressUpdateFailed( "Plane stress iteration did not converge" );
    }
  }

  state2D_.stress               = ContinuumMechanics::Voigt::reduce3DVoigt< VoigtSize::TwoD >( state.stress );
  state2D_.elasticEnergyDensity = state.elasticEnergyDensity;
  state2D_.dissipation          = state.dissipation;
  state2D_.stateVars            = stateVars.data();

  dStress_dStrain2D_ = ContinuumMechanics::PlaneStress::getPlaneStressTangent( dStress_dStrain3D );
}

void Marmot::Materials::MaterialHypoElastic::computeUniaxialStress( state1D&        state1D_,
                                                                    double&         dStress_dStrain1D_,
                                                                    const double    dStrain1D_,
                                                                    const timeInfo& timeInfo ) const
{
  using namespace Marmot;
  using namespace ContinuumMechanics::Voigt;

  Map< const Matrix< double, 1, 1 > > dStrain1D( &dStrain1D_ );
  Map< Matrix< double, 1, 1 > >       stress1D( &state1D_.stress );
  Map< VectorXd >                     stateVars( state1D_.stateVars, stateLayout.totalSize() );

  Marmot::Math::Matrix6d dStress_dStrain3D;

  VectorXd               stateVarsOld  = stateVars;
  Marmot::Math::Vector6d dStrain3DTemp = Marmot::ContinuumMechanics::Voigt::make3DVoigt< VoigtSize::OneD >( dStrain1D );

  state3D state( Marmot::ContinuumMechanics::Voigt::make3DVoigt< VoigtSize::OneD >( stress1D ),
                 state1D_.elasticEnergyDensity,
                 state1D_.dissipation,
                 stateVars.data() );

  int count = 1;
  while ( true ) {
    stateVars = stateVarsOld;

    state.stress               = Marmot::ContinuumMechanics::Voigt::make3DVoigt< VoigtSize::OneD >( stress1D );
    state.elasticEnergyDensity = state1D_.elasticEnergyDensity;
    state.dissipation          = state1D_.dissipation;

    computeStress( state, dStress_dStrain3D, dStrain3DTemp, timeInfo );

    const double residual = state.stress.array().abs().segment( 1, 2 ).sum();

    if ( residual < 1.e-13 || ( count > 7 && residual < 1e-10 ) ) {
      break;
    }

    dStrain3DTemp.segment< 2 >( 1 ) -= dStress_dStrain3D.block< 2, 2 >( 1, 1 ).colPivHouseholderQr().solve(
      state.stress.segment< 2 >( 1 ) );

    count += 1;
    if ( count > 13 ) {
      Marmot::Utilities::Journal::warningToMSG( "UniaxialStressWrapper requires cutback" );
      throw Marmot::Utilities::StressUpdateFailed( "uniaxial stress iteration did not converge" );
    }
  }

  state1D_.stress               = ContinuumMechanics::Voigt::reduce3DVoigt< VoigtSize::OneD >( state.stress )[0];
  state1D_.elasticEnergyDensity = state.elasticEnergyDensity;
  state1D_.dissipation          = state.dissipation;
  state1D_.stateVars            = stateVars.data();

  dStress_dStrain1D_ = ContinuumMechanics::UniaxialStress::getUniaxialStressTangent( dStress_dStrain3D );
}
void Marmot::Materials::MaterialHypoElastic::computeBeamStress( stateBeam&                    stateBeam_,
                                                                Marmot::Math::Matrix3d&       dStress_dStrainBeam,
                                                                const Marmot::Math::Vector3d& dStrainBeam,
                                                                const timeInfo&               timeInfo ) const
{
  using namespace Marmot::Math;
  // prescribed: eps11, gamma12, gamma13; free: sigma22 = sigma33 = sigma23 = 0
  static constexpr std::array< int, 3 > p = { 0, 3, 4 };
  static constexpr std::array< int, 3 > f = { 1, 2, 5 };

  Map< VectorXd > stateVars( stateBeam_.stateVars, stateLayout.totalSize() );
  const VectorXd  stateVarsOld = stateVars;

  Vector6d dStrain = Vector6d::Zero();
  for ( int i = 0; i < 3; i++ )
    dStrain( p[i] ) = dStrainBeam( i );

  const double stressScale = std::max( 1.0, stateBeam_.stress.cwiseAbs().maxCoeff() );
  Matrix6d     C;
  state3D      state;
  for ( int iteration = 0;; iteration++ ) {
    stateVars    = stateVarsOld;
    state.stress = Vector6d::Zero();
    for ( int i = 0; i < 3; i++ )
      state.stress( p[i] ) = stateBeam_.stress( i );
    state.elasticEnergyDensity = stateBeam_.elasticEnergyDensity;
    state.dissipation          = stateBeam_.dissipation;
    state.stateVars            = stateVars.data();

    C.setZero();
    computeStress( state, C, dStrain, timeInfo );

    Vector3d r, sP;
    Matrix3d Cff, Cfp, Cpf, Cpp;
    for ( int i = 0; i < 3; i++ ) {
      r( i )  = state.stress( f[i] );
      sP( i ) = state.stress( p[i] );
      for ( int j = 0; j < 3; j++ ) {
        Cff( i, j ) = C( f[i], f[j] );
        Cfp( i, j ) = C( f[i], p[j] );
        Cpf( i, j ) = C( p[i], f[j] );
        Cpp( i, j ) = C( p[i], p[j] );
      }
    }
    const auto   CffLU    = Cff.fullPivLu();
    const double scale    = std::max( stressScale, sP.cwiseAbs().maxCoeff() );
    const double residual = r.cwiseAbs().maxCoeff();
    if ( residual <= 1e-12 * scale || ( iteration > 7 && residual <= 1e-9 * scale ) ) {
      stateBeam_.stress               = sP;
      stateBeam_.elasticEnergyDensity = state.elasticEnergyDensity;
      stateBeam_.dissipation          = state.dissipation;
      dStress_dStrainBeam             = Cpp - Cpf * CffLU.solve( Cfp );
      return;
    }
    if ( iteration >= 15 || !CffLU.isInvertible() ) {
      Marmot::Utilities::Journal::warningToMSG( "BeamStressWrapper requires cutback" );
      throw Marmot::Utilities::StressUpdateFailed( "beam stress iteration did not converge" );
    }
    const Vector3d ddStrain = CffLU.solve( r );
    for ( int i = 0; i < 3; i++ )
      dStrain( f[i] ) -= ddStrain( i );
  }
}

double Marmot::Materials::MaterialHypoElastic::getMaximumWaveSpeed( const state3D& state ) const
{
  const int nStateVars = getNumberOfRequiredStateVars();

  std::vector< double > stateVarsCopy( nStateVars, 0.0 );
  if ( state.stateVars != nullptr && nStateVars > 0 ) {
    std::copy_n( state.stateVars, nStateVars, stateVarsCopy.begin() );
  }

  state3D stateCopy( state.stress, state.elasticEnergyDensity, state.dissipation, stateVarsCopy.data() );

  Marmot::Math::Matrix6d dStress_dStrain = Marmot::Math::Matrix6d::Zero();
  const auto             dStrain         = Marmot::Math::Vector6d::Zero();
  const timeInfo         timeInfo{ 0.0, 1.0 };

  computeStress( stateCopy, dStress_dStrain, dStrain, timeInfo );

  const double maxStiffnessDiagonal = std::max( { dStress_dStrain( 0, 0 ),
                                                  dStress_dStrain( 1, 1 ),
                                                  dStress_dStrain( 2, 2 ),
                                                  dStress_dStrain( 3, 3 ),
                                                  dStress_dStrain( 4, 4 ),
                                                  dStress_dStrain( 5, 5 ) } );
  const double density              = getDensity( stateCopy.stateVars );

  return density > 0.0 ? std::sqrt( std::max( 0.0, maxStiffnessDiagonal ) / density ) : 0.0;
}
