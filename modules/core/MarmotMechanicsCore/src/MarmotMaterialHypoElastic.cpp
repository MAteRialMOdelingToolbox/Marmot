#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotVoigt.h"
#include <algorithm>
#include <array>

using namespace Eigen;

void Marmot::Materials::MaterialHypoElastic::setCharacteristicElementLength( double length )
{
  characteristicElementLength = length;
}

template < int nP, int nF >
void Marmot::Materials::MaterialHypoElastic::computeStressWithFreeComponents(
  state3D&                         state,
  Eigen::Matrix< double, nP, nP >& dStressP_dStrainP,
  Marmot::Math::Vector6d&          dStrain,
  const std::array< int, nP >&     p,
  const std::array< int, nF >&     f,
  const timeInfo&                  timeInfo ) const
{
  using namespace Marmot::Math;
  using VectorP = Eigen::Matrix< double, nP, 1 >;
  using VectorF = Eigen::Matrix< double, nF, 1 >;

  // the start of the increment: every iteration evaluates the material from it
  Map< VectorXd > stateVars( state.stateVars, stateLayout.totalSize() );
  const VectorXd  stateVarsOld            = stateVars;
  const Vector6d  stressOld               = state.stress;
  const double    elasticEnergyDensityOld = state.elasticEnergyDensity;
  const double    dissipationOld          = state.dissipation;
  VectorP         stressPOld;
  for ( int i = 0; i < nP; i++ )
    stressPOld( i ) = stressOld( p[i] );
  const double stressScaleOld = std::max( 1.0, stressPOld.cwiseAbs().maxCoeff() );

  constexpr double tolerance = 1e-12, relaxedTolerance = 1e-9; // relative to the stress scale
  constexpr int    nRelaxed = 8, nMaxIterations = 15;

  Matrix6d dStress_dStrain;
  for ( int iteration = 0;; iteration++ ) {
    stateVars = stateVarsOld;
    state.stress.setZero();
    for ( int i = 0; i < nP; i++ )
      state.stress( p[i] ) = stressOld( p[i] );
    state.elasticEnergyDensity = elasticEnergyDensityOld;
    state.dissipation          = dissipationOld;

    dStress_dStrain.setZero();
    computeStress( state, dStress_dStrain, dStrain, timeInfo );

    VectorF                         stressF;
    VectorP                         stressP;
    Eigen::Matrix< double, nF, nF > Cff;
    Eigen::Matrix< double, nF, nP > Cfp;
    Eigen::Matrix< double, nP, nF > Cpf;
    Eigen::Matrix< double, nP, nP > Cpp;
    for ( int i = 0; i < nF; i++ ) {
      stressF( i ) = state.stress( f[i] );
      for ( int j = 0; j < nF; j++ )
        Cff( i, j ) = dStress_dStrain( f[i], f[j] );
      for ( int j = 0; j < nP; j++ )
        Cfp( i, j ) = dStress_dStrain( f[i], p[j] );
    }
    for ( int i = 0; i < nP; i++ ) {
      stressP( i ) = state.stress( p[i] );
      for ( int j = 0; j < nF; j++ )
        Cpf( i, j ) = dStress_dStrain( p[i], f[j] );
      for ( int j = 0; j < nP; j++ )
        Cpp( i, j ) = dStress_dStrain( p[i], p[j] );
    }

    const auto   CffLU       = Cff.fullPivLu();
    const double stressScale = std::max( stressScaleOld, stressP.cwiseAbs().maxCoeff() );
    const double residual    = stressF.cwiseAbs().maxCoeff();
    if ( residual <= tolerance * stressScale ||
         ( iteration >= nRelaxed && residual <= relaxedTolerance * stressScale ) ) {
      // static condensation of the free components
      dStressP_dStrainP = Cpp - Cpf * CffLU.solve( Cfp );
      return;
    }
    if ( iteration >= nMaxIterations || !CffLU.isInvertible() ) {
      stateVars = stateVarsOld;
      throw Marmot::Utilities::StressUpdateFailed( Marmot::Utilities::MakeString()
                                                   << __PRETTY_FUNCTION__ << ": the reduced stress iteration (" << nF
                                                   << " free components) did not converge, residual " << residual );
    }
    const VectorF dStrainF = CffLU.solve( stressF );
    for ( int i = 0; i < nF; i++ )
      dStrain( f[i] ) -= dStrainF( i );
  }
}

// the reduced stress states of the interface: plane stress, uniaxial stress, beam stress
template void Marmot::Materials::MaterialHypoElastic::computeStressWithFreeComponents< 3, 1 >(
  state3D&,
  Eigen::Matrix< double, 3, 3 >&,
  Marmot::Math::Vector6d&,
  const std::array< int, 3 >&,
  const std::array< int, 1 >&,
  const timeInfo& ) const;
template void Marmot::Materials::MaterialHypoElastic::computeStressWithFreeComponents< 1, 2 >(
  state3D&,
  Eigen::Matrix< double, 1, 1 >&,
  Marmot::Math::Vector6d&,
  const std::array< int, 1 >&,
  const std::array< int, 2 >&,
  const timeInfo& ) const;
template void Marmot::Materials::MaterialHypoElastic::computeStressWithFreeComponents< 3, 3 >(
  state3D&,
  Eigen::Matrix< double, 3, 3 >&,
  Marmot::Math::Vector6d&,
  const std::array< int, 3 >&,
  const std::array< int, 3 >&,
  const timeInfo& ) const;

void Marmot::Materials::MaterialHypoElastic::computePlaneStress( state2D&                      state2D,
                                                                 Marmot::Math::Matrix3d&       dStress_dStrain2D,
                                                                 const Marmot::Math::Vector3d& dStrain2D,
                                                                 const timeInfo&               timeInfo ) const
{
  using namespace Marmot::ContinuumMechanics::Voigt;
  // prescribed: eps11, eps22, gamma12; free: sigma33 = 0 (the out-of-plane shear strains stay zero)
  static constexpr std::array< int, 3 > p = { 0, 1, 3 };
  static constexpr std::array< int, 1 > f = { 2 };

  Marmot::Math::Vector6d dStrain = make3DVoigt< VoigtSize::TwoD >( dStrain2D );
  dStrain( 2 )                   = -dStrain2D( 0 ) - dStrain2D( 1 ); // initial guess: isochoric
  state3D state( make3DVoigt< VoigtSize::TwoD >( state2D.stress ),
                 state2D.elasticEnergyDensity,
                 state2D.dissipation,
                 state2D.stateVars );

  computeStressWithFreeComponents< 3, 1 >( state, dStress_dStrain2D, dStrain, p, f, timeInfo );

  state2D.stress               = reduce3DVoigt< VoigtSize::TwoD >( state.stress );
  state2D.elasticEnergyDensity = state.elasticEnergyDensity;
  state2D.dissipation          = state.dissipation;
}

void Marmot::Materials::MaterialHypoElastic::computeUniaxialStress( state1D&        state1D,
                                                                    double&         dStress_dStrain1D,
                                                                    const double    dStrain1D,
                                                                    const timeInfo& timeInfo ) const
{
  // prescribed: eps11; free: sigma22 = sigma33 = 0 (the shear strains stay zero)
  static constexpr std::array< int, 1 > p = { 0 };
  static constexpr std::array< int, 2 > f = { 1, 2 };

  Marmot::Math::Vector6d dStrain = Marmot::Math::Vector6d::Zero();
  dStrain( 0 )                   = dStrain1D;
  Marmot::Math::Vector6d stress  = Marmot::Math::Vector6d::Zero();
  stress( 0 )                    = state1D.stress;
  state3D state( stress, state1D.elasticEnergyDensity, state1D.dissipation, state1D.stateVars );

  Eigen::Matrix< double, 1, 1 > tangent;
  computeStressWithFreeComponents< 1, 2 >( state, tangent, dStrain, p, f, timeInfo );

  state1D.stress               = state.stress( 0 );
  state1D.elasticEnergyDensity = state.elasticEnergyDensity;
  state1D.dissipation          = state.dissipation;
  dStress_dStrain1D            = tangent( 0, 0 );
}

void Marmot::Materials::MaterialHypoElastic::computeBeamStress( stateBeam&                    stateBeam,
                                                                Marmot::Math::Matrix3d&       dStress_dStrainBeam,
                                                                const Marmot::Math::Vector3d& dStrainBeam,
                                                                const timeInfo&               timeInfo ) const
{
  // prescribed: eps11, gamma12, gamma13; free: sigma22 = sigma33 = sigma23 = 0
  static constexpr std::array< int, 3 > p = { 0, 3, 4 };
  static constexpr std::array< int, 3 > f = { 1, 2, 5 };

  Marmot::Math::Vector6d dStrain = Marmot::Math::Vector6d::Zero();
  Marmot::Math::Vector6d stress  = Marmot::Math::Vector6d::Zero();
  for ( int i = 0; i < 3; i++ ) {
    dStrain( p[i] ) = dStrainBeam( i );
    stress( p[i] )  = stateBeam.stress( i );
  }
  state3D state( stress, stateBeam.elasticEnergyDensity, stateBeam.dissipation, stateBeam.stateVars );

  computeStressWithFreeComponents< 3, 3 >( state, dStress_dStrainBeam, dStrain, p, f, timeInfo );

  for ( int i = 0; i < 3; i++ )
    stateBeam.stress( i ) = state.stress( p[i] );
  stateBeam.elasticEnergyDensity = state.elasticEnergyDensity;
  stateBeam.dissipation          = state.dissipation;
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
