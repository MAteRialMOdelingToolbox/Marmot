#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"
#include "Marmot/MarmotTesting.h"
#include <cmath>
#include <memory>
#include <string>
#include <vector>

using namespace Marmot;
using namespace Marmot::Testing;
using namespace Marmot::Math;
using Marmot::Utilities::MakeString;
using MHE = Marmot::Materials::MaterialHypoElastic;

namespace {

  const double E = 210000., nu = 0.3, G = E / ( 2 * ( 1 + nu ) );

  std::unique_ptr< MHE > makeMaterial( const std::string& name, const std::vector< double >& properties )
  {
    return std::unique_ptr< MHE >(
      Marmot::Registry::MaterialHypoElasticFactory::createMaterial( name,
                                                                    properties.data(),
                                                                    int( properties.size() ),
                                                                    1 ) );
  }

  /// a material point with its own state, for beam stress updates from a committed state
  struct Point {
    std::vector< double >  properties; ///< the material keeps a pointer to its properties
    std::unique_ptr< MHE > material;
    std::vector< double >  stateVars;
    MHE::stateBeam         state;

    Point( const std::string& name, const std::vector< double >& properties_ )
      : properties( properties_ ),
        material( makeMaterial( name, properties ) ),
        stateVars( material->getNumberOfRequiredStateVars(), 0.0 )
    {
      material->initializeYourself( stateVars.data(), int( stateVars.size() ) );
      state.stateVars = stateVars.data();
    }

    /// stress and tangent for the strain increment, from the committed state; commits if requested
    std::pair< Vector3d, Matrix3d > update( const Vector3d& dStrain, bool commit )
    {
      auto           backupVars  = stateVars;
      auto           backupState = state;
      Matrix3d       C;
      MHE::stateBeam trial = state;
      trial.stateVars      = stateVars.data();
      material->computeBeamStress( trial, C, dStrain, { 0.0, 1.0 } );
      if ( commit )
        state = trial, state.stateVars = stateVars.data();
      else
        stateVars = backupVars, state = backupState, state.stateVars = stateVars.data();
      return { trial.stress, C };
    }
  };

} // namespace

/// linear elastic: sigma11 = E eps11, sigma1k = G gamma1k, tangent diag(E, G, G)
void testLinearElasticClosedForm()
{
  Point          p( "LINEARELASTIC", { E, nu } );
  const Vector3d dStrain( 1e-3, 2e-3, -0.5e-3 );
  const auto [stress, C] = p.update( dStrain, true );
  throwExceptionOnFailure( checkIfEqual< double >( stress, Vector3d( E * 1e-3, G * 2e-3, -G * 0.5e-3 ), 1e-9 ),
                           MakeString() << "linear elastic beam stress: " << stress.transpose() );
  throwExceptionOnFailure( checkIfEqual< double >( C, Vector3d( E, G, G ).asDiagonal().toDenseMatrix(), 1e-6 ),
                           MakeString() << "linear elastic beam tangent:\n"
                                        << C );
}

/// without shear strains, the beam stress equals the uniaxial stress (von Mises, plastic)
void testEqualsUniaxialWithoutShear()
{
  const std::vector< double > vm = { E, nu, 500., 1000., 0., 1., 7.85e-9 };
  Point                       beam( "VONMISES", vm );
  auto                        uni = makeMaterial( "VONMISES", vm );
  std::vector< double >       uniVars( uni->getNumberOfRequiredStateVars(), 0.0 );
  MHE::state1D                s1;
  s1.stateVars = uniVars.data();
  for ( int i = 0; i < 5; i++ ) {
    double C1 = 0;
    uni->computeUniaxialStress( s1, C1, 1e-3, { 0.0, 1.0 } );
    const auto [stress, C] = beam.update( Vector3d( 1e-3, 0, 0 ), true );
    throwExceptionOnFailure( std::abs( stress( 0 ) - s1.stress ) < 1e-8 * std::abs( s1.stress ) &&
                               std::abs( C( 0, 0 ) - C1 ) < 1e-6 * std::abs( C1 ),
                             MakeString() << "beam stress != uniaxial stress: " << stress( 0 ) << " " << s1.stress );
  }
}

/// von Mises in combined tension and shear: the condensed tangent matches central differences in a plastic state, and
/// perfect plasticity reaches the yield surface sigma11^2 + 3 (sigma12^2 + sigma13^2) = fy^2
void testVonMisesTangentAndYieldSurface()
{
  Point p( "VONMISES", { E, nu, 500., 2000., 0., 1., 7.85e-9 } );
  p.update( Vector3d( 2e-3, 1e-3, 0.5e-3 ), true );
  const Vector3d dStrain( 1e-3, 1.5e-3, -0.2e-3 );
  const auto [stress, C] = p.update( dStrain, false );
  Matrix3d     numC;
  const double h = 1e-9;
  for ( int j = 0; j < 3; j++ ) {
    Vector3d dp = dStrain, dm = dStrain;
    dp( j ) += h;
    dm( j ) -= h;
    numC.col( j ) = ( p.update( dp, false ).first - p.update( dm, false ).first ) / ( 2 * h );
  }
  throwExceptionOnFailure( checkIfEqual< double >( C, numC, 1e-5 * numC.cwiseAbs().maxCoeff() ),
                           MakeString() << "von Mises beam tangent:\n"
                                        << C << "\nnumerical:\n"
                                        << numC );
  throwExceptionOnFailure( p.state.dissipation >= 0, "dissipation" );

  Point    perfect( "VONMISES", { E, nu, 500., 0., 0., 1., 7.85e-9 } );
  Vector3d s;
  for ( int i = 0; i < 20; i++ )
    s = perfect.update( Vector3d( 1e-3, 1e-3, 0.5e-3 ), true ).first;
  const double equivalent = std::sqrt( s( 0 ) * s( 0 ) + 3 * ( s( 1 ) * s( 1 ) + s( 2 ) * s( 2 ) ) );
  throwExceptionOnFailure( std::abs( equivalent - 500. ) < 1e-6 * 500.,
                           MakeString() << "von Mises yield surface in beam stress: " << equivalent );
}

int main()
{
  executeTestsAndCollectExceptions( {
    testLinearElasticClosedForm,
    testEqualsUniaxialWithoutShear,
    testVonMisesTangentAndYieldSurface,
  } );
  return 0;
}
