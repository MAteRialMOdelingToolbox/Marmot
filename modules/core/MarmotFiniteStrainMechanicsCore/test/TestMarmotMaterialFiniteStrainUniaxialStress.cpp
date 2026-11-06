#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include "Marmot/MarmotMaterialFiniteStrainFactory.h"
#include "Marmot/MarmotTesting.h"
#include <cmath>
#include <functional>
#include <memory>
#include <string>
#include <vector>

using namespace Marmot;
using namespace Marmot::Testing;
using namespace Marmot::Math::Tensors::FastorTensors;

namespace {

  using Material = Marmot::Materials::MaterialFiniteStrain;

  /// the materials keep a pointer to their properties, which must outlive them
  std::vector< std::vector< double > > propertyStorage;

  std::unique_ptr< Material > createMaterial( const std::string& name, const std::vector< double >& properties_ )
  {
    const auto& properties = propertyStorage.emplace_back( properties_ );
    return std::unique_ptr< Material >(
      Marmot::Registry::MaterialFiniteStrainFactory::createMaterial( name,
                                                                     properties.data(),
                                                                     static_cast< int >( properties.size() ),
                                                                     1 ) );
  }

  /// One uniaxial stress update at the axial stretch lambda, starting from the given state
  std::tuple< double, double, Material::Deformation< 3 > > uniaxialUpdate( const Material&        material,
                                                                           double                 lambda,
                                                                           std::vector< double >& stateVars,
                                                                           const Tensor33d&       tauOld,
                                                                           double                 lateralGuess )
  {
    Material::ConstitutiveResponse< 3 > response( tauOld, 0.0, 0.0, stateVars.data() );
    Material::Deformation< 3 >          deformation{ Tensor33d( 0.0 ) };
    deformation.F( 0, 0 ) = lambda;
    deformation.F( 1, 1 ) = lateralGuess;
    deformation.F( 2, 2 ) = lateralGuess;
    double dTau11_dF11    = 0.0;
    material.computeUniaxialStress( response, dTau11_dF11, deformation, { 0.0, 1.0 } );
    throwExceptionOnFailure( checkIfEqual( response.tau( 1, 1 ), 0.0, 1e-8 ) &&
                               checkIfEqual( response.tau( 2, 2 ), 0.0, 1e-8 ),
                             Marmot::Utilities::MakeString()
                               << "lateral stresses must vanish under uniaxial stress, got " << response.tau( 1, 1 )
                               << ", " << response.tau( 2, 2 ) );
    return { response.tau( 0, 0 ), dTau11_dF11, deformation };
  }

  /// Central difference of the condensed axial stress w.r.t. the axial stretch, from the same old state
  double numericalTangent( const Material&              material,
                           double                       lambda,
                           const std::vector< double >& stateVarsOld,
                           const Tensor33d&             tauOld,
                           double                       h = 1e-7 )
  {
    auto svPlus              = stateVarsOld;
    auto svMinus             = stateVarsOld;
    auto [tauPlus, _p, _dp]  = uniaxialUpdate( material, lambda + h, svPlus, tauOld, 1.0 );
    auto [tauMinus, _m, _dm] = uniaxialUpdate( material, lambda - h, svMinus, tauOld, 1.0 );
    return ( tauPlus - tauMinus ) / ( 2 * h );
  }

} // namespace

/// Compressible Neo-Hooke: vanishing lateral stress, transversely isotropic lateral stretch, and a condensed
/// tangent that matches the central difference of the converged axial stress
void testNeoHookeUniaxialStressAndTangent()
{
  const double          K = 3000., G = 1000.;
  auto                  material = createMaterial( "COMPRESSIBLENEOHOOKE", { K, G, 1.0 } );
  std::vector< double > stateVars( std::max( 1, material->getNumberOfRequiredStateVars() ), 0.0 );

  for ( const double lambda : { 0.7, 0.95, 1.0, 1.05, 1.6 } ) {
    auto [tau11, dTau11_dF11, deformation] = uniaxialUpdate( *material, lambda, stateVars, Tensor33d( 0.0 ), 1.0 );

    throwExceptionOnFailure( checkIfEqual( deformation.F( 1, 1 ), deformation.F( 2, 2 ), 1e-10 ),
                             "lateral stretches of an isotropic material must be equal" );

    const double numTangent = numericalTangent( *material, lambda, stateVars, Tensor33d( 0.0 ) );
    throwExceptionOnFailure( checkIfEqual( dTau11_dF11, numTangent, 1e-5 * std::abs( numTangent ) ),
                             Marmot::Utilities::MakeString() << "condensed tangent " << dTau11_dF11 << " != numerical "
                                                             << numTangent << " at lambda=" << lambda );
    if ( lambda == 1.0 ) {
      // the small strain limit: Young's modulus
      const double E = 9. * K * G / ( 3. * K + G );
      throwExceptionOnFailure( checkIfEqual( tau11, 0.0, 1e-10 ) && checkIfEqual( dTau11_dF11, E, 1e-6 * E ),
                               Marmot::Utilities::MakeString()
                                 << "undeformed uniaxial tangent " << dTau11_dF11 << " != E=" << E );
    }
    else
      throwExceptionOnFailure( ( tau11 > 0 ) == ( lambda > 1 ), "axial stress must have the sign of the strain" );
  }
}

/// J2 plasticity without hardening: after yielding, the uniaxial Kirchhoff stress equals the yield stress, the state
/// is only updated once, and the condensed elastoplastic tangent matches the numerical one
void testJ2PlasticityUniaxialStress()
{
  const double K = 175000., G = 80769., fy = 300.;
  //                                                           K  G  fy  fyInf  eta  H  impl  rho
  auto material = createMaterial( "FINITESTRAINJ2PLASTICITY", { K, G, fy, fy, 0.0, 1000., 1, 7.8e-9 } );

  const int             nStateVars = material->getNumberOfRequiredStateVars();
  std::vector< double > stateVars( nStateVars, 0.0 );
  material->initializeYourself( stateVars.data(), nStateVars );

  // elastic: E (lambda-1) to leading order
  const double E       = 9. * K * G / ( 3. * K + G );
  const double lamEl   = 1.0 + 0.5 * fy / E;
  auto         svEl    = stateVars;
  auto [tauEl, CEl, _] = uniaxialUpdate( *material, lamEl, svEl, Tensor33d( 0.0 ), 1.0 );
  throwExceptionOnFailure( checkIfEqual( tauEl, 0.5 * fy, 1e-3 * fy ) && checkIfEqual( CEl, E, 1e-2 * E ),
                           Marmot::Utilities::MakeString() << "elastic uniaxial response " << tauEl << ", " << CEl );

  // plastic, in one large step from the virgin state
  const double lamPl               = 1.01;
  auto         svPl                = stateVars;
  auto [tauPl, CPl, deformationPl] = uniaxialUpdate( *material, lamPl, svPl, Tensor33d( 0.0 ), 1.0 );
  const double alphaP              = *material->getStateView( "alphaP", svPl.data() ).stateLocation;
  throwExceptionOnFailure( alphaP > 0, "the plastic step must yield" );
  const double expectedYieldStress = fy + 1000. * alphaP;
  throwExceptionOnFailure( checkIfEqual( tauPl, expectedYieldStress, 1e-6 * fy ),
                           Marmot::Utilities::MakeString() << "uniaxial Kirchhoff stress " << tauPl
                                                           << " != current yield stress " << expectedYieldStress );

  // plastic flow is isochoric: lateral contraction larger than the elastic one
  throwExceptionOnFailure( deformationPl.F( 1, 1 ) < 1.0 - 0.3 * ( lamPl - 1.0 ), "plastic lateral contraction" );

  const double numTangent = numericalTangent( *material, lamPl, stateVars, Tensor33d( 0.0 ) );
  throwExceptionOnFailure( checkIfEqual( CPl, numTangent, 1e-4 * std::abs( numTangent ) + 1e-3 ),
                           Marmot::Utilities::MakeString()
                             << "elastoplastic condensed tangent " << CPl << " != numerical " << numTangent );
}

/// a wrong initial guess must still converge, and a non-positive lateral stretch must never be tried
void testRobustInitialGuess()
{
  auto                  material = createMaterial( "COMPRESSIBLENEOHOOKE", { 3000., 1000., 1.0 } );
  std::vector< double > stateVars( std::max( 1, material->getNumberOfRequiredStateVars() ), 0.0 );
  auto [tauGood, CGood, defGood] = uniaxialUpdate( *material, 1.3, stateVars, Tensor33d( 0.0 ), 1.0 );
  auto [tauBad, CBad, defBad]    = uniaxialUpdate( *material, 1.3, stateVars, Tensor33d( 0.0 ), 3.0 );
  throwExceptionOnFailure( checkIfEqual( tauGood, tauBad, 1e-8 ) && checkIfEqual( CGood, CBad, 1e-6 ),
                           "the converged uniaxial state must not depend on the initial guess" );
}

int main()
{
  const std::vector< std::function< void() > > tests = {
    testNeoHookeUniaxialStressAndTangent,
    testJ2PlasticityUniaxialStress,
    testRobustInitialGuess,
  };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
