#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"
#include "Marmot/MarmotMaterialPointSolverHypoElastic.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotTypedefs.h"
#include <cmath>
#include <memory>
#include <utility>

using namespace Marmot::Testing;
using namespace Marmot::Materials;

void testVonMisesCoordinateInvariance()
{
  // material properties
  std::vector< double > materialProperties = { 210000., 0.3, 200., 2100., 20., 20 };
  auto                  solveropts         = Marmot::Materials::MaterialPointSolverHypoElastic::SolverOptions();
  std::string           matName            = "VONMISES";

  // create material point solver instance
  auto solver = Marmot::Materials::MaterialPointSolverHypoElastic( matName,
                                                                   &materialProperties[0],
                                                                   static_cast< int >( materialProperties.size() ),
                                                                   solveropts );

  // define a step
  Marmot::Materials::MaterialPointSolverHypoElastic::Step step;

  // define step parameters
  step.isStrainComponentControlled = { true, true, true, true, true, true };
  step.isStressComponentControlled = { false, false, false, false, false, false };
  step.stressIncrementTarget       = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  step.strainIncrementTarget       = { 0.01, 0., 0., 0., 0.03, 0.0 };
  step.dTStart                     = 1.0;

  // add step to solver
  solver.addStep( step );
  // check coordinate invariance with Turbokreisel
  throwExceptionOnFailure( spinTurbokreisel( solver, 1e-10, 1e-8 ), "Turbokreisel failed!" );
}

void testVonMises()
{
  // material properties
  std::vector< double > materialProperties = { 210000., 0.3, 200., 2100., 20., 20 };
  auto                  solveropts         = Marmot::Materials::MaterialPointSolverHypoElastic::SolverOptions();
  std::string           matName            = "VONMISES";

  // create material point solver instance
  auto solver = Marmot::Materials::MaterialPointSolverHypoElastic( matName,
                                                                   &materialProperties[0],
                                                                   static_cast< int >( materialProperties.size() ),
                                                                   solveropts );

  // define a step
  Marmot::Materials::MaterialPointSolverHypoElastic::Step step;

  // define step parameters
  step.isStrainComponentControlled = { true, true, true, true, true, true };
  step.isStressComponentControlled = { false, false, false, false, false, false };
  step.stressIncrementTarget       = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  step.strainIncrementTarget       = { 0.00839244, 0.00089344, -0.00703916, 0.00013635, 0.00160548, 0.00572825 };
  step.dTStart                     = 1.0;

  // add step to solver
  solver.addStep( step );

  solver.solve();

  // read history
  auto history = solver.getHistory();
  // reference solution
  const Marmot::Math::Vector6d stressTarget =
    { 511.262747695, 395.408929527, 272.856322779, 1.0532516407, 12.4017194288, 44.2485420671 };

  // compare solutions
  throwExceptionOnFailure( checkIfEqual< double >( history.back().stress, stressTarget, 1e-9 ),
                           "comparison with reference solution failed" );
}

/// a single, very large strain increment of a perfectly plastic material (trial stress ~ 1e5 f_y): the return mapping
/// converges (relative tolerance; an absolute one is below the round-off of such trial stresses) and the stress lies on
/// the yield surface
void testLargeIncrementPerfectPlasticity()
{
  using MHE = Marmot::Materials::MaterialHypoElastic;
  using Marmot::Utilities::MakeString;
  const double                fy         = 235.;
  const std::vector< double > properties = { 210000., 0.3, fy, 0., 0., 1., 7.85e-9 };
  auto                        makePoint  = [&]() {
    auto material = std::unique_ptr< MHE >(
      Marmot::Registry::MaterialHypoElasticFactory::createMaterial( "VONMISES",
                                                                    properties.data(),
                                                                    static_cast< int >( properties.size() ),
                                                                    1 ) );
    std::vector< double > stateVars( material->getNumberOfRequiredStateVars(), 0.0 );
    material->initializeYourself( stateVars.data(), static_cast< int >( stateVars.size() ) );
    return std::make_pair( std::move( material ), std::move( stateVars ) );
  };
  auto vonMisesStress = []( const Marmot::Math::Vector6d& s ) {
    const double p = ( s( 0 ) + s( 1 ) + s( 2 ) ) / 3;
    return std::sqrt(
      1.5 * ( ( s( 0 ) - p ) * ( s( 0 ) - p ) + ( s( 1 ) - p ) * ( s( 1 ) - p ) + ( s( 2 ) - p ) * ( s( 2 ) - p ) ) +
      3 * ( s( 3 ) * s( 3 ) + s( 4 ) * s( 4 ) + s( 5 ) * s( 5 ) ) );
  };
  Marmot::Math::Vector6d dStrain;
  dStrain << 1.0, -0.3, 0.2, 0.5, -0.4, 0.3;

  {
    auto [material, stateVars] = makePoint();
    MHE::state3D           state;
    Marmot::Math::Matrix6d C;
    state.stateVars = stateVars.data();
    material->computeStress( state, C, dStrain, { 0.0, 1.0 } );
    throwExceptionOnFailure( std::abs( vonMisesStress( state.stress ) / fy - 1 ) < 1e-10,
                             MakeString()
                               << "on the yield surface, sigma_eq / f_y = " << vonMisesStress( state.stress ) / fy );
  }
  // the explicit return mapping (fixed 5 Newton iterations, then the same convergence check)
  {
    auto [material, stateVars] = makePoint();
    MHE::state3D state;
    state.stateVars = stateVars.data();
    material->computeStressExplicit( state, dStrain, { 0.0, 1.0 } );
    throwExceptionOnFailure( std::abs( vonMisesStress( state.stress ) / fy - 1 ) < 1e-10,
                             MakeString() << "explicit: on the yield surface, sigma_eq / f_y = "
                                          << vonMisesStress( state.stress ) / fy );
  }
}

int main()
{
  std::vector< std::function< void( void ) > > tests = { testVonMises,
                                                         testVonMisesCoordinateInvariance,
                                                         testLargeIncrementPerfectPlasticity };

  executeTestsAndCollectExceptions( tests );
  return 0;
}
