#include "Marmot/ADLinearElastic.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMaterialPointSolverHypoElastic.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotTypedefs.h"

using namespace Marmot::Testing;
using namespace Marmot::Materials;

auto initializeSolverADLinearElastic()
{
  std::vector< double > materialProperties = { 210000., 0.3 };
  auto                  solveropts         = Marmot::Materials::MaterialPointSolverHypoElastic::SolverOptions();
  std::string           matName            = "ADLINEARELASTIC";
  return std::unique_ptr< Marmot::Materials::MaterialPointSolverHypoElastic >(
    new Marmot::Materials::MaterialPointSolverHypoElastic( matName,
                                                           &materialProperties[0],
                                                           static_cast< int >( materialProperties.size() ),
                                                           solveropts ) );
}

void testADLinearElasticObjectivity()
{

  std::vector< double > materialProperties = { 210000., 0.3 };
  auto                  solveropts         = Marmot::Materials::MaterialPointSolverHypoElastic::SolverOptions();
  std::string           matName            = "ADLINEARELASTIC";
  auto                  solver             = Marmot::Materials::MaterialPointSolverHypoElastic( matName,
                                                                   &materialProperties[0],
                                                                   static_cast< int >( materialProperties.size() ),
                                                                   solveropts );

  Marmot::Materials::MaterialPointSolverHypoElastic::Step step;

  step.isStrainComponentControlled = { true, true, true, true, true, true };
  step.isStressComponentControlled = { false, false, false, false, false, false };
  step.stressIncrementTarget       = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  step.strainIncrementTarget       = { 0.00839244, 0.00089344, -0.00703916, 0.00013635, 0.00160548, 0.00572825 };

  step.dTStart = 1.0;

  solver.addStep( step );
  // check coordinate invariance
  throwExceptionOnFailure( spinTurbokreisel( solver, 1e-10, 1e-8 ),
                           Marmot::Utilities::MakeString() << __PRETTY_FUNCTION__ << "Turbokreisel failed!" );
}

void testADLinearElastic()
{

  std::vector< double > materialProperties = { 210000., 0.3 };
  auto                  solveropts         = Marmot::Materials::MaterialPointSolverHypoElastic::SolverOptions();
  std::string           matName            = "ADLINEARELASTIC";
  auto                  solver             = Marmot::Materials::MaterialPointSolverHypoElastic( matName,
                                                                   &materialProperties[0],
                                                                   static_cast< int >( materialProperties.size() ),
                                                                   solveropts );

  Marmot::Materials::MaterialPointSolverHypoElastic::Step step;

  step.isStrainComponentControlled = { true, true, true, true, true, true };
  step.isStressComponentControlled = { false, false, false, false, false, false };
  step.stressIncrementTarget       = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  step.strainIncrementTarget       = { 0.00839244, 0.00089344, -0.00703916, 0.00013635, 0.00160548, 0.00572825 };

  step.dTStart = 1.0;

  solver.addStep( step );

  // solve
  solver.solve();

  // reference solution
  const Marmot::Math::Vector6d stressTarget =
    { 1627.90061538, 416.523692308, -864.896307692, 11.0128846154, 129.673384615, 462.666346154 };

  auto history = solver.getHistory();

  // compare solutions
  throwExceptionOnFailure( checkIfEqual< double >( history.back().stress, stressTarget, 1e-8 ),
                           Marmot::Utilities::MakeString()
                             << __PRETTY_FUNCTION__ << "comparison with reference solution failed" );
}

void testADLinearElasticEnergyDensityOfNonProportionalStrainPath()
{
  // The strain energy of linear elasticity depends only on the final strain, W = 1/2 eps : C : eps
  // = 1/2 ( lambda tr(eps)^2 + 2 mu eps : eps ), whatever the path. A non-proportional path in several increments
  // checks the accumulation from a non-zero stress at the start of each increment.
  std::vector< double > materialProperties = { 210000., 0.3 };
  ADLinearElastic       material( materialProperties.data(), static_cast< int >( materialProperties.size() ), 1 );
  MaterialHypoElastic&  materialBase = material;
  const double          E = materialProperties[0], nu = materialProperties[1];
  const double          lambda = E * nu / ( ( 1. + nu ) * ( 1. - 2. * nu ) );
  const double          mu     = E / ( 2. * ( 1. + nu ) );
  MaterialHypoElastic::timeInfo timeInfo{ 0., 1. };

  Marmot::Math::Vector6d legA, legB; // engineering shear strains
  legA << 1e-3, -2e-4, 3e-4, 4e-4, 0., -1e-4;
  legB << -5e-4, 6e-4, 0., -2e-4, 5e-4, 3e-4;

  MaterialHypoElastic::state3D state;
  Marmot::Math::Matrix6d       tangent;
  Marmot::Math::Vector6d       strain = Marmot::Math::Vector6d::Zero();
  for ( const auto& [leg, nIncrements] : { std::pair{ legA, 4 }, std::pair{ legB, 7 } } )
    for ( int i = 0; i < nIncrements; i++ ) {
      const Marmot::Math::Vector6d dStrain = leg / nIncrements;
      materialBase.computeStress( state, tangent, dStrain, timeInfo );
      strain += dStrain;
    }

  const double trace           = strain.head< 3 >().sum();
  const double strainDotStrain = strain.head< 3 >().squaredNorm() + 0.5 * strain.tail< 3 >().squaredNorm();
  const double expected        = 0.5 * ( lambda * trace * trace + 2. * mu * strainDotStrain );

  throwExceptionOnFailure( checkIfEqual( state.elasticEnergyDensity, expected, 1e-10 ),
                           Marmot::Utilities::MakeString()
                             << __PRETTY_FUNCTION__
                             << ": elastic energy density must be 1/2 eps : C : eps after a non-proportional path" );
  throwExceptionOnFailure( checkIfEqual( state.dissipation, 0.0, 1e-15 ),
                           Marmot::Utilities::MakeString()
                             << __PRETTY_FUNCTION__ << ": linear elasticity must not dissipate" );
}

int main()
{

  std::vector< std::function< void( void ) > > tests = {
    testADLinearElastic,
    testADLinearElasticObjectivity,
    testADLinearElasticEnergyDensityOfNonProportionalStrainPath,
  };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
