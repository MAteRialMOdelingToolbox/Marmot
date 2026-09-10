#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotTesting.h"
#include <utility>
#include <vector>

using namespace Marmot;
using namespace Marmot::Testing;

namespace {

  // A minimal test-only material whose out-of-plane/lateral stress components are positive
  // constants, independent of strain, with correspondingly zero tangent entries for those
  // components. No real, registered hypoelastic material behaves like this (dStress/dStrain is
  // always non-singular for a physically sound law), so this is the only practical way to
  // deterministically exercise computePlaneStress()'s and computeUniaxialStress()'s robustness
  // branches: the near-singular-tangent compliance cap, and the cutback/failure throw.
  class NeverConvergingHypoElasticMaterial : public MarmotMaterialHypoElastic {
  public:
    NeverConvergingHypoElasticMaterial() : MarmotMaterialHypoElastic( nullptr, 0, 1 )
    {
      stateLayout.add( "dummy", 1 );
      stateLayout.finalize();
    }

    void computeStress( state3D&          state,
                        Marmot::Matrix6d& dStress_dStrain,
                        const Marmot::Vector6d&,
                        const timeInfo& ) const override
    {
      state.stress      = Marmot::Vector6d::Zero();
      state.stress( 1 ) = 5.0; // never zero, and independent of dStrain
      state.stress( 2 ) = 5.0;

      dStress_dStrain.setZero(); // dSigma_yy/dEps_yy == dSigma_zz/dEps_zz == 0
    }

    double getDensity( const double* ) const override { return 1.0; }
  };

  std::pair< NeverConvergingHypoElasticMaterial, std::vector< double > > makeNeverConvergingMaterial()
  {
    NeverConvergingHypoElasticMaterial mat;
    std::vector< double >              stateVars( mat.getNumberOfRequiredStateVars(), 0.0 );
    return { std::move( mat ), std::move( stateVars ) };
  }

} // namespace

// ---------------------------------------------------------------------------------------------
// computePlaneStress(): must throw StressUpdateFailed after exhausting its cutback attempts when
// the out-of-plane stress can never be driven to zero, hitting the near-singular tangent
// (dStress_dStrain(2,2) == 0) "cap the compliance instead of dividing by zero" branch first.
// ---------------------------------------------------------------------------------------------
void testComputePlaneStressThrowsWhenItCannotConverge()
{
  auto [mat, stateVars] = makeNeverConvergingMaterial();

  Marmot::Vector3d dStrain2D = Marmot::Vector3d::Zero();
  dStrain2D( 0 )             = 1e-3;
  MarmotMaterialHypoElastic::timeInfo timeInfo{ 0.0, 1.0 };
  MarmotMaterialHypoElastic::state2D  state2D;
  state2D.stateVars = stateVars.data();
  Marmot::Matrix3d tan2D;

  bool threw = false;
  try {
    mat.computePlaneStress( state2D, tan2D, dStrain2D, timeInfo );
  }
  catch ( const Marmot::StressUpdateFailed& ) {
    threw = true;
  }
  throwExceptionOnFailure( threw,
                           "computePlaneStress() must throw StressUpdateFailed when it cannot converge in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// ---------------------------------------------------------------------------------------------
// computeUniaxialStress(): same robustness requirement as computePlaneStress(), but condensing
// both lateral stress components (indices 1 and 2) simultaneously via a 2x2 solve.
// ---------------------------------------------------------------------------------------------
void testComputeUniaxialStressThrowsWhenItCannotConverge()
{
  auto [mat, stateVars] = makeNeverConvergingMaterial();

  const double                        dStrain1D = 1e-3;
  MarmotMaterialHypoElastic::timeInfo timeInfo{ 0.0, 1.0 };
  MarmotMaterialHypoElastic::state1D  state1D;
  state1D.stateVars = stateVars.data();
  double tan1D;

  bool threw = false;
  try {
    mat.computeUniaxialStress( state1D, tan1D, dStrain1D, timeInfo );
  }
  catch ( const Marmot::StressUpdateFailed& ) {
    threw = true;
  }
  throwExceptionOnFailure( threw,
                           "computeUniaxialStress() must throw StressUpdateFailed when it cannot converge in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// ---------------------------------------------------------------------------------------------
// getMaximumWaveSpeed(): when the material has state variables and a valid (non-null) stateVars
// pointer is supplied, the perturbed evaluation must be performed on a *copy* of the state.
// ---------------------------------------------------------------------------------------------
void testGetMaximumWaveSpeedCopiesNonEmptyStateVars()
{
  auto [mat, stateVars] = makeNeverConvergingMaterial();
  stateVars[0]          = -1.0;

  MarmotMaterialHypoElastic::state3D state;
  state.stateVars = stateVars.data();

  const double waveSpeed = mat.getMaximumWaveSpeed( state );

  // dStress_dStrain is identically zero for this material, so the max stiffness diagonal is 0.
  throwExceptionOnFailure( waveSpeed == 0.0,
                           "getMaximumWaveSpeed() did not return 0 for a zero tangent in " +
                             std::string( __PRETTY_FUNCTION__ ) );
  throwExceptionOnFailure( stateVars[0] == -1.0,
                           "getMaximumWaveSpeed() must not mutate the caller's state variables in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

int main()
{
  const std::vector< std::function< void() > > tests = {
    testComputePlaneStressThrowsWhenItCannotConverge,
    testComputeUniaxialStressThrowsWhenItCannotConverge,
    testGetMaximumWaveSpeedCopiesNonEmptyStateVars,
  };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
