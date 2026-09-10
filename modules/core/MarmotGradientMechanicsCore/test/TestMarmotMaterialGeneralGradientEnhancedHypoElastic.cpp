#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotMaterialGeneralGradientEnhancedHypoElastic.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotTypedefs.h"
#include <Eigen/Dense>
#include <utility>
#include <vector>

using namespace Marmot;
using namespace Marmot::Testing;

namespace {

  using Mat1 = MarmotMaterialGeneralGradientEnhancedHypoElastic< 1 >;

  // ---------------------------------------------------------------------------------------------
  // A minimal test-only material whose out-of-plane stress component (index 2) is a positive
  // constant, independent of strain, with a correspondingly zero tangent for that component. No
  // real, registered material behaves like this (dStress/dStrain is always non-singular for a
  // physically sound hypoelastic law), so this stand-in is the only practical way to deterministically
  // exercise computePlaneStress()'s robustness branches: the near-singular-tangent compliance cap,
  // and the cutback/failure throw when the iteration cannot drive sigma_zz to zero.
  // ---------------------------------------------------------------------------------------------
  class NeverConvergingPlaneStressMaterial : public Mat1 {
  public:
    using Mat1::Mat1;

    void initializeStateLayout()
    {
      stateLayout.add( "dummy", 1 );
      stateLayout.finalize();
    }

    void computeStress( response& res, tangents& tan, const increment& inc ) const override
    {
      res.stress      = Marmot::Vector6d::Zero();
      res.stress( 2 ) = 5.0; // never zero, and independent of inc.dStrain
      res.KLocal( 0 ) = 0.0;
      res.c( 0 )      = 1.0;

      tan.dStressddStrain.setZero(); // dSigma_zz/dEps_zz == 0 -> forces the near-singular branch
      tan.dStressddK.setZero();
      tan.dKLocalddStrain.setZero();
      tan.dKLocalddK.setZero();
      tan.dcddK.setZero();
      tan.d2cddK2.setZero();
    }

    double getDensity( const double* ) const override { return 0.0; }

    std::vector< double > getNonlocalViscosity( const double* ) const override { return { 0.0 }; }
  };

  std::pair< NeverConvergingPlaneStressMaterial, std::vector< double > > makeNeverConvergingMaterial()
  {
    NeverConvergingPlaneStressMaterial mat( nullptr, 0, 1 );
    mat.initializeStateLayout();
    std::vector< double > stateVars( mat.getNumberOfRequiredStateVars(), 0.0 );
    return { std::move( mat ), std::move( stateVars ) };
  }

} // namespace

// ---------------------------------------------------------------------------------------------
// computePlaneStress(): must throw StressUpdateFailed after exhausting its cutback attempts when
// the out-of-plane stress can never be driven to zero. Along the way, the near-singular tangent
// (dStressddStrain(2,2) == 0) must hit the "cap the compliance instead of dividing by zero" branch
// rather than propagating a NaN/inf strain correction.
// ---------------------------------------------------------------------------------------------
void testComputePlaneStressThrowsWhenItCannotConverge()
{
  auto [mat, stateVars] = makeNeverConvergingMaterial();

  Mat1::increment inc;
  inc.dStrain      = Marmot::Vector6d::Zero();
  inc.dStrain( 0 ) = 1e-3;
  inc.K( 0 )       = 0.0;
  inc.dK( 0 )      = 0.0;
  inc.time         = 0.0;
  inc.dT           = 1.0;

  Mat1::response res;
  Mat1::tangents tan;
  res.stress    = Marmot::Vector6d::Zero();
  res.KLocal    = Eigen::Vector< double, 1 >::Zero();
  res.c         = Eigen::Vector< double, 1 >::Zero();
  res.stateVars = stateVars.data();

  bool threw = false;
  try {
    mat.computePlaneStress( res, tan, inc );
  }
  catch ( const Marmot::StressUpdateFailed& ) {
    threw = true;
  }
  throwExceptionOnFailure( threw,
                           "computePlaneStress() must throw StressUpdateFailed when it cannot converge in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// ---------------------------------------------------------------------------------------------
// getMaximumWaveSpeed(): must return 0 (rather than dividing by zero / returning NaN) when the
// material reports zero density.
// ---------------------------------------------------------------------------------------------
void testGetMaximumWaveSpeedReturnsZeroForZeroDensity()
{
  auto [mat, stateVars] = makeNeverConvergingMaterial();

  Mat1::response res;
  res.stress    = Marmot::Vector6d::Zero();
  res.KLocal    = Eigen::Vector< double, 1 >::Zero();
  res.c         = Eigen::Vector< double, 1 >::Zero();
  res.stateVars = stateVars.data();

  const double waveSpeed = mat.getMaximumWaveSpeed( res );
  throwExceptionOnFailure( waveSpeed == 0.0,
                           "getMaximumWaveSpeed() must return 0 for zero density in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// ---------------------------------------------------------------------------------------------
// getMaximumWaveSpeed(): must tolerate a null state-variable pointer / zero state-variable count
// (e.g. a material with no state variables at all) rather than attempting to copy through it.
// ---------------------------------------------------------------------------------------------
void testGetMaximumWaveSpeedToleratesNullStateVars()
{
  NeverConvergingPlaneStressMaterial mat( nullptr, 0, 1 );
  // Deliberately skip initializeStateLayout(): zero state variables, matching a null stateVars
  // pointer in the response, which is the scenario the null-pointer guard in getMaximumWaveSpeed()
  // is there to handle.
  mat.stateLayout.finalize();

  Mat1::response res;
  res.stress    = Marmot::Vector6d::Zero();
  res.KLocal    = Eigen::Vector< double, 1 >::Zero();
  res.c         = Eigen::Vector< double, 1 >::Zero();
  res.stateVars = nullptr;

  const double waveSpeed = mat.getMaximumWaveSpeed( res );
  throwExceptionOnFailure( waveSpeed == 0.0,
                           "getMaximumWaveSpeed() must tolerate a null state-variable pointer in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// ---------------------------------------------------------------------------------------------
int main()
{
  const std::vector< std::function< void() > > tests = {
    testComputePlaneStressThrowsWhenItCannotConverge,
    testGetMaximumWaveSpeedReturnsZeroForZeroDensity,
    testGetMaximumWaveSpeedToleratesNullStateVars,
  };

  executeTestsAndCollectExceptions( tests );
  return 0;
}
