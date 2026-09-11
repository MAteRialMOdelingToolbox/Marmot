#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMaterialGeneralGradientEnhancedHypoElastic.h"
#include "Marmot/MarmotTesting.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <vector>

using namespace Marmot;
using namespace Marmot::Testing;

namespace {

  /* A material whose stiffness is degraded by the NON-LOCAL FIELD and by nothing else, which is
   * how gradient-enhanced damage models are built: `omega` is a function of the non-local variable
   * handed to computeStress(), not a state the material carries. GCDP is exactly this at its
   * default weighting m = 1.
   *
   * That is what makes the wave-speed query worth a test of its own. The field is an INPUT, so a
   * query that does not pass it does not ask "what is the wave speed now", it asks "what would it
   * be at a field of zero" -- the virgin answer, however damaged the point is.
   */
  class FieldDegradedMaterial : public MarmotMaterialGeneralGradientEnhancedHypoElastic< 1 > {
  public:
    static constexpr double stiffness = 30000.0;
    static constexpr double density_  = 3e-9;

    FieldDegradedMaterial( const double* props, int nProps, int matNumber )
      : MarmotMaterialGeneralGradientEnhancedHypoElastic< 1 >( props, nProps, matNumber )
    {
    }

    void computeStress( response& res, tangents& tan, const increment& inc ) const override
    {
      const double omega = std::min( 0.99, std::max( 0.0, inc.K( 0 ) ) );

      tan.dStressddStrain.setZero();
      for ( int i = 0; i < 6; i++ )
        tan.dStressddStrain( i, i ) = ( 1.0 - omega ) * stiffness;

      res.stress.setZero();
      res.KLocal.setZero();
      res.c.setZero();
    }

    double getDensity( const double* ) const override { return density_; }

    std::vector< double > getNonlocalViscosity( const double* ) const override { return { 1e-4 }; }
  };

  double expectedWaveSpeed( double omega )
  {
    return std::sqrt( ( 1.0 - omega ) * FieldDegradedMaterial::stiffness / FieldDegradedMaterial::density_ );
  }

  FieldDegradedMaterial::response virginResponse( std::vector< double >& stateVars )
  {
    FieldDegradedMaterial::response res;
    res.stress.setZero();
    res.KLocal.setZero();
    res.c.setZero();
    res.stateVars            = stateVars.data();
    res.elasticEnergyDensity = 0.0;
    res.dissipation          = 0.0;
    return res;
  }

  /* Without an explicit field the query must return the UNDAMAGED wave speed. This is the
   * reference the bulk viscosity's optional damage degradation is measured against, so it has to
   * be the virgin one by construction rather than by whenever it happened to be captured -- on a
   * restart, at a refinement, or on an element that enters explicit dynamics already damaged.
   */
  void testDefaultQueryIsTheUndamagedWaveSpeed()
  {
    const std::vector< double > props{ 0.0 };
    FieldDegradedMaterial       mat( props.data(), 1, 1 );
    std::vector< double >       stateVars( mat.getNumberOfRequiredStateVars(), 0.0 );
    auto                        res = virginResponse( stateVars );

    throwExceptionOnFailure( checkIfEqual( mat.getMaximumWaveSpeed( res ), expectedWaveSpeed( 0.0 ), 1e-10 ),
                             MakeString()
                               << __PRETTY_FUNCTION__ << ": the default query is not the undamaged wave speed" );
  }

  /* And with a field it must follow it. This is the regression: the default implementation used to
   * zero inc.K unconditionally, so every query returned the virgin speed, the degradation factor
   * c/c_0 was identically 1.0, and the whole 'bulk viscosity damage degradation' property was
   * inert in any analysis whose damage is driven by the non-local field.
   */
  void testQueryFollowsTheNonlocalField()
  {
    const std::vector< double > props{ 0.0 };
    FieldDegradedMaterial       mat( props.data(), 1, 1 );
    std::vector< double >       stateVars( mat.getNumberOfRequiredStateVars(), 0.0 );
    auto                        res = virginResponse( stateVars );

    const double c0 = mat.getMaximumWaveSpeed( res );

    for ( const double omega : { 0.1, 0.5, 0.75, 0.99 } ) {
      Eigen::Vector< double, 1 > K;
      K( 0 ) = omega;

      const double c = mat.getMaximumWaveSpeed( res, K );

      throwExceptionOnFailure( checkIfEqual( c, expectedWaveSpeed( omega ), 1e-10 ),
                               MakeString() << __PRETTY_FUNCTION__ << ": the wave speed at a non-local field of "
                                            << omega << " is " << c << ", expected " << expectedWaveSpeed( omega ) );

      throwExceptionOnFailure( c < c0,
                               MakeString() << __PRETTY_FUNCTION__ << ": a damaged point reports the virgin wave "
                                            << "speed, so the bulk viscosity's degradation factor would be 1.0" );

      // (c/c_0)^2 is exactly (1 - omega) for this material, which is what the degradation
      // exponent n = 2 reproduces.
      throwExceptionOnFailure( checkIfEqual( ( c / c0 ) * ( c / c0 ), 1.0 - omega, 1e-10 ),
                               MakeString() << __PRETTY_FUNCTION__ << ": (c/c_0)^2 is not 1 - omega" );
    }
  }

  /* The reference and the current query differ in the non-local field and in nothing else, so
   * their ratio measures damage alone -- including on a material whose state variables have been
   * carried in from somewhere else, which a restart does.
   */
  void testReferenceIsIndependentOfTheCarriedState()
  {
    const std::vector< double > props{ 0.0 };
    FieldDegradedMaterial       mat( props.data(), 1, 1 );

    std::vector< double > virginStateVars( mat.getNumberOfRequiredStateVars(), 0.0 );
    std::vector< double > carriedStateVars( mat.getNumberOfRequiredStateVars(), 7.0 );

    auto virginRes  = virginResponse( virginStateVars );
    auto carriedRes = virginResponse( carriedStateVars );

    throwExceptionOnFailure( checkIfEqual( mat.getMaximumWaveSpeed( virginRes ),
                                           mat.getMaximumWaveSpeed( carriedRes ),
                                           1e-10 ),
                             MakeString()
                               << __PRETTY_FUNCTION__ << ": the undamaged reference depends on the carried state" );
  }

} // namespace

int main()
{
  auto tests = std::vector< std::function< void() > >{
    testDefaultQueryIsTheUndamagedWaveSpeed,
    testQueryFollowsTheNonlocalField,
    testReferenceIsIndependentOfTheCarriedState,
  };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
