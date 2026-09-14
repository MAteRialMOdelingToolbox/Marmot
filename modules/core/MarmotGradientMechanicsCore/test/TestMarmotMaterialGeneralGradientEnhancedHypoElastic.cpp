#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMaterialGeneralGradientEnhancedHypoElastic.h"
#include "Marmot/MarmotTesting.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace Marmot;
using namespace Marmot::Testing;

namespace {

  /* A material degraded by the NON-LOCAL FIELD alone, which is how gradient-enhanced damage
   * models are built -- `omega` is a function of the variable handed to computeStress(), not a
   * state the material carries, and GCDP is exactly this at its default m = 1. The field being an
   * INPUT is what makes the query worth a test: without it, the answer is the virgin one however
   * damaged the point is.
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

  /* A material that reads a micro-inertia out of its properties and returns it through the base
   * class's validation, which is how every material in this stack is expected to answer.
   */
  class MicroInertiaMaterial : public MarmotMaterialGeneralGradientEnhancedHypoElastic< 1 > {
  public:
    static constexpr double viscosity = 1e-4;

    MicroInertiaMaterial( const double* props, int nProps, int matNumber )
      : MarmotMaterialGeneralGradientEnhancedHypoElastic< 1 >( props, nProps, matNumber )
    {
    }

    void computeStress( response&, tangents&, const increment& ) const override {}

    double getDensity( const double* ) const override { return 3e-9; }

    std::vector< double > getNonlocalViscosity( const double* ) const override { return { viscosity }; }

    std::vector< double > getNonlocalMicroInertia( const double* stateVars ) const override
    {
      return validatedNonlocalMicroInertia( { materialProperties[0] }, stateVars );
    }
  };

  /* The micro-inertia is zero unless a material provides one, and that is not merely a default but
   * the whole compatibility story: every material written before this interface existed keeps a
   * first-order non-local field, integrated by its viscosity alone, without being touched.
   */
  void testMicroInertiaIsZeroForAMaterialThatProvidesNone()
  {
    const std::vector< double > props{ 0.0 };
    FieldDegradedMaterial       mat( props.data(), 1, 1 );

    const auto m_k = mat.getNonlocalMicroInertia( nullptr );

    throwExceptionOnFailure( m_k.size() == 1,
                             MakeString()
                               << __PRETTY_FUNCTION__ << ": the default must return one entry per non-local variable" );
    throwExceptionOnFailure( m_k[0] == 0.0,
                             MakeString()
                               << __PRETTY_FUNCTION__ << ": the default micro-inertia is " << m_k[0] << ", not zero" );
  }

  /* m_k and eta are ONE parameter: above eta^2/4 the zeroth-order reaction mode is underdamped and
   * the regularisation rings, which looks exactly like the mesh-scale oscillation the gradient
   * enhancement exists to remove. Nothing could check this while the micro-inertia was an element
   * property and the viscosity a material one -- the two never met in the same object.
   */
  void testMicroInertiaIsBoundedByTheViscosity()
  {
    const double admissible = 0.25 * MicroInertiaMaterial::viscosity * MicroInertiaMaterial::viscosity;

    const auto accepts = []( double m_k ) {
      const std::vector< double > props{ m_k };
      MicroInertiaMaterial        mat( props.data(), 1, 1 );
      try {
        return mat.getNonlocalMicroInertia( nullptr )[0] == m_k;
      }
      catch ( const std::invalid_argument& ) {
        return false;
      }
    };

    // The recommended value is the largest admissible one, so the bound has to accept it exactly.
    throwExceptionOnFailure( accepts( admissible ),
                             MakeString() << __PRETTY_FUNCTION__ << ": eta^2/4 = " << admissible
                                          << " was rejected, which is the value the documentation recommends" );
    throwExceptionOnFailure( accepts( 0.5 * admissible ),
                             MakeString() << __PRETTY_FUNCTION__ << ": an admissible micro-inertia was rejected" );
    throwExceptionOnFailure( accepts( 0.0 ),
                             MakeString() << __PRETTY_FUNCTION__ << ": a zero micro-inertia was rejected" );

    throwExceptionOnFailure( !accepts( 2.0 * admissible ),
                             MakeString() << __PRETTY_FUNCTION__
                                          << ": a micro-inertia of twice eta^2/4 was accepted; the reaction mode "
                                             "rings there" );
    throwExceptionOnFailure( !accepts( -admissible ),
                             MakeString() << __PRETTY_FUNCTION__ << ": a negative micro-inertia was accepted" );
    throwExceptionOnFailure( !accepts( std::numeric_limits< double >::quiet_NaN() ),
                             MakeString() << __PRETTY_FUNCTION__ << ": a non-finite micro-inertia was accepted" );
  }

} // namespace

int main()
{
  auto tests = std::vector< std::function< void() > >{
    testDefaultQueryIsTheUndamagedWaveSpeed,
    testQueryFollowsTheNonlocalField,
    testReferenceIsIndependentOfTheCarriedState,
    testMicroInertiaIsZeroForAMaterialThatProvidesNone,
    testMicroInertiaIsBoundedByTheViscosity,
  };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
