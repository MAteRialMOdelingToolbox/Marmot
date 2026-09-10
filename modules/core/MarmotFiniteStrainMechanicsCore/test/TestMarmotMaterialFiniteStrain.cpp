#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include "Marmot/MarmotTesting.h"
#include <cmath>
#include <vector>

using namespace Marmot;
using namespace Marmot::Testing;
using namespace Marmot::FastorStandardTensors;

namespace {

  // A minimal, linear-elastic (isotropic Kirchhoff = C:small-strain), test-only material with a
  // single state variable. No shipped material derived from MarmotMaterialFiniteStrain both (a)
  // has state variables and (b) leaves getMaximumWaveSpeed() at its base-class default, so this is
  // needed to exercise that default's non-null/non-empty state-variable-copying branch.
  class LinearElasticFiniteStrainMaterialWithStateVar : public MarmotMaterialFiniteStrain {
  public:
    LinearElasticFiniteStrainMaterialWithStateVar( const double* matProperties_,
                                                   int           nMaterialProperties_,
                                                   int           materialNumber_ )
      : MarmotMaterialFiniteStrain( matProperties_, nMaterialProperties_, materialNumber_ )
    {
      stateLayout.add( "dummy", 1 );
      stateLayout.finalize();
    }

    void computeStress( ConstitutiveResponse< 3 >& response,
                        AlgorithmicModuli< 3 >&    tangents,
                        const Deformation< 3 >&    deformation,
                        const TimeIncrement& ) const override
    {
      const double& E      = materialProperties[0];
      const double& nu     = materialProperties[1];
      const double  mu     = E / ( 2. * ( 1. + nu ) );
      const double  lambda = E * nu / ( ( 1. + nu ) * ( 1. - 2. * nu ) );

      Tensor33d eps( 0.0 );
      for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ )
          eps( i, j ) = 0.5 * ( deformation.F( i, j ) + deformation.F( j, i ) ) - ( i == j ? 1.0 : 0.0 );

      const double trEps = eps( 0, 0 ) + eps( 1, 1 ) + eps( 2, 2 );

      Tensor33d tau( 0.0 );
      for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ )
          tau( i, j ) = 2. * mu * eps( i, j ) + ( i == j ? lambda * trEps : 0.0 );

      response.tau                  = tau;
      response.elasticEnergyDensity = 0.0;
      response.dissipation          = 0.0;

      tangents.dTau_dF = Tensor3333d( 0.0 );
      for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ )
          for ( int k = 0; k < 3; k++ )
            for ( int l = 0; l < 3; l++ )
              tangents.dTau_dF( i, j, k, l ) = mu * ( ( i == k ? 1. : 0. ) * ( j == l ? 1. : 0. ) +
                                                      ( i == l ? 1. : 0. ) * ( j == k ? 1. : 0. ) ) +
                                               ( i == j && k == l ? lambda : 0.0 );

      if ( response.stateVars != nullptr )
        response.stateVars[0] = trEps;
    }

    double getDensity( const double* ) const override { return materialProperties[2]; }
  };

} // namespace

// ---------------------------------------------------------------------------------------------
// getMaximumWaveSpeed() (base class default): when the material has state variables and a valid
// (non-null) stateVars pointer is supplied, the perturbed evaluations must be performed on
// *copies* of the state, and must reproduce the P-wave modulus sqrt((lambda+2*mu)/rho) for this
// isotropic linear-elastic material.
// ---------------------------------------------------------------------------------------------
void testGetMaximumWaveSpeedCopiesNonEmptyStateVars()
{
  const double                                  E = 20000., nu = 0.25, rho = 2400.;
  const std::vector< double >                   props = { E, nu, rho };
  LinearElasticFiniteStrainMaterialWithStateVar mat( props.data(), 3, 1 );

  std::vector< double > stateVars( mat.getNumberOfRequiredStateVars(), -1.0 );

  MarmotMaterialFiniteStrain::ConstitutiveResponse< 3 > res;
  res.stateVars = stateVars.data();

  const Tensor33d F = Marmot::FastorStandardTensors::Spatial3D::I;

  const double waveSpeed = mat.getMaximumWaveSpeed( stateVars.data(), F );

  const double mu       = E / ( 2. * ( 1. + nu ) );
  const double lambda   = E * nu / ( ( 1. + nu ) * ( 1. - 2. * nu ) );
  const double expected = std::sqrt( ( lambda + 2. * mu ) / rho );

  throwExceptionOnFailure( checkIfEqual( waveSpeed, expected, 1e-5 ),
                           "getMaximumWaveSpeed() does not match sqrt((lambda+2*mu)/rho) for an isotropic "
                           "linear-elastic material in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  // getMaximumWaveSpeed() must operate on a *copy* of the state variables, not mutate the
  // caller's array.
  throwExceptionOnFailure( stateVars[0] == -1.0,
                           "getMaximumWaveSpeed() must not mutate the caller's state variables in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

int main()
{
  const std::vector< std::function< void() > > tests = {
    testGetMaximumWaveSpeedCopiesNonEmptyStateVars,
  };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
