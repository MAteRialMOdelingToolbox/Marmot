#include "Marmot/LinearViscoElasticWiechert.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotWiechert.h"

#include <Eigen/Core>

#include <stdexcept>

namespace Marmot::Materials {

  using namespace Eigen;
  using namespace Marmot;

  LinearViscoElasticWiechert::LinearViscoElasticWiechert( const double* materialProperties,
                                                          int           nMaterialProperties,
                                                          int           materialNumber )
    : MarmotMaterialHypoElastic( materialProperties, nMaterialProperties, materialNumber ),
      // clang-format off
      // elasticity parameters
      E( materialProperties[0] ),
      nu( materialProperties[1] ),
      // viscoelastic parameters for 1 Maxwell element, relaxation time
      m( materialProperties[2] ),
      // viscoelastic parameters for 1 Maxwell element, elastic modulus
      n( materialProperties[3] ),
      nMaxwell( static_cast< size_t >( materialProperties[4] ) ),
      minTau( materialProperties[5] ),
      timeToDays( materialProperties[6] ),
      relaxationTimes( static_cast< int >( nMaxwell ) ),
      elasticModuli( static_cast< int >( nMaxwell ) ),
      zerothWiechertStiffness( 0.0 )
  // clang-format on
  {
    initializeStateLayout();
    relaxationTimes = Marmot::Materials::Wiechert::initializeRelaxationTimes( nMaxwell, m );
    elasticModuli   = Marmot::Materials::Wiechert::initializeElasticModuli( nMaxwell, n );
  }

  void LinearViscoElasticWiechert::computeStress( state3D&        state,
                                                  Matrix6d&       dStressDDStrain,
                                                  const Vector6d& dStrain,
                                                  const timeInfo& timeInfo ) const
  {
    mVector6d     nomStress( state.stress.data() );
    mMatrix6d     D( dStressDDStrain.data() );
    const auto    dE = Map< const Vector6d >( dStrain.data() );
    const double& dT = timeInfo.dT;

    if ( ( dE.array() == 0 ).all() && dT == 0 ) {
      D = ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( E, nu );
      return;
    }

    Eigen::Map< Wiechert::StateVarMatrix > creepStateVars( state.stateVars, 6, nMaxwell );
    const double                           dTimeDays = dT * timeToDays;

    Matrix6d CelUnit = ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( 1.0, nu );

    Vector6d creepStressIncrement = Vector6d::Zero();
    double   creepStiffness       = 0;

    Wiechert::evaluateWiechert( dTimeDays,
                                elasticModuli,
                                relaxationTimes,
                                creepStateVars,
                                creepStiffness,
                                creepStressIncrement,
                                1.0 );

    double effectiveStiffness = E + zerothWiechertStiffness + creepStiffness;

    D                    = ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( effectiveStiffness, nu );
    Vector6d deltaStress = D * dE - creepStressIncrement;
    nomStress            = nomStress + deltaStress;

    Wiechert::updateStateVarMatrix( dTimeDays, elasticModuli, relaxationTimes, creepStateVars, deltaStress, CelUnit );

    return;
  }

  double LinearViscoElasticWiechert::getDensity( const double* stateVars ) const
  {
    if ( nMaterialProperties <= 7 + 1 ) {
      throw std::runtime_error(
        std::string( MakeString() << __PRETTY_FUNCTION__ << ": Density not specified for this material." ) );
    }
    return this->materialProperties[7];
  }
} // namespace Marmot::Materials
