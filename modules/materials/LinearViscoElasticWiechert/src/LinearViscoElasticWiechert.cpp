#include "Marmot/LinearViscoElasticWiechert.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotUtility.h"
#include "Marmot/MarmotViscoelasticity.h"
#include "Marmot/MarmotVoigt.h"
#include "Marmot/MarmotWiechert.h"

#include "Fastor/Fastor.h"
#include "Marmot/MarmotInterfaceMaterialHelperFunctions.h"
#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/Core/util/Constants.h>
#include <Fastor/expressions/linalg_ops/unary_norm_op.h>
#include <Fastor/tensor/TensorMap.h>

#include "autodiff/forward/real.hpp"
#include <iostream>
#include <map>
#include <string>

using namespace Marmot;
using namespace Eigen;

using Tensor1D = Fastor::Tensor< double, 3 >;
using Tensor2D = Fastor::Tensor< double, 3, 3 >;
using Tensor3D = Fastor::Tensor< double, 3, 3, 3 >;
using Tensor4D = Fastor::Tensor< double, 3, 3, 3, 3 >;

namespace Marmot::Materials {

  LinearViscoElasticWiechert::LinearViscoElasticWiechert( const double* materialProperties,
                                                          int           nMaterialProperties,
                                                          int           materialNumber )
    : MarmotMaterialHypoElastic( materialProperties, nMaterialProperties, materialNumber ),
      // clang-format off
      // elasticity parameters
      E( materialProperties[0] ),
      nu( materialProperties[1] ),
      m( materialProperties[2] ),
      n( materialProperties[3] ),
      nMaxwell( static_cast< size_t >( materialProperties[4] ) ),
      minTau( materialProperties[5] ),
      timeToDays( materialProperties[6] )
  // clang-format on
  {
    relaxationTimes = Marmot::Materials::Wiechert::initializeRelaxationTimes( nMaxwell, m );
    elasticModuli   = Marmot::Materials::Wiechert::initializeElasticModuli( nMaxwell, n );

    using namespace Marmot::ContinuumMechanics::Viscoelasticity;
    auto phi_ = [&]( autodiff::Real< powerLawApproximationOrder, double > tau ) {
      return ComplianceFunctions::powerLaw( tau, m, n );
    };

    // elasticModuli =
    // Marmot::Materials::Wiechert::computeElasticModuli<powerLawApproximationOrder>(phi,
    // relaxationTimes);

    zerothWiechertStiffness = 0.0; // m_Ru*(1. - n_Ru )*pow( 2., n_Ru )*pow(minTau_Ru/sqrt(10.), n_Ru);
  }

  void LinearViscoElasticWiechert::computeStress( double*       stress,
                                                  double*       dStressDDStrain,
                                                  const double* dStrain,
                                                  const double* timeOld,
                                                  const double  dT,
                                                  double&       pNewDT )

  {
    mVector6d nomStress( stress );
    Vector6d  dE( dStrain );
    mMatrix6d D( dStressDDStrain );

    if ( ( dE.array() == 0 ).all() && timeOld == 0 ) {
      D = ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( E, nu );
      return;
    }

    Eigen::Ref< Wiechert::mapStateVarMatrix > creepStateVars( stateVarManager->MaxwellStateVars );

    const double dTimeDays = dT * timeToDays;

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

    using namespace Marmot::ContinuumMechanics::Viscoelasticity;
    double effectiveStiffness = E + zerothWiechertStiffness + creepStiffness;
    D                         = ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( effectiveStiffness, nu );
    Vector6d deltaStress      = D * (dE)-creepStressIncrement;
    nomStress                 = nomStress + deltaStress;

    Wiechert::updateStateVarMatrix( dTimeDays, elasticModuli, relaxationTimes, creepStateVars, dE, CelUnit );

    return;
  }

  void LinearViscoElasticWiechert::assignStateVars( double* stateVars_, int nStateVars )
  {
    if ( nStateVars < getNumberOfRequiredStateVars() )
      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": Not sufficient stateVars!" );

    this->stateVarManager = std::make_unique< LinearViscoElasticWiechertStateVarManager >( stateVars_, nMaxwell );

    MarmotMaterial::assignStateVars( stateVars_, nStateVars );
  }

  StateView LinearViscoElasticWiechert::getStateView( const std::string& stateName )
  {
    return stateVarManager->getStateView( stateName );
  }

  int LinearViscoElasticWiechert::getNumberOfRequiredStateVars()
  {
    return LinearViscoElasticWiechertStateVarManager::layout.nRequiredStateVars + 2 * 3 * nMaxwell;
  }
} // namespace Marmot::Materials
