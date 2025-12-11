#include "Marmot/ADLinearElastic.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotTypedefs.h"
#include "autodiff/forward/dual/eigen.hpp"

namespace Marmot::Materials {

  using namespace Marmot;
  using namespace autodiff;
  using namespace Eigen;
  using namespace ContinuumMechanics::Elasticity;

  ADLinearElastic::ADLinearElastic( const double* materialProperties, int nMaterialProperties, int materialNumber )
    : MarmotMaterialHypoElasticAD::MarmotMaterialHypoElasticAD( materialProperties,
                                                                nMaterialProperties,
                                                                materialNumber ),
      E( materialProperties[0] ),
      nu( materialProperties[1] )
  {
    assert( nMaterialProperties == 2 );
  }
  void ADLinearElastic::computeStressAD( state3DAD&            state,
                                         const autodiff::dual* dStrain,
                                         const timeInfo&       timeInfo ) const
  {
    mVector6dual            s( state.stress );
    const mVector6dualConst dE( dStrain );

    const MatrixXdual C( ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( E, nu ) );

    s = s + C * dE;
  }
} // namespace Marmot::Materials
