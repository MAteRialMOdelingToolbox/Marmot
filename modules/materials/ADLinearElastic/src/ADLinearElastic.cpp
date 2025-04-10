#include "Marmot/ADLinearElastic.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotUtility.h"
#include "Marmot/MarmotVoigt.h"
#include "autodiff/forward/dual.hpp"
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
  void ADLinearElastic::computeStressAD( autodiff::dual*       stress,
                                       const autodiff::dual* dStrain,
                                       const double*         timeOld,
                                       const double          dT,
                                       double&               pNewDT )
  {
    using Vector6dual     = Eigen::Matrix< dual, 6, 1 >;
    using mVector6dual   = Eigen::Map< Vector6dual >;
    using mVector6dualConst   = Eigen::Map< const Vector6dual >;



    mVector6dual            s( stress );
    const mVector6dualConst dE( dStrain );

    
    const MatrixXdual C( ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( E, nu ) );

    s = s + C * dE;
  }
} // namespace Marmot::Materials
