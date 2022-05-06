#include "Marmot/B4Shrinkage.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotUtility.h"
#include "Marmot/MarmotVoigt.h"
#include <iostream>

using namespace Marmot;
using namespace std;

namespace Marmot::Materials {
  namespace Shrinkage::B4 {
    Vector6d computeShrinkageStrainIncrement( const double tStartDays,
                                              const double dTDays,
                                              const double ultimateAutogenousShrinkageStrain,
                                              const double autogenousShrinkageHalfTime,
                                              const double alpha,
                                              const double rt,
                                              const double ultimateDryingShrinkageStrain,
                                              const double dryingShrinkageHalfTime,
                                              const double kHum,
                                              const double dryingStart )
    {

      const double tEndDays = tStartDays + dTDays;
      if ( tEndDays == 0 )
        return Vector6d::Zero();

      // autogenous shrinkage
      double deltaAutogenous = 0;

      if ( tEndDays > 5e-2 ) {
        deltaAutogenous = pow( 1. + pow( autogenousShrinkageHalfTime / ( tStartDays + dTDays ), alpha ), rt ) -
                          pow( 1. + pow( autogenousShrinkageHalfTime / tStartDays, alpha ), rt );
      }

      // drying shrinkage

      const double end         = Math::macauly( tEndDays - dryingStart );
      double       deltaDrying = 0;
      double       start       = Math::macauly( tStartDays - dryingStart );

      if ( end != 0 )
        deltaDrying = tanh( sqrt( end / dryingShrinkageHalfTime ) ) - tanh( sqrt( start / dryingShrinkageHalfTime ) );

      return ContinuumMechanics::VoigtNotation::I * ( ultimateAutogenousShrinkageStrain * deltaAutogenous +
                                                      ultimateDryingShrinkageStrain * kHum * deltaDrying );
    }
  } // namespace Shrinkage::B4
} // namespace Marmot::Materials
