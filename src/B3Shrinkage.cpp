#include "Marmot/B3Shrinkage.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotUtility.h"
#include "Marmot/MarmotVoigt.h"
#include <iostream>

using namespace Marmot;
using namespace std;

namespace Marmot::Materials {
  namespace Shrinkage::B3 {
    Vector6d computeShrinkageStrainIncrement( const double tStartDays,
                                              const double dTDays,
                                              const double ultimateShrinkageStrain,
                                              const double shrinkageHalfTime,
                                              const double kHum )
    {
      const double tEndDays = tStartDays + dTDays;
      double       delta    = 0;

      if ( tEndDays < 5e-2 )
        return Vector6d::Zero();
      delta = tanh( sqrt( ( tStartDays + dTDays ) / shrinkageHalfTime ) ) -
              tanh( sqrt( tStartDays / shrinkageHalfTime ) );

      return ContinuumMechanics::VoigtNotation::I * ultimateShrinkageStrain * kHum * delta;
    }
  } // namespace Shrinkage::B3
} // namespace Marmot::Materials
