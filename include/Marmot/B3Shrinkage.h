#include "Marmot/MarmotTypedefs.h"

namespace Marmot::Materials {
  namespace Shrinkage {
    namespace B3 {

      Marmot::Vector6d computeShrinkageStrainIncrement( const double tStartDays,
                                                        const double dTDays,
                                                        const double ultimateShrinkageStrain,
                                                        const double shrinkageHalfTime,
                                                        const double kHum );
    }
  } // namespace Shrinkage
} // namespace Marmot::Materials
