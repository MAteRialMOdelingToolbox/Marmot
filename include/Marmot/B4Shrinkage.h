#include "Marmot/MarmotTypedefs.h"

namespace Marmot::Materials {
  namespace Shrinkage {
    namespace B4 {

      Marmot::Vector6d computeShrinkageStrainIncrement( const double tStartDays,
                                                        const double dTDays,
                                                        const double ultimateAutogenousShrinkageStrain,
                                                        const double autogenousShrinkageHalfTime,
                                                        const double alpha,
                                                        const double rt,
                                                        const double ultimateDryingShrinkageStrain,
                                                        const double dryingShrinkageHalfTime,
                                                        const double kHum,
                                                        const double dryingStart );
    }
  } // namespace Shrinkage
} // namespace Marmot::Materials
