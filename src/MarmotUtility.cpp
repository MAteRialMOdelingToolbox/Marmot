#include "Marmot/MarmotUtility.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"

namespace Marmot {

  void discardTheIncrement( double& pNewDT, double value, const std::string& message )
  {
    pNewDT = value;
    MarmotJournal::warningToMSG( message );
    return;
  }

} // namespace Marmot
