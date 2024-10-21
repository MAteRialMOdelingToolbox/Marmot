#include "Marmot/MenetreyWillam.h"
#include "Marmot/MarmotConstants.h"
#include <cmath>
#include <sstream>

namespace Marmot {
  namespace ContinuumMechanics::CommonConstitutiveModels {
    using namespace Constants;
    using namespace ContinuumMechanics::HaighWestergaard;

    MenetreyWillam::MenetreyWillam( const double ft, const MenetreyWillamType& type, const double fc )
    {
      setParameters( ft, fc, type );
    }

    void MenetreyWillam::setParameters( const double ft, const double fc, const MenetreyWillamType& type )
    {
      switch ( type ) {
      case MenetreyWillamType::Mises:
        param.Af = 0;
        param.Bf = sqrt3_2 / ft;
        param.Cf = 0;
        param.m  = 1;
        param.e  = 1;
        break;
      case MenetreyWillamType::Rankine:
        param.Af = 0;
        param.Bf = 1. / ( sqrt6 * ft );
        param.Cf = 1. / ( sqrt3 * ft );
        param.m  = 1;
        param.e  = 0.51;
        break;
      case MenetreyWillamType::DruckerPrager:
        param.Af = 0;
        param.Bf = sqrt3_8 * ( fc + ft ) / ( fc * ft );
        param.Cf = 3. / 2 * ( fc - ft ) / ( fc * ft );
        param.m  = 1;
        param.e  = 1;
        break;
      case MenetreyWillamType::MohrCoulomb:
        param.Af = 0;
        param.Bf = 1. / sqrt6 * ( fc + 2. * ft ) / ( fc * ft );
        param.Cf = 1. / sqrt3 * ( fc - ft ) / ( fc * ft );
        param.m  = 1.;
        param.e  = ( fc + 2 * ft ) / ( 2 * fc + ft );
        break;
      default: throw std::invalid_argument( "Requested MenetreyWillamType not found." );
      }
    }

  } // namespace ContinuumMechanics::CommonConstitutiveModels
} // namespace Marmot
