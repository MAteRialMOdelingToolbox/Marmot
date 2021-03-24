#include "Marmot/HaighWestergaard.h"
#include "Marmot/MarmotConstants.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTensor.h"
#include "Marmot/MarmotVoigt.h"
#include <iostream>

using namespace Eigen;

namespace Marmot {
  namespace ContinuumMechanics::HaighWestergaard {

    using namespace Constants;
    using namespace ContinuumMechanics::VoigtNotation;
    using namespace ContinuumMechanics::VoigtNotation::Invariants;

    HaighWestergaardCoordinates haighWestergaard( const Vector6d& stress )
    {
      HaighWestergaardCoordinates hw;
      const double                J2_ = J2( stress );
      hw.xi                           = I1( stress ) / sqrt3;
      hw.rho                          = sqrt( 2 * J2_ );

      if ( hw.rho != 0 ) {
        const double J3_ = J3( stress );
        const double x   = 3 * ( sqrt3 / 2 ) * J3_ / ( pow( J2_, 3. / 2 ) );
        if ( x <= -1 )
          hw.theta = 1. / 3 * Pi;
        else if ( x >= 1 )
          hw.theta = 0;
        else if ( x != x )
          hw.theta = 1. / 3 * Pi;
        else
          hw.theta = 1. / 3 * acos( x );
      }
      else
        hw.theta = 0;

      return hw;
    }

    HaighWestergaardCoordinates haighWestergaardFromStrain( const Vector6d& strain )
    {
      HaighWestergaardCoordinates hw;
      const double                J2_ = J2Strain( strain );

      hw.xi  = I1( strain ) / sqrt3; // will be changed to eM if used
      hw.rho = sqrt( 2 * J2_ );

      if ( hw.rho != 0 ) {
        const double J3_ = J3Strain( strain );
        const double x   = 3. * ( sqrt3 / 2. ) * J3_ / ( pow( J2_, 3. / 2 ) );
        if ( x <= -1 )
          hw.theta = 1. / 3 * Constants::Pi;
        else if ( x >= 1 )
          hw.theta = 0;
        else if ( x != x )
          hw.theta = 1. / 3 * Constants::Pi;
        else
          hw.theta = 1. / 3 * acos( x );
      }
      else
        hw.theta = 0;

      return hw;
    }

    namespace Complex {

      // implementation for complex data types
      using namespace Constants;
      using namespace ContinuumMechanics::VoigtNotation::Invariants::Complex;

      HaighWestergaardCoordinates haighWestergaard( const Vector6cd& stress )
      {
        HaighWestergaardCoordinates hw;
        const complexDouble         J2_ = J2( stress );
        hw.xi                           = I1( stress ) / sqrt3;
        hw.rho                          = sqrt( 2. * J2_ );

        if ( hw.rho.real() != 0. ) {
          const complexDouble J3_       = J3( stress );
          const complexDouble cos3Theta = ( 3. * ( sqrt3 / 2. ) * J3_ / ( pow( J2_, 3. / 2 ) ) );
          if ( cos3Theta.real() < -1 )
            hw.theta = 1. / 3 * Pi;
          else if ( cos3Theta.real() > 1 )
            hw.theta = 0.;
          else if ( cos3Theta != cos3Theta )
            hw.theta = 1. / 3 * Pi;
          else
            hw.theta = 1. / 3 * acos( cos3Theta );
        }
        else
          hw.theta = 0.;

        return hw;
      }

    } // namespace Complex

  } // namespace ContinuumMechanics::HaighWestergaard
} // namespace Marmot
