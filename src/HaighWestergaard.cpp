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

        Vector3d haighWestergaard( const Vector6d& stress )
        {
            Vector3d     hw;
            const double J2_ = J2( stress );
            hw( 0 )          = I1( stress ) / sqrt3;
            hw( 1 )          = sqrt( 2 * J2_ );

            if ( hw( 1 ) != 0 ) {
                const double J3_ = J3( stress );
                const double x   = 3 * ( sqrt3 / 2 ) * J3_ / ( pow( J2_, 3. / 2 ) );
                if ( x <= -1 )
                    hw( 2 ) = 1. / 3 * Pi;
                else if ( x >= 1 )
                    hw( 2 ) = 0;
                else if ( x != x )
                    hw( 2 ) = 1. / 3 * Pi;
                else
                    hw( 2 ) = 1. / 3 * acos( x );
            }
            else
                hw( 2 ) = 0;

            return hw;
        }

        Vector3d haighWestergaardStrain( const Vector6d& strain )
        {
            Vector3d     hw;
            const double J2_ = J2strain( strain );

            hw( 0 ) = I1( strain ) / sqrt3; // will be changed to eM if used
            hw( 1 ) = sqrt( 2 * J2_ );

            if ( hw( 1 ) != 0 ) {
                const double J3_ = J3strain( strain );
                const double x   = 3. * ( sqrt3 / 2. ) * J3_ / ( pow( J2_, 3. / 2 ) );
                if ( x <= -1 )
                    hw( 2 ) = 1. / 3 * Constants::Pi;
                else if ( x >= 1 )
                    hw( 2 ) = 0;
                else if ( x != x )
                    hw( 2 ) = 1. / 3 * Constants::Pi;
                else
                    hw( 2 ) = 1. / 3 * acos( x );
            }
            else
                hw( 2 ) = 0;

            return hw;
        }

        Vector3d principalStrainsHW( const Vector6d& voigtStrain )
        {
            // if you wanna sort your eigenvalues after size
            Vector3d hw = haighWestergaardStrain( voigtStrain );
            Vector3d strainPrinc;

            strainPrinc << hw( 0 ) / sqrt3 + sqrt2_3 * hw( 1 ) * std::cos( hw( 2 ) ),
                hw( 0 ) / sqrt3 + sqrt2_3 * hw( 1 ) * ( -std::sin( Constants::Pi / 6. - hw( 2 ) ) ),
                hw( 0 ) / sqrt3 + sqrt2_3 * hw( 1 ) * ( -std::sin( Constants::Pi / 6. + hw( 2 ) ) );

            return strainPrinc;
        }

    } // namespace ContinuumMechanics::HaighWestergaard

} // namespace Marmot
