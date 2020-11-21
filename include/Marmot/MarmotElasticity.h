#pragma once

namespace Marmot {

    namespace mechanics {
        namespace Elasticity {

            double constexpr shearModulus( const double E, const double nu ) { return E / ( 2 * ( 1 + nu ) ); }

            double constexpr lameParameter( const double E, const double nu )
            {
                return E * nu / ( ( 1 + nu ) * ( 1 - 2 * nu ) );
            }

        } // namespace Elasticity

    } // namespace mechanics

} // namespace Marmot
