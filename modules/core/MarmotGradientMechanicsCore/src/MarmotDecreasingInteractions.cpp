#include "Marmot/MarmotDecreasingInteractions.h"

namespace Marmot::GradientDamage::DecreasingInteractions {

  namespace FirstOrderDerived {

    std::tuple< double, double > poh( double omega, double eta, double R )
    {

      const double g         = Marmot::GradientDamage::DecreasingInteractions::poh( omega, eta, R );
      const double dg_domega = -( 1. - R ) * eta * exp( -eta * omega ) / ( 1. - exp( -eta ) );

      if ( omega > 0 ) {
        return { g, dg_domega };
      }
      else {
        return { 1, 0 };
      }
    }
  } // namespace FirstOrderDerived

  namespace SecondOrderDerived {

    std::tuple< double, double, double > poh( double omega, double eta, double R )
    {
      const double g           = Marmot::GradientDamage::DecreasingInteractions::poh( omega, eta, R );
      const double dg_domega   = -( 1. - R ) * eta * exp( -eta * omega ) / ( 1. - exp( -eta ) );
      const double d2g_domega2 = ( 1. - R ) * eta * eta * exp( -eta * omega ) / ( 1. - exp( -eta ) );
      return { g, dg_domega, d2g_domega2 };
    }

  } // namespace SecondOrderDerived

} // namespace Marmot::GradientDamage::DecreasingInteractions
