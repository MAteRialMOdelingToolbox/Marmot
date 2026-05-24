#include "Marmot/MarmotPhaseFieldEnergyDegradation.h"
#include <autodiff/forward/dual.hpp>
#include <cmath>

namespace Marmot::PhaseField {

  namespace EnergyDegradationFunctions {

    namespace SecondOrderDerived {

      std::tuple< double, double, double > quadratic( const double pf )
      {

        const double s      = 1. - pf;
        const double ds_dpf = -1.;
        return { s * s, ds_dpf * 2. * s, 2. };
      }

      std::tuple< double, double, double > cubic( const double pf )
      {

        const double s      = 1. - pf;
        const double ds_dpf = -1.;

        return { 3 * s * s - 2 * s * s * s, ds_dpf * 6 * s - ds_dpf * 6 * s * s, 6 - 12 * s };
      }
      std::tuple< double, double, double > quartic( const double pf )
      {

        const double s      = 1. - pf;
        const double ds_dpf = -1.;

        return { 4 * s * s * s - 3 * s * s * s * s,
                 ds_dpf * 12 * s * s - ds_dpf * 12 * s * s * s,
                 24 * s - 36 * s * s };
      }
      std::tuple< double, double, double > generic( const double pf,
                                                    const double p,
                                                    const double a1,
                                                    const double a2,
                                                    const double a3 )
      {

        using namespace autodiff;
        autodiff::dual2nd pf_dual( pf );
        autodiff::seed< 1 >( pf_dual, 1 );
        autodiff::seed< 2 >( pf_dual, 1 );

        const dual2nd result = Marmot::PhaseField::EnergyDegradationFunctions::generic( pf_dual, p, a1, a2, a3 );

        return { result.val.val, derivative< 1 >( result ), derivative< 2 >( result ) };
      }
    }; // namespace SecondOrderDerived
  }    // namespace EnergyDegradationFunctions
};     // namespace Marmot::PhaseField
