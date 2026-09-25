#include "Marmot/GradientEnhancedFiniteStrainDruckerPragerDamage.h"
#include <algorithm>
#include <cmath>

namespace Marmot::Materials {

  using Damage = GradientEnhancedFiniteStrainDruckerPragerDamage;

  double Damage::xs( double Rs ) const
  {
    const double As = modelParameters.As;
    return Rs < 1.0 ? 1.0 + As * Rs * Rs : 1.0 + As * ( 4.0 * std::sqrt( Rs ) - 3.0 );
  }

  double Damage::deltaAlphaLocal( const Eigen::Vector3d& dEp ) const
  {
    const double dEpVol = dEp.sum();
    if ( dEpVol <= 0.0 )
      return 0.0; // only dilatant plastic flow drives the damage

    double dEpNeg = 0.0;
    for ( int a = 0; a < 3; a++ )
      dEpNeg += std::max( -dEp( a ), 0.0 );

    return dEpVol / xs( dEpNeg / dEpVol );
  }

  double Damage::omega( double kappa ) const
  {
    if ( kappa <= 0.0 )
      return 0.0;
    return std::min( 1.0 - std::exp( -kappa / modelParameters.softeningModulus ), modelParameters.maxDamage );
  }

  Damage::Result Damage::computeDamage( double alphaLocal, double alphaNonLocal, double kappaOld ) const
  {
    const double m             = modelParameters.weightingParameter;
    const double alphaWeighted = m * alphaNonLocal + ( 1.0 - m ) * alphaLocal;
    const double kappa         = std::max( kappaOld, alphaWeighted );
    return { alphaWeighted, kappa, omega( kappa ) };
  }

} // namespace Marmot::Materials
