#include "Marmot/RankineSoftening.h"
#include "Marmot/MarmotConstants.h"
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotTensor.h"

namespace Marmot::Materials::RankineSoftening {

  using namespace Eigen;

  // index pairs (A,B) of the off-diagonal basis tensors
  static constexpr int pairs[3][2] = { { 0, 1 }, { 0, 2 }, { 1, 2 } };

  void searchActiveSurfaces( YieldSurfFlagArr                                        activeSurfaces,
                             const std::function< bool( const YieldSurfFlagArr& ) >& returnMappingAttempt,
                             const std::string&                                      caller )
  {
    YieldSurfaceCombinationManager yieldSurfCombiManager;
    while ( true ) {
      yieldSurfCombiManager.markYieldFlagCombinationAsUsed( activeSurfaces );

      if ( returnMappingAttempt( activeSurfaces ) )
        return;

      if ( !yieldSurfCombiManager.getAnotherYieldFlagCombination( activeSurfaces ) )
        throw Utilities::StressUpdateFailed(
          Utilities::MakeString() << caller << ": return mapping failed for all combinations of yield surfaces" );
    }
  }

  SpectralBasis spectralBasis( const Math::Matrix3d& V )
  {
    using namespace Math::Tensors;
    SpectralBasis basis;
    for ( int a = 0; a < 3; a++ ) {
      const Math::Matrix3d Ea  = V.col( a ) * V.col( a ).transpose();
      basis.E.col( a )         = flatten( Ea );
      const auto [A, B]        = pairs[a];
      const Math::Matrix3d NAB = Math::Constants::sqrt2 / 2. *
                                 ( V.col( A ) * V.col( B ).transpose() + V.col( B ) * V.col( A ).transpose() );
      basis.N.col( a ) = flatten( NAB );
    }
    return basis;
  }

  Math::Vector3d spinCoefficients( const Math::Vector3d& sigTrial,
                                   const Math::Vector3d& sig,
                                   const Math::Matrix3d& D,
                                   double                tol )
  {
    Math::Vector3d c;
    for ( int i = 0; i < 3; i++ ) {
      const auto [a, b]      = pairs[i];
      const double dSigTrial = sigTrial( a ) - sigTrial( b );
      c( i )                 = std::abs( dSigTrial ) > tol ? ( sig( a ) - sig( b ) ) / dSigTrial
                                                           : 0.5 * ( D( a, a ) - D( a, b ) + D( b, b ) - D( b, a ) );
    }
    return c;
  }

} // namespace Marmot::Materials::RankineSoftening
