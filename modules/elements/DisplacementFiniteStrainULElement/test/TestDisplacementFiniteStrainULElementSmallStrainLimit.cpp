#include "Marmot/MarmotElementTesting.h"
#include "Marmot/MarmotTesting.h"
#include <Eigen/Dense>
#include <string>
#include <vector>

using namespace Marmot::Testing;
using Marmot::Discretizations::Elements::Element;

// ---------------------------------------------------------------------------------------------
// In the small strain limit, at rest, the axisymmetric finite strain elements CX8UL and CX8RUL load
// a face as the small strain elements CAX8 and CAX8R do: over the ring 2 pi r of the face, not over
// the face times a thickness. The thickness 7 given to the finite strain elements must not matter.
// ---------------------------------------------------------------------------------------------

namespace {

  /// A ring section r in [1, 2], z in [0, 1], distorted, nodes of an 8-node quadrilateral.
  const std::vector< double > ringSectionCoordinates =
    { 1.0, 0.1, 2.2, -0.2, 2.0, 1.1, 0.9, 0.8, 1.6, -0.1, 2.15, 0.45, 1.45, 1.0, 0.95, 0.4 };

  /**
   * @brief The nodal loads of a distributed load on a face of an element at rest.
   * @param element The element.
   * @param loadType Pressure or surface traction.
   * @param face The face, 1-based.
   * @param load The pressure, or the traction vector.
   * @return The nodal loads.
   */
  Eigen::VectorXd distributedLoadAtRestOf( Element&                      element,
                                           Element::DistributedLoadTypes loadType,
                                           int                           face,
                                           const std::vector< double >&  load )
  {
    const int             nDof   = element.getNDofPerElement();
    Eigen::VectorXd       P      = Eigen::VectorXd::Zero( nDof );
    Eigen::MatrixXd       K      = Eigen::MatrixXd::Zero( nDof, nDof );
    const Eigen::VectorXd QTotal = Eigen::VectorXd::Zero( nDof );
    element.computeDistributedLoad( loadType, P.data(), K.data(), face, load.data(), QTotal.data(), 0.0, 1.0 );
    return P;
  }

} // namespace

void testAxiSymmetricSurfaceLoadsAreThoseOfTheSmallStrainElement()
{
  for ( const auto& [finiteStrainType, smallStrainType] :
        { std::make_pair( "CX8UL", "CAX8" ), std::make_pair( "CX8RUL", "CAX8R" ) } ) {
    ElementUnderTest
                     finiteStrain( finiteStrainType,
                    ringSectionCoordinates,
                                   { 7.0 },
                    "COMPRESSIBLEHYPERELASTICITY",
                    // NeoHooke, DollSchweizerhof: G, K, alpha, beta, density
                                   { 0.0, 6.0, 1000.0 / ( 2 * ( 1 + 0.3 ) ), 1000.0 / ( 3 * ( 1 - 2 * 0.3 ) ), 1.0, 3.0, 1.0 } );
    ElementUnderTest smallStrain( smallStrainType, ringSectionCoordinates, {}, "LINEARELASTIC", { 1000.0, 0.3, 1.0 } );

    for ( int face = 1; face <= 4; face++ ) {
      const std::string at = std::string( finiteStrainType ) + " and " + smallStrainType + ", face " +
                             std::to_string( face );

      throwExceptionOnFailure( checkIfEqual( Eigen::MatrixXd( distributedLoadAtRestOf( *finiteStrain.element,
                                                                                       Element::Pressure,
                                                                                       face,
                                                                                       { 3.0 } ) ),
                                             Eigen::MatrixXd( distributedLoadAtRestOf( *smallStrain.element,
                                                                                       Element::Pressure,
                                                                                       face,
                                                                                       { 3.0 } ) ),
                                             1e-12 ),
                               at + ": the nodal loads of a pressure differ." );

      const std::vector< double > traction = { 0.3, -1.7 };
      throwExceptionOnFailure( checkIfEqual( Eigen::MatrixXd( distributedLoadAtRestOf( *finiteStrain.element,
                                                                                       Element::SurfaceTraction,
                                                                                       face,
                                                                                       traction ) ),
                                             Eigen::MatrixXd( distributedLoadAtRestOf( *smallStrain.element,
                                                                                       Element::SurfaceTraction,
                                                                                       face,
                                                                                       traction ) ),
                                             1e-12 ),
                               at + ": the nodal loads of a surface traction differ." );
    }
  }
}

int main()
{
  auto tests = std::vector< std::function< void() > >{
    testAxiSymmetricSurfaceLoadsAreThoseOfTheSmallStrainElement,
  };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
