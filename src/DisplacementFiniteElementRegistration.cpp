#include "Marmot/DisplacementFiniteElement.h"
#include "Marmot/Marmot.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotFiniteElementSpatialWrapper.h"

namespace Marmot::Elements::Registration {

  enum DisplacementElementCode {

    /* TAG explanation
     * XXXX
     * ||||___    4: type of element
     * |||____    3: active fields
     * ||_____    2: number of nodes
     * |______    1: number of nodes
     *
     * active fields:   0: displacement,
     *
     * type of element: 1: 1D full integration,
     *                  2: 2D full integration, plane stress
     *                  3: 3D full integration,
     *                  4: 1D red. integration,
     *                  5: 2D red. integration, plane stress
     *                  6: 3D red. integration
     *                  7: 2D full integration, plane strain
     *                  8: 2D red. integration, plane strain
     * */

    // Truss 2D
    T2D2 = 202,
    // Plane stress 2D
    CPS4  = 402,
    CPS8R = 805,

    // Plane Strain 2D
    CPE4  = 407,
    CPE8R = 808,

    // Solid
    C3D8   = 803,
    C3D8R  = 806,
    C3D20  = 2003,
    C3D20R = 2006
  };

  template < class T,
             Marmot::FiniteElement::Quadrature::IntegrationTypes integrationType,
             typename T::SectionType                             sectionType >
  MarmotLibrary::MarmotElementFactory::elementFactoryFunction makeFactoryFunction()
  {
    return []( int elementID ) -> MarmotElement* { return new T( elementID, integrationType, sectionType ); };
  }

  using namespace MarmotLibrary;
  using namespace Marmot::FiniteElement::Quadrature;

  const static bool CPS4_isRegistered = MarmotElementFactory::
    registerElement( "CPS4",
                     DisplacementElementCode::CPS4,
                     makeFactoryFunction< DisplacementFiniteElement< 2, 4 >,
                                          FullIntegration,
                                          DisplacementFiniteElement< 2, 4 >::PlaneStress >() );

  const static bool CPE4_isRegistered = MarmotElementFactory::
    registerElement( "CPE4",
                     DisplacementElementCode::CPE4,
                     makeFactoryFunction< DisplacementFiniteElement< 2, 4 >,
                                          FullIntegration,
                                          DisplacementFiniteElement< 2, 4 >::PlaneStrain >() );

  const static bool CPS8R_isRegistered = MarmotElementFactory::
    registerElement( "CPS8R",
                     DisplacementElementCode::CPS8R,
                     makeFactoryFunction< DisplacementFiniteElement< 2, 8 >,
                                          ReducedIntegration,
                                          DisplacementFiniteElement< 2, 8 >::PlaneStress >() );

  const static bool CPE8R_isRegistered = MarmotElementFactory::
    registerElement( "CPE8R",
                     DisplacementElementCode::CPE8R,
                     makeFactoryFunction< DisplacementFiniteElement< 2, 8 >,
                                          ReducedIntegration,
                                          DisplacementFiniteElement< 2, 8 >::PlaneStrain >() );

  const static bool C3D8_isRegistered = MarmotLibrary::MarmotElementFactory::
    registerElement( "C3D8", DisplacementElementCode::C3D8, []( int elementID ) -> MarmotElement* {
      return new DisplacementFiniteElement< 3,
                                            8 >( elementID,
                                                 Marmot::FiniteElement::Quadrature::IntegrationTypes::FullIntegration,
                                                 DisplacementFiniteElement< 3, 8 >::SectionType::Solid );
    } );

  const static bool C3D20_isRegistered = MarmotLibrary::MarmotElementFactory::
    registerElement( "C3D20", DisplacementElementCode::C3D20, []( int elementID ) -> MarmotElement* {
      return new DisplacementFiniteElement< 3,
                                            20 >( elementID,
                                                  Marmot::FiniteElement::Quadrature::IntegrationTypes::FullIntegration,
                                                  DisplacementFiniteElement< 3, 20 >::SectionType::Solid );
    } );

  const static bool C3D20R_isRegistered = MarmotLibrary::MarmotElementFactory::
    registerElement( "C3D20R", DisplacementElementCode::C3D20R, []( int elementID ) -> MarmotElement* {
      return new DisplacementFiniteElement< 3, 20 >( elementID,
                                                     Marmot::FiniteElement::Quadrature::IntegrationTypes::
                                                       ReducedIntegration,
                                                     DisplacementFiniteElement< 3, 20 >::SectionType::Solid );
    } );

  MarmotElement* generateT2D2( int elementID )
  {
    auto uelT2D2 = std::unique_ptr< MarmotElement >(
      new DisplacementFiniteElement< 1, 2 >( elementID,
                                             Marmot::FiniteElement::Quadrature::IntegrationTypes::FullIntegration,
                                             DisplacementFiniteElement< 1, 2 >::SectionType::UniaxialStress ) );
    constexpr static int indicesToBeWrapped[] = { 0, 1 };
    constexpr static int nIndicesToBeWrapped  = 2;
    return new MarmotElementSpatialWrapper( 2, 1, 2, 2, indicesToBeWrapped, nIndicesToBeWrapped, std::move( uelT2D2 ) );
  }
  const static bool
    T2D2_isRegistered = MarmotLibrary::MarmotElementFactory::registerElement( "T2D2",
                                                                              DisplacementElementCode::T2D2,
                                                                              generateT2D2 );
} // namespace Marmot::Elements::Registration
