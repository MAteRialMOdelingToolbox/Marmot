#include "Marmot/DisplacementFiniteElement.h"
#include "Marmot/MarmotElementFactory.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotFiniteElementSpatialWrapper.h"

namespace {
  using namespace Marmot::Discretizations::Elements;

  template < class T,
             Marmot::Discretizations::Elements::Quadrature::IntegrationTypes integrationType,
             typename T::SectionType                             sectionType >
  Marmot::Registry::ElementFactory::elementFactoryFunction makeFactoryFunction()
  {
    return []( int elementID ) -> Marmot::Discretizations::Elements::Element* { return new T( elementID, integrationType, sectionType ); };
  }

  using namespace Marmot::Discretizations::Elements::Quadrature;

  const static bool CPS4_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "CPS4",
                     makeFactoryFunction< DisplacementFiniteElement< 2, 4 >,
                                          FullIntegration,
                                          DisplacementFiniteElement< 2, 4 >::PlaneStress >() );

  const static bool CPE4_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "CPE4",
                     makeFactoryFunction< DisplacementFiniteElement< 2, 4 >,
                                          FullIntegration,
                                          DisplacementFiniteElement< 2, 4 >::PlaneStrain >() );

  const static bool CPS8R_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "CPS8R",
                     makeFactoryFunction< DisplacementFiniteElement< 2, 8 >,
                                          ReducedIntegration,
                                          DisplacementFiniteElement< 2, 8 >::PlaneStress >() );

  const static bool CPE8R_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "CPE8R",
                     makeFactoryFunction< DisplacementFiniteElement< 2, 8 >,
                                          ReducedIntegration,
                                          DisplacementFiniteElement< 2, 8 >::PlaneStrain >() );

  const static bool CPE8_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "CPE8",
                     makeFactoryFunction< DisplacementFiniteElement< 2, 8 >,
                                          FullIntegration,
                                          DisplacementFiniteElement< 2, 8 >::PlaneStrain >() );

  const static bool C3D8_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "C3D8",
                     makeFactoryFunction< DisplacementFiniteElement< 3, 8 >,
                                          FullIntegration,
                                          DisplacementFiniteElement< 3, 8 >::SectionType::Solid >() );

  const static bool C3D20_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "C3D20",
                     makeFactoryFunction< DisplacementFiniteElement< 3, 20 >,
                                          FullIntegration,
                                          DisplacementFiniteElement< 3, 20 >::SectionType::Solid >() );

  const static bool C3D20R_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "C3D20R",
                     makeFactoryFunction< DisplacementFiniteElement< 3, 20 >,
                                          ReducedIntegration,
                                          DisplacementFiniteElement< 3, 20 >::SectionType::Solid >() );

  Marmot::Discretizations::Elements::Element* generateT2D2( int elementID )
  {
    auto uelT2D2 = std::unique_ptr< Marmot::Discretizations::Elements::Element >(
      new DisplacementFiniteElement< 1, 2 >( elementID,
                                             Marmot::Discretizations::Elements::Quadrature::IntegrationTypes::FullIntegration,
                                             DisplacementFiniteElement< 1, 2 >::SectionType::UniaxialStress ) );
    constexpr static int indicesToBeWrapped[] = { 0, 1 };
    constexpr static int nIndicesToBeWrapped  = 2;
    return new Marmot::Discretizations::Elements::ElementSpatialWrapper( 2, 1, 2, 2, indicesToBeWrapped, nIndicesToBeWrapped, std::move( uelT2D2 ) );
  }
  const static bool T2D2_isRegistered = Marmot::Registry::ElementFactory::registerElement( "T2D2", generateT2D2 );
} // namespace
