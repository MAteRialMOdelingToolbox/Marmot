#include "Marmot/DisplacementFiniteStrainULElement.h"
#include "Marmot/MarmotElementFactory.h"
#include "Marmot/MarmotFiniteElement.h"

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

  const static bool CX8RUL_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "CX8RUL",
                     makeFactoryFunction< AxiSymmetricDisplacementFiniteStrainULElement< 8 >,
                                          ReducedIntegration,
                                          AxiSymmetricDisplacementFiniteStrainULElement< 8 >::PlaneStrain >() );

  const static bool CX8UL_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "CX8UL",
                     makeFactoryFunction< AxiSymmetricDisplacementFiniteStrainULElement< 8 >,
                                          FullIntegration,
                                          AxiSymmetricDisplacementFiniteStrainULElement< 8 >::PlaneStrain >() );

  const static bool CPE8RUL_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "CPE8RUL",
                     makeFactoryFunction< DisplacementFiniteStrainULElement< 2, 8 >,
                                          ReducedIntegration,
                                          DisplacementFiniteStrainULElement< 2, 8 >::PlaneStrain >() );

  const static bool CPE8UL_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "CPE8UL",
                     makeFactoryFunction< DisplacementFiniteStrainULElement< 2, 8 >,
                                          FullIntegration,
                                          DisplacementFiniteStrainULElement< 2, 8 >::PlaneStrain >() );

  const static bool C3D8UL_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "C3D8UL",
                     makeFactoryFunction< DisplacementFiniteStrainULElement< 3, 8 >,
                                          FullIntegration,
                                          DisplacementFiniteStrainULElement< 3, 8 >::SectionType::Solid >() );

  const static bool C3D20RUL_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "C3D20RUL",
                     makeFactoryFunction< DisplacementFiniteStrainULElement< 3, 20 >,
                                          ReducedIntegration,
                                          DisplacementFiniteStrainULElement< 3, 20 >::SectionType::Solid >() );

  const static bool C3D20UL_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "C3D20UL",
                     makeFactoryFunction< DisplacementFiniteStrainULElement< 3, 20 >,
                                          FullIntegration,
                                          DisplacementFiniteStrainULElement< 3, 20 >::SectionType::Solid >() );

} // namespace
