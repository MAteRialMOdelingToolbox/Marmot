#include "Marmot/DisplacementFiniteElement.h"
#include "Marmot/MarmotElementFactory.h"
#include "Marmot/MarmotFiniteElement.h"

namespace {
  using namespace Marmot::Discretizations::Elements;

  template < class T,
             Marmot::Discretizations::Elements::Quadrature::IntegrationTypes integrationType,
             typename T::SectionType                                         sectionType >
  Marmot::Registry::ElementFactory::elementFactoryFunction makeFactoryFunction()
  {
    return []( int elementID ) -> Marmot::Discretizations::Elements::Element* {
      return new T( elementID, integrationType, sectionType );
    };
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

} // namespace
