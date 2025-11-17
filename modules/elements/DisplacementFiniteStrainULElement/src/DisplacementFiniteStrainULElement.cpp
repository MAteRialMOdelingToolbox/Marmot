#include "Marmot/DisplacementFiniteStrainULElement.h"
#include "Marmot/MarmotElementFactory.h"
#include "Marmot/MarmotFiniteElement.h"

namespace Marmot::Elements::Registration {

  template < class T,
             Marmot::FiniteElement::Quadrature::IntegrationTypes integrationType,
             typename T::SectionType                             sectionType >
  MarmotLibrary::MarmotElementFactory::elementFactoryFunction makeFactoryFunction()
  {
    return []( int elementID ) -> MarmotElement* { return new T( elementID, integrationType, sectionType ); };
  }

  using namespace MarmotLibrary;
  using namespace Marmot::FiniteElement::Quadrature;

  const static bool CX8RUL_isRegistered = MarmotElementFactory::
    registerElement( "CX8RUL",
                     makeFactoryFunction< AxiSymmetricDisplacementFiniteStrainULElement< 8 >,
                                          ReducedIntegration,
                                          AxiSymmetricDisplacementFiniteStrainULElement< 8 >::PlaneStrain >() );

  const static bool CX8UL_isRegistered = MarmotElementFactory::
    registerElement( "CX8UL",
                     makeFactoryFunction< AxiSymmetricDisplacementFiniteStrainULElement< 8 >,
                                          FullIntegration,
                                          AxiSymmetricDisplacementFiniteStrainULElement< 8 >::PlaneStrain >() );

  const static bool CPE8RGradientEnhancedMicropolar_isRegistered = MarmotElementFactory::
    registerElement( "CPE8RUL",
                     makeFactoryFunction< DisplacementFiniteStrainULElement< 2, 8 >,
                                          ReducedIntegration,
                                          DisplacementFiniteStrainULElement< 2, 8 >::PlaneStrain >() );

  const static bool C3D8UL_isRegistered = MarmotLibrary::MarmotElementFactory::
    registerElement( "C3D8UL",
                     makeFactoryFunction< DisplacementFiniteStrainULElement< 3, 8 >,
                                          FullIntegration,
                                          DisplacementFiniteStrainULElement< 3, 8 >::SectionType::Solid >() );

  const static bool C3D20RUL_isRegistered = MarmotLibrary::MarmotElementFactory::
    registerElement( "C3D20RUL",
                     makeFactoryFunction< DisplacementFiniteStrainULElement< 3, 20 >,
                                          ReducedIntegration,
                                          DisplacementFiniteStrainULElement< 3, 20 >::SectionType::Solid >() );

  const static bool C3D20UL_isRegistered = MarmotLibrary::MarmotElementFactory::
    registerElement( "C3D20UL",
                     makeFactoryFunction< DisplacementFiniteStrainULElement< 3, 20 >,
                                          FullIntegration,
                                          DisplacementFiniteStrainULElement< 3, 20 >::SectionType::Solid >() );

} // namespace Marmot::Elements::Registration
