#include "Marmot/DisplacementFiniteStrainFBarElement.h"
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

  const static bool C3D4FB_isRegistered = MarmotElementFactory::
    registerElement( "C3D4FB",
                     makeFactoryFunction< DisplacementFiniteStrainFBarElement< 3, 4 >,
                                          FullIntegration,
                                          DisplacementFiniteStrainFBarElement< 3, 4 >::SectionType::Solid >() );

  const static bool C3D8FB_isRegistered = MarmotElementFactory::
    registerElement( "C3D8FB",
                     makeFactoryFunction< DisplacementFiniteStrainFBarElement< 3, 8 >,
                                          FullIntegration,
                                          DisplacementFiniteStrainFBarElement< 3, 8 >::SectionType::Solid >() );

} // namespace Marmot::Elements::Registration
