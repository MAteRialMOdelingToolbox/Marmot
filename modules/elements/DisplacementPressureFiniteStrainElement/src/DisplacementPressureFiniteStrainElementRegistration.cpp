#include "Marmot/DisplacementPressureFiniteStrainElement.h"
#include "Marmot/MarmotElementFactory.h"
#include "Marmot/MarmotFiniteElement.h"

namespace Marmot::Elements::Registration {

  template < class T, Marmot::FiniteElement::Quadrature::IntegrationTypes integrationType >
  MarmotLibrary::MarmotElementFactory::elementFactoryFunction makeFactoryFunction()
  {
    return []( int elementID ) -> MarmotElement* { return new T( elementID, integrationType ); };
  }

  using namespace MarmotLibrary;
  using namespace Marmot::FiniteElement::Quadrature;

  const static bool C3D10MP_isRegistered = MarmotElementFactory::
    registerElement( "C3D10MP",
                     makeFactoryFunction< DisplacementPressureFiniteStrainElement< 10, 4 >, FullIntegration >() );

  const static bool C3D20MP_isRegistered = MarmotElementFactory::
    registerElement( "C3D20MP",
                     makeFactoryFunction< DisplacementPressureFiniteStrainElement< 20, 8 >, FullIntegration >() );

} // namespace Marmot::Elements::Registration
