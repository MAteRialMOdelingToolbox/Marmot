#include "Marmot/InterfaceFiniteElement.h"
#include "Marmot/MarmotElementFactory.h"
#include "Marmot/MarmotFiniteElement.h"

namespace Marmot::Elements::Registration {

  /**
   * @brief Create the factory function of an element type.
   *
   * @tparam T element type
   * @tparam integrationType integration type of the created elements
   * @tparam sectionType section type of the created elements
   * @return a function that creates an element of type T from an element number
   */
  template < class T,
             Marmot::FiniteElement::Quadrature::IntegrationTypes integrationType,
             typename T::SectionType                             sectionType >
  MarmotLibrary::MarmotElementFactory::elementFactoryFunction makeFactoryFunction()
  {
    return []( int elementID ) -> MarmotElement* { return new T( elementID, integrationType, sectionType ); };
  }

  using namespace MarmotLibrary;
  using namespace Marmot::FiniteElement::Quadrature;

  const static bool ILINE2_isRegistered = MarmotElementFactory::
    registerElement( "ILINE2",
                     makeFactoryFunction< InterfaceFiniteElement< 2, 4 >,
                                          FullIntegration,
                                          InterfaceFiniteElement< 2, 4 >::SectionType::Interface >() );

  const static bool IQUAD4_isRegistered = MarmotElementFactory::
    registerElement( "IQUAD4",
                     makeFactoryFunction< InterfaceFiniteElement< 3, 8 >,
                                          FullIntegration,
                                          InterfaceFiniteElement< 3, 8 >::SectionType::Interface >() );

} // namespace Marmot::Elements::Registration