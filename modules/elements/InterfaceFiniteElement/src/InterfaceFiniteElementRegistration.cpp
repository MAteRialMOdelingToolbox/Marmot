#include "Marmot/InterfaceFiniteElement.h"
#include "Marmot/MarmotElementFactory.h"
#include "Marmot/MarmotFiniteElement.h"

namespace {
  using namespace Marmot::Discretizations::Elements;

  /**
   * @brief Create the factory function of an element type.
   *
   * @tparam T element type
   * @tparam quadratureRule quadrature rule of the created elements
   * @tparam sectionType section type of the created elements
   * @return a function that creates an element of type T from an element number
   */
  template < class T,
             Marmot::Discretizations::Elements::Quadrature::QuadratureRule quadratureRule,
             typename T::SectionType                                       sectionType >
  Marmot::Registry::ElementFactory::elementFactoryFunction makeFactoryFunction()
  {
    return []( int elementID ) -> Marmot::Discretizations::Elements::Element* {
      return new T( elementID, quadratureRule, sectionType );
    };
  }

  const static bool ILINE2_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "ILINE2",
                     makeFactoryFunction< InterfaceFiniteElement< 2, 4 >,
                                          Quadrature::gaussExactToDegree( 2 ),
                                          InterfaceFiniteElement< 2, 4 >::SectionType::Interface >() );

  const static bool IQUAD4_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "IQUAD4",
                     makeFactoryFunction< InterfaceFiniteElement< 3, 8 >,
                                          Quadrature::gaussExactToDegree( 2 ),
                                          InterfaceFiniteElement< 3, 8 >::SectionType::Interface >() );

} // namespace