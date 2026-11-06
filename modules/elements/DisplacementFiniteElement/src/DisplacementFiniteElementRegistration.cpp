#include "Marmot/DisplacementFiniteElement.h"
#include "Marmot/MarmotElementFactory.h"
#include "Marmot/MarmotFiniteElement.h"

namespace {
  using namespace Marmot::Discretizations::Elements;

  template < class T,
             Marmot::Discretizations::Elements::Quadrature::QuadratureRule quadratureRule,
             typename T::SectionType                                       sectionType >
  Marmot::Registry::ElementFactory::elementFactoryFunction makeFactoryFunction()
  {
    return []( int elementID ) -> Marmot::Discretizations::Elements::Element* {
      return new T( elementID, quadratureRule, sectionType );
    };
  }

  // linear triangles: the strain is constant in the element, so a single Gauss point integrates them exactly
  const static bool CPS3_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "CPS3",
                     makeFactoryFunction< DisplacementFiniteElement< 2, 3 >,
                                          Quadrature::gaussExactToDegree( 1 ),
                                          DisplacementFiniteElement< 2, 3 >::PlaneStress >() );

  const static bool CPE3_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "CPE3",
                     makeFactoryFunction< DisplacementFiniteElement< 2, 3 >,
                                          Quadrature::gaussExactToDegree( 1 ),
                                          DisplacementFiniteElement< 2, 3 >::PlaneStrain >() );

  const static bool CPS4_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "CPS4",
                     makeFactoryFunction< DisplacementFiniteElement< 2, 4 >,
                                          Quadrature::gaussExactToDegree( 2 ),
                                          DisplacementFiniteElement< 2, 4 >::PlaneStress >() );

  const static bool CPE4_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "CPE4",
                     makeFactoryFunction< DisplacementFiniteElement< 2, 4 >,
                                          Quadrature::gaussExactToDegree( 2 ),
                                          DisplacementFiniteElement< 2, 4 >::PlaneStrain >() );

  const static bool CPS8R_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "CPS8R",
                     makeFactoryFunction< DisplacementFiniteElement< 2, 8 >,
                                          Quadrature::gaussExactToDegree( 3 ),
                                          DisplacementFiniteElement< 2, 8 >::PlaneStress >() );

  const static bool CPE8R_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "CPE8R",
                     makeFactoryFunction< DisplacementFiniteElement< 2, 8 >,
                                          Quadrature::gaussExactToDegree( 3 ),
                                          DisplacementFiniteElement< 2, 8 >::PlaneStrain >() );

  const static bool CPE8_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "CPE8",
                     makeFactoryFunction< DisplacementFiniteElement< 2, 8 >,
                                          Quadrature::gaussExactToDegree( 4 ),
                                          DisplacementFiniteElement< 2, 8 >::PlaneStrain >() );

  const static bool C3D4_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "C3D4",
                     makeFactoryFunction< DisplacementFiniteElement< 3, 4 >,
                                          Quadrature::gaussExactToDegree( 1 ),
                                          DisplacementFiniteElement< 3, 4 >::SectionType::Solid >() );

  const static bool C3D8_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "C3D8",
                     makeFactoryFunction< DisplacementFiniteElement< 3, 8 >,
                                          Quadrature::gaussExactToDegree( 2 ),
                                          DisplacementFiniteElement< 3, 8 >::SectionType::Solid >() );

  const static bool C3D20_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "C3D20",
                     makeFactoryFunction< DisplacementFiniteElement< 3, 20 >,
                                          Quadrature::gaussExactToDegree( 4 ),
                                          DisplacementFiniteElement< 3, 20 >::SectionType::Solid >() );

  const static bool C3D20R_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "C3D20R",
                     makeFactoryFunction< DisplacementFiniteElement< 3, 20 >,
                                          Quadrature::gaussExactToDegree( 3 ),
                                          DisplacementFiniteElement< 3, 20 >::SectionType::Solid >() );

} // namespace
