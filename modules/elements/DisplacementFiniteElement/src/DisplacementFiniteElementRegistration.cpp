#include "Marmot/DisplacementFiniteElement.h"
#include "Marmot/MarmotElementFactory.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotFiniteElementSpatialWrapper.h"

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

  using namespace Marmot::Discretizations::Elements::Quadrature;

  // linear triangles: the strain is constant in the element, so a single Gauss point integrates them exactly
  const static bool CPS3_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "CPS3",
                     makeFactoryFunction< DisplacementFiniteElement< 2, 3 >,
                                          gauss( 1 ),
                                          DisplacementFiniteElement< 2, 3 >::PlaneStress >() );

  const static bool CPE3_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "CPE3",
                     makeFactoryFunction< DisplacementFiniteElement< 2, 3 >,
                                          gauss( 1 ),
                                          DisplacementFiniteElement< 2, 3 >::PlaneStrain >() );

  const static bool CPS4_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "CPS4",
                     makeFactoryFunction< DisplacementFiniteElement< 2, 4 >,
                                          gauss( 2 ),
                                          DisplacementFiniteElement< 2, 4 >::PlaneStress >() );

  const static bool CPE4_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "CPE4",
                     makeFactoryFunction< DisplacementFiniteElement< 2, 4 >,
                                          gauss( 2 ),
                                          DisplacementFiniteElement< 2, 4 >::PlaneStrain >() );

  const static bool CPS8R_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "CPS8R",
                     makeFactoryFunction< DisplacementFiniteElement< 2, 8 >,
                                          gauss( 3 ),
                                          DisplacementFiniteElement< 2, 8 >::PlaneStress >() );

  const static bool CPE8R_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "CPE8R",
                     makeFactoryFunction< DisplacementFiniteElement< 2, 8 >,
                                          gauss( 3 ),
                                          DisplacementFiniteElement< 2, 8 >::PlaneStrain >() );

  const static bool CPE8_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "CPE8",
                     makeFactoryFunction< DisplacementFiniteElement< 2, 8 >,
                                          gauss( 4 ),
                                          DisplacementFiniteElement< 2, 8 >::PlaneStrain >() );

  const static bool C3D4_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "C3D4",
                     makeFactoryFunction< DisplacementFiniteElement< 3, 4 >,
                                          gauss( 1 ),
                                          DisplacementFiniteElement< 3, 4 >::SectionType::Solid >() );

  const static bool C3D8_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "C3D8",
                     makeFactoryFunction< DisplacementFiniteElement< 3, 8 >,
                                          gauss( 2 ),
                                          DisplacementFiniteElement< 3, 8 >::SectionType::Solid >() );

  const static bool C3D20_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "C3D20",
                     makeFactoryFunction< DisplacementFiniteElement< 3, 20 >,
                                          gauss( 4 ),
                                          DisplacementFiniteElement< 3, 20 >::SectionType::Solid >() );

  const static bool C3D20R_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "C3D20R",
                     makeFactoryFunction< DisplacementFiniteElement< 3, 20 >,
                                          gauss( 3 ),
                                          DisplacementFiniteElement< 3, 20 >::SectionType::Solid >() );

  Marmot::Discretizations::Elements::Element* generateT2D2( int elementID )
  {
    auto uelT2D2 = std::unique_ptr< Marmot::Discretizations::Elements::Element >(
      new DisplacementFiniteElement< 1, 2 >( elementID,
                                             Marmot::Discretizations::Elements::Quadrature::gauss( 2 ),
                                             DisplacementFiniteElement< 1, 2 >::SectionType::UniaxialStress ) );
    constexpr static int indicesToBeWrapped[] = { 0, 1 };
    constexpr static int nIndicesToBeWrapped  = 2;
    return new Marmot::Discretizations::Elements::ElementSpatialWrapper( 2,
                                                                         1,
                                                                         2,
                                                                         2,
                                                                         indicesToBeWrapped,
                                                                         nIndicesToBeWrapped,
                                                                         std::move( uelT2D2 ) );
  }
  const static bool T2D2_isRegistered = Marmot::Registry::ElementFactory::registerElement( "T2D2", generateT2D2 );
} // namespace
