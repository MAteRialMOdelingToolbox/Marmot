#include "Marmot/GeneralGradientEnhancedDisplacementFiniteElement.h"
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

  using namespace Marmot::Discretizations::Elements::Quadrature;

  const static bool GCPS4_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "GCPS4",
                     makeFactoryFunction< GeneralGradientEnhancedDisplacementFiniteElement< 2, 4 >,
                                          gauss( 2 ),
                                          GeneralGradientEnhancedDisplacementFiniteElement< 2, 4 >::PlaneStress >() );

  const static bool GCPS8_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "GCPS8",
                     makeFactoryFunction< GeneralGradientEnhancedDisplacementFiniteElement< 2, 8 >,
                                          gauss( 4 ),
                                          GeneralGradientEnhancedDisplacementFiniteElement< 2, 8 >::PlaneStress >() );

  const static bool GCPE4_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "GCPE4",
                     makeFactoryFunction< GeneralGradientEnhancedDisplacementFiniteElement< 2, 4 >,
                                          gauss( 2 ),
                                          GeneralGradientEnhancedDisplacementFiniteElement< 2, 4 >::PlaneStrain >() );

  const static bool GCPE8_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "GCPE8",
                     makeFactoryFunction< GeneralGradientEnhancedDisplacementFiniteElement< 2, 8 >,
                                          gauss( 4 ),
                                          GeneralGradientEnhancedDisplacementFiniteElement< 2, 8 >::PlaneStrain >() );
  const static bool GCPS8R_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "GCPS8R",
                     makeFactoryFunction< GeneralGradientEnhancedDisplacementFiniteElement< 2, 8 >,
                                          gauss( 3 ),
                                          GeneralGradientEnhancedDisplacementFiniteElement< 2, 8 >::PlaneStress >() );
  const static bool GCPE8R_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "GCPE8R",
                     makeFactoryFunction< GeneralGradientEnhancedDisplacementFiniteElement< 2, 8 >,
                                          gauss( 3 ),
                                          GeneralGradientEnhancedDisplacementFiniteElement< 2, 8 >::PlaneStrain >() );

  const static bool GC3D8_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "GC3D8",
                     makeFactoryFunction<
                       GeneralGradientEnhancedDisplacementFiniteElement< 3, 8 >,
                       gauss( 2 ),
                       GeneralGradientEnhancedDisplacementFiniteElement< 3, 8 >::SectionType::Solid >() );

  const static bool GC3D20_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "GC3D20",
                     makeFactoryFunction<
                       GeneralGradientEnhancedDisplacementFiniteElement< 3, 20 >,
                       gauss( 4 ),
                       GeneralGradientEnhancedDisplacementFiniteElement< 3, 20 >::SectionType::Solid >() );

  const static bool GC3D20R_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "GC3D20R",
                     makeFactoryFunction<
                       GeneralGradientEnhancedDisplacementFiniteElement< 3, 20 >,
                       gauss( 3 ),
                       GeneralGradientEnhancedDisplacementFiniteElement< 3, 20 >::SectionType::Solid >() );

  const static bool G2GCPS4_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "G2GCPS4",
                     makeFactoryFunction<
                       GeneralGradientEnhancedDisplacementFiniteElement< 2, 4, 2 >,
                       gauss( 2 ),
                       GeneralGradientEnhancedDisplacementFiniteElement< 2, 4, 2 >::PlaneStress >() );

  const static bool G2GCPS8_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "G2GCPS8",
                     makeFactoryFunction<
                       GeneralGradientEnhancedDisplacementFiniteElement< 2, 8, 2 >,
                       gauss( 4 ),
                       GeneralGradientEnhancedDisplacementFiniteElement< 2, 8, 2 >::PlaneStress >() );

  const static bool G2GCPS8R_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "G2GCPS8R",
                     makeFactoryFunction<
                       GeneralGradientEnhancedDisplacementFiniteElement< 2, 8, 2 >,
                       gauss( 3 ),
                       GeneralGradientEnhancedDisplacementFiniteElement< 2, 8, 2 >::PlaneStress >() );

  const static bool G2GC3D8_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "G2GC3D8",
                     makeFactoryFunction<
                       GeneralGradientEnhancedDisplacementFiniteElement< 3, 8, 2 >,
                       gauss( 2 ),
                       GeneralGradientEnhancedDisplacementFiniteElement< 3, 8, 2 >::SectionType::Solid >() );

  const static bool G2GC3D20_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "G2GC3D20",
                     makeFactoryFunction<
                       GeneralGradientEnhancedDisplacementFiniteElement< 3, 20, 2 >,
                       gauss( 4 ),
                       GeneralGradientEnhancedDisplacementFiniteElement< 3, 20, 2 >::SectionType::Solid >() );

  const static bool G2GC3D20R_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "G2GC3D20R",
                     makeFactoryFunction<
                       GeneralGradientEnhancedDisplacementFiniteElement< 3, 20, 2 >,
                       gauss( 3 ),
                       GeneralGradientEnhancedDisplacementFiniteElement< 3, 20, 2 >::SectionType::Solid >() );

  const static bool G6GCPS8_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "G6GCPS8",
                     makeFactoryFunction<
                       GeneralGradientEnhancedDisplacementFiniteElement< 2, 8, 6 >,
                       gauss( 4 ),
                       GeneralGradientEnhancedDisplacementFiniteElement< 2, 8, 6 >::PlaneStress >() );

  const static bool G6GC3D8_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "G6GC3D8",
                     makeFactoryFunction<
                       GeneralGradientEnhancedDisplacementFiniteElement< 3, 8, 6 >,
                       gauss( 2 ),
                       GeneralGradientEnhancedDisplacementFiniteElement< 3, 8, 6 >::SectionType::Solid >() );

  const static bool G6GC3D20_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "G6GC3D20",
                     makeFactoryFunction<
                       GeneralGradientEnhancedDisplacementFiniteElement< 3, 20, 6 >,
                       gauss( 4 ),
                       GeneralGradientEnhancedDisplacementFiniteElement< 3, 20, 6 >::SectionType::Solid >() );

  const static bool G6GC3D20R_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "G6GC3D20R",
                     makeFactoryFunction<
                       GeneralGradientEnhancedDisplacementFiniteElement< 3, 20, 6 >,
                       gauss( 3 ),
                       GeneralGradientEnhancedDisplacementFiniteElement< 3, 20, 6 >::SectionType::Solid >() );

  const static bool G6GCPS8M_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "G6GCPS8M",
                     makeFactoryFunction<
                       GeneralGradientEnhancedDisplacementFiniteElement< 2, 8, 6, 4 >,
                       gauss( 4 ),
                       GeneralGradientEnhancedDisplacementFiniteElement< 2, 8, 6, 4 >::PlaneStress >() );

  const static bool G6GC3D20RM_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "G6GC3D20RM",
                     makeFactoryFunction<
                       GeneralGradientEnhancedDisplacementFiniteElement< 3, 20, 6, 8 >,
                       gauss( 3 ),
                       GeneralGradientEnhancedDisplacementFiniteElement< 3, 20, 6, 8 >::SectionType::Solid >() );

  const static bool G6GC3D20M_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "G6GC3D20M",
                     makeFactoryFunction<
                       GeneralGradientEnhancedDisplacementFiniteElement< 3, 20, 6, 8 >,
                       gauss( 4 ),
                       GeneralGradientEnhancedDisplacementFiniteElement< 3, 20, 6, 8 >::SectionType::Solid >() );

} // namespace
