#include "Marmot/GradientEnhancedFiniteStrainDisplacementElement.h"
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

  // 2D plane strain, 8-node
  const static bool GCPE8UL_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "GCPE8UL",
                     makeFactoryFunction< GradientEnhancedFiniteStrainDisplacementElement< 2, 8 >,
                                          FullIntegration,
                                          GradientEnhancedFiniteStrainDisplacementElement< 2, 8 >::PlaneStrain >() );

  const static bool GCPE8RUL_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "GCPE8RUL",
                     makeFactoryFunction< GradientEnhancedFiniteStrainDisplacementElement< 2, 8 >,
                                          ReducedIntegration,
                                          GradientEnhancedFiniteStrainDisplacementElement< 2, 8 >::PlaneStrain >() );

  // 3D 8-node
  const static bool GC3D8UL_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "GC3D8UL",
                     makeFactoryFunction<
                       GradientEnhancedFiniteStrainDisplacementElement< 3, 8 >,
                       FullIntegration,
                       GradientEnhancedFiniteStrainDisplacementElement< 3, 8 >::SectionType::Solid >() );

  // 3D 20-node
  const static bool GC3D20UL_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "GC3D20UL",
                     makeFactoryFunction<
                       GradientEnhancedFiniteStrainDisplacementElement< 3, 20 >,
                       FullIntegration,
                       GradientEnhancedFiniteStrainDisplacementElement< 3, 20 >::SectionType::Solid >() );

  const static bool GC3D20RUL_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "GC3D20RUL",
                     makeFactoryFunction<
                       GradientEnhancedFiniteStrainDisplacementElement< 3, 20 >,
                       ReducedIntegration,
                       GradientEnhancedFiniteStrainDisplacementElement< 3, 20 >::SectionType::Solid >() );

} // namespace
