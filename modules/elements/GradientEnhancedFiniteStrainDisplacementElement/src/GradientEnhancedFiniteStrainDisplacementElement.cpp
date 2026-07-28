#include "Marmot/GradientEnhancedFiniteStrainDisplacementElement.h"
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

  // 2D plane strain, 8-node
  const static bool GCPE8UL_isRegistered = MarmotElementFactory::registerElement(
    "GCPE8UL",
    makeFactoryFunction< GradientEnhancedFiniteStrainDisplacementElement< 2, 8 >,
                         FullIntegration,
                         GradientEnhancedFiniteStrainDisplacementElement< 2, 8 >::PlaneStrain >() );

  const static bool GCPE8RUL_isRegistered = MarmotElementFactory::registerElement(
    "GCPE8RUL",
    makeFactoryFunction< GradientEnhancedFiniteStrainDisplacementElement< 2, 8 >,
                         ReducedIntegration,
                         GradientEnhancedFiniteStrainDisplacementElement< 2, 8 >::PlaneStrain >() );

  // 3D 8-node
  const static bool GC3D8UL_isRegistered = MarmotElementFactory::registerElement(
    "GC3D8UL",
    makeFactoryFunction< GradientEnhancedFiniteStrainDisplacementElement< 3, 8 >,
                         FullIntegration,
                         GradientEnhancedFiniteStrainDisplacementElement< 3, 8 >::SectionType::Solid >() );

  // 3D 20-node
  const static bool GC3D20UL_isRegistered = MarmotElementFactory::registerElement(
    "GC3D20UL",
    makeFactoryFunction< GradientEnhancedFiniteStrainDisplacementElement< 3, 20 >,
                         FullIntegration,
                         GradientEnhancedFiniteStrainDisplacementElement< 3, 20 >::SectionType::Solid >() );

  const static bool GC3D20RUL_isRegistered = MarmotElementFactory::registerElement(
    "GC3D20RUL",
    makeFactoryFunction< GradientEnhancedFiniteStrainDisplacementElement< 3, 20 >,
                         ReducedIntegration,
                         GradientEnhancedFiniteStrainDisplacementElement< 3, 20 >::SectionType::Solid >() );

} // namespace Marmot::Elements::Registration
