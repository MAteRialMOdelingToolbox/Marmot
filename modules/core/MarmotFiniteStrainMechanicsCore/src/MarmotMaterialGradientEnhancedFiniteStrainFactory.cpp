#include "Marmot/MarmotMaterialGradientEnhancedFiniteStrainFactory.h"
#include "Marmot/MarmotJournal.h"

namespace Marmot::Registry {


Marmot::Materials::MaterialGradientEnhancedFiniteStrain* MaterialGradientEnhancedFiniteStrainFactory::createMaterial(
  const std::string& materialName,
  const double*      materialProperties,
  int                nMaterialProperties,
  int                materialNumber )
{
  auto& map = materialFactoryFunctionByName();
  auto  it  = map.find( materialName );
  if ( it == map.end() ) {
    std::string reg = "Registered materials are: ";
    for ( const auto& pair : map ) {
      reg += pair.first + ", ";
    }
    throw std::invalid_argument( Marmot::Utilities::MakeString()
                                 << __PRETTY_FUNCTION__ << " Material " + materialName + " not registered!" + reg );
  }

  return it->second( materialProperties, nMaterialProperties, materialNumber );
}

MaterialGradientEnhancedFiniteStrainFactory::MaterialFactoryMap& MaterialGradientEnhancedFiniteStrainFactory::
  materialFactoryFunctionByName()
{
  static MaterialFactoryMap map;
  return map;
}

} // namespace Marmot::Registry
