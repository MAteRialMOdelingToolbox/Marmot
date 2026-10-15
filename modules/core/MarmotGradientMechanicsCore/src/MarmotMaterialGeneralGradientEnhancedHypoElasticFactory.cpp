#include "Marmot/MarmotMaterialGeneralGradientEnhancedHypoElasticFactory.h"
#include "Marmot/MarmotJournal.h"

namespace Marmot::Registry {


template <>
MaterialGradientEnhancedHypoElasticFactory< 1 >::MaterialFactoryMap& MaterialGradientEnhancedHypoElasticFactory<
  1 >::materialFactoryFunctionByName()
{
  static MaterialFactoryMap map;
  return map;
}

template <>
Marmot::Materials::MaterialGradientEnhancedHypoElastic< 1 >* MaterialGradientEnhancedHypoElasticFactory<
  1 >::createMaterial( const std::string& materialName,
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

template <>
MaterialGradientEnhancedHypoElasticFactory< 2 >::MaterialFactoryMap& MaterialGradientEnhancedHypoElasticFactory<
  2 >::materialFactoryFunctionByName()
{
  static MaterialFactoryMap map;
  return map;
}

template <>
Marmot::Materials::MaterialGradientEnhancedHypoElastic< 2 >* MaterialGradientEnhancedHypoElasticFactory<
  2 >::createMaterial( const std::string& materialName,
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

template <>
MaterialGradientEnhancedHypoElasticFactory< 6 >::MaterialFactoryMap& MaterialGradientEnhancedHypoElasticFactory<
  6 >::materialFactoryFunctionByName()
{
  static MaterialFactoryMap map;
  return map;
}

template <>
Marmot::Materials::MaterialGradientEnhancedHypoElastic< 6 >* MaterialGradientEnhancedHypoElasticFactory<
  6 >::createMaterial( const std::string& materialName,
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

} // namespace Marmot::Registry
