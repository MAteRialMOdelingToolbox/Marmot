#include "Marmot/MarmotMaterialGeneralGradientEnhancedHypoElasticFactory.h"
#include "Marmot/MarmotJournal.h"

using namespace Marmot::Factory;

template <>
MarmotMaterialGeneralGradientEnhancedHypoElasticFactory< 1 >::MaterialFactoryMap& MarmotMaterialGeneralGradientEnhancedHypoElasticFactory<
  1 >::materialFactoryFunctionByName()
{
  static MaterialFactoryMap map;
  return map;
}

template <>
Marmot::MarmotMaterialGeneralGradientEnhancedHypoElastic< 1 >* MarmotMaterialGeneralGradientEnhancedHypoElasticFactory<
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
    throw std::invalid_argument( Marmot::MakeString()
                                 << __PRETTY_FUNCTION__ << " Material " + materialName + " not registered!" + reg );
  }

  return it->second( materialProperties, nMaterialProperties, materialNumber );
}

template <>
MarmotMaterialGeneralGradientEnhancedHypoElasticFactory< 2 >::MaterialFactoryMap& MarmotMaterialGeneralGradientEnhancedHypoElasticFactory<
  2 >::materialFactoryFunctionByName()
{
  static MaterialFactoryMap map;
  return map;
}

template <>
Marmot::MarmotMaterialGeneralGradientEnhancedHypoElastic< 2 >* MarmotMaterialGeneralGradientEnhancedHypoElasticFactory<
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
    throw std::invalid_argument( Marmot::MakeString()
                                 << __PRETTY_FUNCTION__ << " Material " + materialName + " not registered!" + reg );
  }

  return it->second( materialProperties, nMaterialProperties, materialNumber );
}

template <>
MarmotMaterialGeneralGradientEnhancedHypoElasticFactory< 6 >::MaterialFactoryMap& MarmotMaterialGeneralGradientEnhancedHypoElasticFactory<
  6 >::materialFactoryFunctionByName()
{
  static MaterialFactoryMap map;
  return map;
}

template <>
Marmot::MarmotMaterialGeneralGradientEnhancedHypoElastic< 6 >* MarmotMaterialGeneralGradientEnhancedHypoElasticFactory<
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
    throw std::invalid_argument( Marmot::MakeString()
                                 << __PRETTY_FUNCTION__ << " Material " + materialName + " not registered!" + reg );
  }

  return it->second( materialProperties, nMaterialProperties, materialNumber );
}
