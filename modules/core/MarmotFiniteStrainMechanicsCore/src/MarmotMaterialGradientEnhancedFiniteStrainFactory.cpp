#include "Marmot/MarmotMaterialGradientEnhancedFiniteStrainFactory.h"
#include "Marmot/MarmotJournal.h"

namespace Marmot::Registry {

  template < int nNonlocal >
  typename MaterialGradientEnhancedFiniteStrainFactory< nNonlocal >::Material* MaterialGradientEnhancedFiniteStrainFactory<
    nNonlocal >::createMaterial( const std::string& materialName,
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

  template < int nNonlocal >
  typename MaterialGradientEnhancedFiniteStrainFactory< nNonlocal >::MaterialFactoryMap& MaterialGradientEnhancedFiniteStrainFactory<
    nNonlocal >::materialFactoryFunctionByName()
  {
    static MaterialFactoryMap map;
    return map;
  }

  template class MaterialGradientEnhancedFiniteStrainFactory< 1 >;
  template class MaterialGradientEnhancedFiniteStrainFactory< 6 >;

} // namespace Marmot::Registry
