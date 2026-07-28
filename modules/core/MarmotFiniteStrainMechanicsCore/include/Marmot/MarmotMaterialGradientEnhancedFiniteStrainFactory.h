/* ---------------------------------------------------------------------
 *                                       _
 *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_
 * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
 * | | | | | | (_| | |  | | | | | | (_) | |_
 * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
 *
 * Unit of Strength of Materials and Structural Analysis
 * University of Innsbruck.
 *
 * Thomas Mader thomas.mader@boku.ac.at
 *
 * This file is part of the MAteRialMOdellingToolbox (marmot).
 * LGPL v2.1+, see LICENSE.md at the top level directory of marmot.
 * ---------------------------------------------------------------------
 */
#pragma once
#include "Marmot/MarmotMaterialGradientEnhancedFiniteStrain.h"
#include <cassert>
#include <functional>
#include <string>
#include <unordered_map>

namespace MarmotLibrary {

  /**
   * @class MarmotMaterialGradientEnhancedFiniteStrainFactory
   * @brief Factory for creating gradient-enhanced finite-strain material instances by name.
   */
  class MarmotMaterialGradientEnhancedFiniteStrainFactory {
  public:
    using materialFactoryFunction = std::function< MarmotMaterialGradientEnhancedFiniteStrain*(
      const double* materialProperties,
      int           nMaterialProperties,
      int           materialNumber ) >;

    MarmotMaterialGradientEnhancedFiniteStrainFactory() = delete;

    static MarmotMaterialGradientEnhancedFiniteStrain* createMaterial( const std::string& materialName,
                                                                       const double*      materialProperties,
                                                                       int                nMaterialProperties,
                                                                       int                materialNumber );

    template < class T >
    static bool registerMaterial( const std::string& materialName )
    {
      auto& map = materialFactoryFunctionByName();

      assert( map.find( materialName ) == map.end() && "Material already registered!" );

      map[materialName] = []( const double* materialProperties, int nMaterialProperties, int materialNumber )
        -> MarmotMaterialGradientEnhancedFiniteStrain* {
        return new T( materialProperties, nMaterialProperties, materialNumber );
      };
      return true;
    }

  private:
    using MaterialFactoryMap = std::unordered_map< std::string, materialFactoryFunction >;
    static MaterialFactoryMap& materialFactoryFunctionByName();
  };
} // namespace MarmotLibrary
