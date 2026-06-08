/* ---------------------------------------------------------------------
 *                                       _
 *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_
 * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
 * | | | | | | (_| | |  | | | | | | (_) | |_
 * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
 *
 * Unit of Strength of Materials and Structural Analysis
 * University of Innsbruck,
 * 2020 - today
 *
 * festigkeitslehre@uibk.ac.at
 *
 * Alexandros Stathas alexandros.stathas@boku.ac.at
 * Matthias Neuner matthias.neuner@uibk.ac.at
 *
 * This file is part of the MAteRialMOdellingToolbox (marmot).
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of marmot.
 * ---------------------------------------------------------------------
 */

#pragma once

#include "Fastor/Fastor.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotStateHelpers.h"
#include "Marmot/MarmotTypedefs.h"

#include <cassert>
#include <functional>
#include <string>
#include <unordered_map>

/**
 *
 * Abstract base class for hypoelastic interface materials.
 *
 * This class remains an independent interface-material base class
 * because interface materials have their own stress-update signature.
 */
class MarmotInterfaceMaterialHypoElastic {

protected:
  const double* materialProperties;
  const int     nMaterialProperties;

public:
  const int materialNumber;

  MarmotInterfaceMaterialHypoElastic( const double* matProperties_, int nMaterialProperties_, int materialNumber_ )
    : materialProperties( matProperties_ ),
      nMaterialProperties( nMaterialProperties_ ),
      materialNumber( materialNumber_ )
  {
  }

  /// Default destructor
  virtual ~MarmotInterfaceMaterialHypoElastic() = default;

  /// Layout of the state variables
  MarmotStateLayoutDynamic stateLayout;

  /// Characteristic element length
  double characteristicElementLength;

  /**
   * Set the characteristic element length at the considered quadrature point.
   * It is needed for the regularization of materials with softening behavior
   * based on the mesh-adjusted softening modulus.
   *
   * @param[in] length characteristic length; will be assigned to
   * @ref characteristicElementLength
   */
  void setCharacteristicElementLength( double length );

  struct State {
    double* force;
    double* surfaceStress;
    double* stateVars;
  };

  struct Tangents {
    double* Q_ij;
    double* Z_ijkl;
    double* H_ijk;
    double* Y_ijkl;
  };

  struct Deformation {
    const double* dU;
    const double* dSurfaceStrain;
    const double* normal;
  };

  struct TimeIncrement {
    const double* timeOld;
    double        dT;
    double&       pNewDT;
  };

  /**
   * For a given interface displacement jump increment and surface strain
   * increment, compute the conjugate interface quantities and algorithmic
   * tangent terms.
   */
  virtual void computeStress( State&               state,
                              Tangents&            tangents,
                              const Deformation&   deformation,
                              const TimeIncrement& timeIncrement ) = 0;

  /**
   * @brief Initialize the layout of the state variables.
   *
   * This method has to be implemented in derived classes.
   *
   * @warning This method has to be called in the constructor of the derived class.
   */
  virtual void initializeStateLayout() = 0;

  /**
   * @brief Get a view to the state variables.
   *
   * @param stateName Name of the state variable.
   * @param stateVars Pointer to the state variable array.
   * @return StateView to access the requested state variable.
   */
  StateView getStateView( const std::string& stateName, double* stateVars ) const
  {
    return stateLayout.getStateView( stateVars, stateName );
  }

  /**
   * @brief Get the total number of required state variables.
   *
   * @return Total number of required state variables.
   */
  virtual int getNumberOfRequiredStateVars() const { return stateLayout.totalSize(); }

  /**
   * @brief Initialize the state variables at a material point.
   *
   * The default implementation initializes all state variables to zero.
   */
  virtual void initializeYourself( double* stateVars, int nStateVars )
  {
    for ( int i = 0; i < nStateVars; ++i ) {
      stateVars[i] = 0.0;
    }
  }

  virtual double getDensity() { return -1; }
};

namespace MarmotLibrary {

  /**
   * @class MarmotInterfaceMaterialHypoElasticFactory
   * @brief Factory class for creating hypoelastic interface-material instances by name.
   */
  class MarmotInterfaceMaterialHypoElasticFactory {
  public:
    using materialFactoryFunction = std::function<
      MarmotInterfaceMaterialHypoElastic*( const double* materialProperties,
                                           int           nMaterialProperties,
                                           int           materialNumber ) >;

    MarmotInterfaceMaterialHypoElasticFactory() = delete;
    static MarmotInterfaceMaterialHypoElastic* createMaterial( const std::string& materialName,
                                                               const double*      materialProperties,
                                                               int                nMaterialProperties,
                                                               int                materialNumber )
    {
      auto& map = materialFactoryFunctionByName();
      auto  it  = map.find( materialName );

      if ( it == map.end() ) {
        return nullptr;
      }

      return it->second( materialProperties, nMaterialProperties, materialNumber );
    }

    template < class T >
    static bool registerMaterial( const std::string& materialName )
    {
      auto& map = materialFactoryFunctionByName();

      assert( map.find( materialName ) == map.end() && "Interface material already registered!" );

      map[materialName] = []( const double* materialProperties,
                              int           nMaterialProperties,
                              int           materialNumber ) -> MarmotInterfaceMaterialHypoElastic* {
        return new T( materialProperties, nMaterialProperties, materialNumber );
      };

      return true;
    }

  private:
    using MaterialFactoryMap = std::unordered_map< std::string, materialFactoryFunction >;

    static MaterialFactoryMap& materialFactoryFunctionByName()
    {
      static MaterialFactoryMap map;
      return map;
    }
  };

} // namespace MarmotLibrary
