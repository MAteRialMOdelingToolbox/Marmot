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
 * This class follows the structural pattern of MarmotMaterialHypoElastic in
 * Marmot v26.05, but remains an independent interface-material base class
 * because interface materials have their own stress-update signature.
 */
class MarmotMaterialHypoElasticInterface {

protected:
  const double* materialProperties;
  const int     nMaterialProperties;

public:
  const int materialNumber;

  MarmotMaterialHypoElasticInterface( const double* matProperties_, int nMaterialProperties_, int materialNumber_ )
    : materialProperties( matProperties_ ),
      nMaterialProperties( nMaterialProperties_ ),
      materialNumber( materialNumber_ )
  {
  }

  /// Default destructor
  virtual ~MarmotMaterialHypoElasticInterface() = default;

  /// Layout of the state variables
  MarmotStateLayoutDynamic stateLayout;

  using Tensor1D = Fastor::Tensor< double, 3 >;
  using Tensor2D = Fastor::Tensor< double, 3, 3 >;

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

  /**
   * For a given interface displacement jump increment and surface strain
   * increment, compute the conjugate interface quantities and algorithmic
   * tangent terms.
   *
   * @param[in,out] force              conjugate force due to displacement jump
   * @param[in,out] surface_stress     conjugate surface stress due to average surface strain
   * @param[in,out] H_inv_ij           tangent related to the linearized displacement jump
   * @param[in,out] Z_ijkl             tangent related to the linearized average surface strain
   * @param[in,out] H_inv_nF_ijk       tangent of the consistency terms
   * @param[in,out] Yn_H_inv_Fn_ijkl   tangent contribution for average surface strain
   * @param[in]     dU                 linearized displacement increment on top/bottom interface sides
   * @param[in]     dSurface_strain    linearized surface strain increment
   * @param[in]     normal             interface normal, positive toward the top side
   * @param[in]     timeOld            old pseudo-time
   * @param[in]     dT                 pseudo-time increment
   * @param[in,out] pNewDT             suggested new time increment
   */
  virtual void computeStress( double*       force,
                              double*       surface_stress,
                              double*       H_inv_ij,
                              double*       Z_ijkl,
                              double*       H_inv_nF_ijk,
                              double*       Yn_H_inv_Fn_ijkl,
                              const double* dU,
                              const double* dSurface_strain,
                              const double* normal,
                              const double* timeOld,
                              const double  dT,
                              double&       pNewDT ) = 0;

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

  /**
   * @brief Compatibility wrapper for interface materials ported from the old API.
   *
   * The old Marmot API used assignStateVars. Marmot v26.05 uses
   * initializeYourself/stateLayout. Keeping this wrapper allows the old
   * interface-material implementations to be ported incrementally.
   */
  virtual void assignStateVars( double* stateVars, int nStateVars )
  {
    initializeYourself( stateVars, nStateVars );
  }

  virtual double getDensity() { return -1; }
};

namespace MarmotLibrary {

  /**
   * @class MarmotMaterialHypoElasticInterfaceFactory
   * @brief Factory class for creating hypoelastic interface-material instances by name.
   *
   * This follows the v26.05 name-based registration style. It does not use
   * material registration numbers.
   */
  class MarmotMaterialHypoElasticInterfaceFactory {
  public:
    using materialFactoryFunction = std::function<
      MarmotMaterialHypoElasticInterface*( const double* materialProperties, int nMaterialProperties, int materialNumber ) >;

    MarmotMaterialHypoElasticInterfaceFactory() = delete;

    static MarmotMaterialHypoElasticInterface* createMaterial( const std::string& materialName,
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

      map[materialName] = []( const double* materialProperties, int nMaterialProperties, int materialNumber )
        -> MarmotMaterialHypoElasticInterface* {
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
