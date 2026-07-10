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

#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotStateHelpers.h"

#include <memory>
#include <string>
#include <vector>

/**
 *
 * Abstract base class for hypoelastic interface materials.
 *
 * This class remains an independent interface-material base class
 * because interface materials have their own stress-update signature.
 */
class MarmotInterfaceMaterialHypoElastic {

protected:
  const double*                                materialProperties;
  const int                                    nMaterialProperties;
  double                                       h = 0.0;
  std::vector< double >                        baseMaterialProperties;
  std::unique_ptr< MarmotMaterialHypoElastic > baseMaterial;

public:
  using TensorMap3d    = Marmot::FastorStandardTensors::TensorMap3d;
  using TensorMap33d   = Marmot::FastorStandardTensors::TensorMap33d;
  using TensorMap333d  = Marmot::FastorStandardTensors::TensorMap333d;
  using TensorMap3333d = Marmot::FastorStandardTensors::TensorMap3333d;
  using TensorMap6d    = Marmot::FastorStandardTensors::TensorMap6d;
  using TensorMap18d   = Marmot::FastorStandardTensors::TensorMap18d;

  const int materialNumber;

  MarmotInterfaceMaterialHypoElastic( const double* matProperties_, int nMaterialProperties_, int materialNumber_ )
    : materialProperties( matProperties_ ),
      nMaterialProperties( nMaterialProperties_ ),
      h( 0.0 ),
      baseMaterialProperties(),
      baseMaterial(),
      materialNumber( materialNumber_ )
  {
  }

  MarmotInterfaceMaterialHypoElastic( const std::string& materialName,
                                      const double*      matProperties_,
                                      int                nMaterialProperties_,
                                      int                materialNumber_ );

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
    TensorMap3d  force;
    TensorMap33d surfaceStress;
    double*      stateVars;
  };

  struct Tangents {
    TensorMap33d   Q_ij;
    TensorMap3333d Z_ijkl;
    TensorMap333d  H_ijk;
    TensorMap3333d Y_ijkl;
  };

  struct Deformation {
    TensorMap6d  dU;
    TensorMap18d dSurfaceStrain;
    TensorMap3d  normal;

    // Fastor's const TensorMap cannot be used with slicing and norm operations.
    // These views are therefore mutable types but are exposed through const Deformation&.
    Deformation( const double* dU_, const double* dSurfaceStrain_, const double* normal_ )
      : dU( const_cast< double* >( dU_ ) ),
        dSurfaceStrain( const_cast< double* >( dSurfaceStrain_ ) ),
        normal( const_cast< double* >( normal_ ) )
    {
    }
  };

  struct TimeIncrement {
    double timeOld;
    double dT;
  };

  /**
   * For a given interface displacement jump increment and surface strain
   * increment, compute the conjugate interface quantities and algorithmic
   * tangent terms.
   */
  virtual void computeStress( State&               state,
                              Tangents&            tangents,
                              const Deformation&   deformation,
                              const TimeIncrement& timeIncrement );

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
  virtual void initializeYourself( double* stateVars, int nStateVars );

  virtual double getDensity();
};
