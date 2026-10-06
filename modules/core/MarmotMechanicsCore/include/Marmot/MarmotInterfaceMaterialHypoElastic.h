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
  const double*         materialProperties;     ///< material properties [E, nu, h, properties of the base material...]
  const int             nMaterialProperties;    ///< number of material properties
  const double          h;                      ///< thickness of the interface, the third material property
  std::vector< double > baseMaterialProperties; ///< properties handed to the base material
  std::unique_ptr< MarmotMaterialHypoElastic > baseMaterial; ///< hypoelastic material of the bulk

public:
  using Tensor3d    = Marmot::FastorStandardTensors::Tensor3d;    ///< first order tensor in 3D
  using Tensor33d   = Marmot::FastorStandardTensors::Tensor33d;   ///< second order tensor in 3D
  using Tensor333d  = Marmot::FastorStandardTensors::Tensor333d;  ///< third order tensor in 3D
  using Tensor3333d = Marmot::FastorStandardTensors::Tensor3333d; ///< fourth order tensor in 3D
  using Tensor6d    = Marmot::FastorStandardTensors::Tensor6d;    ///< vector of six components, one block per side
  using Tensor18d   = Marmot::FastorStandardTensors::Tensor18d;   ///< vector of 18 components, one block per side

  const int materialNumber; ///< number of the material

  /**
   * @brief Construct an interface material from a registered hypoelastic material.
   *
   * The material properties are [E, nu, h, properties of the base material...], where h is the thickness of the
   * interface. The base material receives E, nu and the properties after h.
   *
   * @param[in] materialName name of the base material, registered in the MarmotMaterialHypoElasticFactory
   * @param[in] matProperties_ material properties
   * @param[in] nMaterialProperties_ number of material properties, at least three
   * @param[in] materialNumber_ number of the material
   *
   * @throws std::invalid_argument if fewer than three properties are given or the base material is not registered
   */
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

  /// Interface response: holds the previous values on entry and the updated ones on exit.
  struct State {
    Tensor3d  force;         ///< traction on the interface
    Tensor33d surfaceStress; ///< surface stress
    double*   stateVars;     ///< pointer to the state variables
  };

  /// Algorithmic tangent terms, set by computeStress().
  struct Tangents {
    Tensor33d   Q_ij;   ///< tangent term acting on the displacement jump
    Tensor3333d Z_ijkl; ///< tangent term acting on the average surface gradient
    Tensor333d  H_ijk;  ///< coupling term between the displacement jump and the average surface gradient
    Tensor3333d Y_ijkl; ///< additional tangent term acting on the average surface gradient
  };

  /// Increment of the interface kinematics.
  struct Deformation {
    Tensor6d  dU;             ///< displacement jump increment of the two interface sides
    Tensor18d dSurfaceStrain; ///< surface displacement gradient increment of the two interface sides
    Tensor3d  normal;         ///< interface normal
  };

  /// Time at the beginning of the increment and time increment.
  struct TimeIncrement {
    double timeOld; ///< time at the beginning of the increment
    double dT;      ///< time increment
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

  /**
   * @brief Get the density of the material.
   *
   * @return the density of the base material
   */
  virtual double getDensity();
};
