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
 * Matthias Neuner matthias.neuner@uibk.ac.at
 * Magdalena Schreter magdalena.schreter@uibk.ac.at
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
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotUtils.h"
#include <string>

/**
 * @class MarmotMaterialPoint
 * @brief Abstract interface of a material point of the material point method (MPM).
 *
 * @details A material point is a Lagrangian integration point that carries mass (volume and density),
 * position and the complete history (deformation, stress, material state variables) through the
 * simulation, while the background cells (MarmotCell) carry the nodal fields and are reset every
 * increment. A material point owns its constitutive model (assigned with assignMaterial()); its state
 * variables live in an externally owned array (assigned with assignStateVars()).
 *
 * **Life cycle** (as called by the host framework, e.g. EdelweissMeshfree):
 *  1. creation by name via MarmotLibrary::MarmotMaterialPointFactory, then assignMaterial(),
 *     getNumberOfRequiredStateVars() and assignStateVars(),
 *  2. initializeYourself() once, before the first increment,
 *  3. per iteration: prepareYourself(), then the cells interpolate the nodal increments to the material
 *     point (MarmotCell::interpolateFieldsToMaterialPoints()), then computeYourself(), which evaluates the
 *     constitutive response that the cells integrate in MarmotCell::computeMaterialPointKernels(),
 *  4. acceptStateAndPosition() after convergence of the increment, which commits the increment
 *     (e.g. updates position and deformation gradient).
 *
 * **State contract:** the state variable array is a *trial* copy. Before every iteration (i.e. before
 * prepareYourself()), the host restores it to the values committed by the last acceptStateAndPosition(), and the
 * cells interpolate the total increment since that state (not the Newton correction). Implementations may
 * therefore update their state in place during an iteration (e.g. add the interpolated increment of a field to its
 * committed value), since each iteration starts again from the committed state; a cutback is simply a restore.
 * EdelweissMeshfree implements this by copying the committed array into the trial array in prepareYourself() and
 * back in acceptStateAndPosition().
 *
 * The kinematic update from the cell (e.g. @c incrementDeformation) and the response/tangent quantities
 * read by the cell are not part of this interface: they are specific to each pair of concrete cell and
 * material point (see e.g. Marmot::MaterialPoints::DisplacementMaterialPoint), and the cell obtains them by
 * a @c dynamic_cast.
 */
class MarmotMaterialPoint {

public:
  /// Virtual destructor; material points are owned through MarmotMaterialPoint pointers.
  virtual ~MarmotMaterialPoint(){};

  /**
   * @brief Assigns the (externally owned) state variable array.
   * @param[in,out] stateVars  State variable array of size getNumberOfRequiredStateVars().
   * @param[in]     nStateVars Size of @p stateVars.
   */
  virtual void assignStateVars( double* stateVars, int nStateVars ) = 0;

  /**
   * @brief Creates and assigns the constitutive model.
   * @param[in] material Material section (material name and properties).
   */
  virtual void assignMaterial( const MarmotMaterialSection& material ) = 0;

  /**
   * @brief Initializes the state (e.g. deformation gradient and material state) before the first increment.
   */
  virtual void initializeYourself() = 0;

  /**
   * @brief Prepares a new iteration: resets the kinematic increment that the cells accumulate.
   * @param[in] timeNew Time at the end of the increment.
   * @param[in] dT      Time increment.
   */
  virtual void prepareYourself( double timeNew, double dT ) = 0;

  /**
   * @brief Evaluates the constitutive response for the current kinematic increment.
   * @param[in] timeNew Time at the end of the increment.
   * @param[in] dT      Time increment.
   */
  virtual void computeYourself( double timeNew, double dT ) = 0;

  /// Commits the converged increment (position, deformation). Does nothing by default.
  virtual void acceptStateAndPosition(){};

  /**
   * @brief Access to a named state (e.g. @c "stress" or a material state variable).
   * @param[in] stateName Name of the state.
   * @return View (location and size) into the state variable array.
   */
  virtual StateView getStateView( const std::string& stateName ) const = 0;

  /**
   * @brief Number of state variables required by the material point and its material.
   * @return Number of state variables; valid after assignMaterial().
   */
  virtual int getNumberOfRequiredStateVars() const = 0;

  /**
   * @brief Coordinates of the material point center in the last accepted state.
   * @details The increment of the current step is not included (it is committed in acceptStateAndPosition()),
   * so the cells locate the material point in the configuration at the beginning of the increment.
   * @param[out] coordinates getDimension() values.
   */
  virtual void getCoordinatesAtCenter( double* coordinates ) const = 0;

  /**
   * @brief Unique number of the material point.
   * @return Material point number.
   */
  virtual int getMaterialPointNumber() const = 0;

  /**
   * @brief Spatial dimension of the material point.
   * @return 2 or 3.
   */
  virtual int getDimension() const = 0;

  /**
   * @brief Number of vertices describing the material point domain.
   * @return Number of vertices (1 for a point).
   */
  virtual int getNumberOfVertices() const = 0;

  /**
   * @brief Shape of the material point for output.
   * @return Ensight Gold shape name, e.g. @c "point".
   */
  virtual std::string getMaterialPointShape() const = 0;

  /**
   * @brief Coordinates of the vertices in the last accepted state.
   * @param[out] coordinates getNumberOfVertices() @f$ \times @f$ getDimension() values, vertex by vertex.
   */
  virtual void getVertexCoordinates( double* coordinates ) const = 0;

  /**
   * @brief Displacement of the material point center in the last accepted state.
   * @param[out] displacement getDimension() values.
   */
  virtual void getCenterDisplacement( double* displacement ) const = 0;

  /**
   * @brief Volume in the reference configuration.
   * @return @f$ V_0 @f$.
   */
  virtual double getVolumeUndeformed() const = 0;

  /**
   * @brief Density in the reference configuration.
   * @return @f$ \rho_0 @f$.
   */
  virtual double getDensityUndeformed() const = 0;

  /**
   * @brief Sets an initial condition.
   * @param[in] conditionName Name of the initial condition.
   * @param[in] value         Values of the initial condition.
   */
  virtual void setInitialCondition( const std::string& conditionName, const double* value ) = 0;
};
