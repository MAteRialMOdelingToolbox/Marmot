/* ---------------------------------------------------------------------
 *                                       _
 *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_
 * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
 * | | | | | | (_| | |  | | | | | | (_) | |_
 * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
 *
 * Unit of Strength of Materials and Structural Analysis
 * University of Innsbruck
 * 2020 - today
 *
 * festigkeitslehre@uibk.ac.at
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
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMaterialPoint.h"
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace Marmot::Discretizations::Cells {

  /**
   * @class Cell
   * @brief Abstract interface of a background cell of the material point method (MPM).
   *
   * @details A cell is a fixed background element with @f$ n_\mathrm{nodes} @f$ nodes that carries the
   * nodal fields (e.g. the displacement) and the interpolation @f$ N_A(\boldsymbol{X}) @f$. It does not own
   * any quadrature points: in every connectivity update, the host framework assigns to each cell the
   * material points (Marmot::Discretizations::MaterialPoints::MaterialPoint) that currently lie inside it, and the cell
   * integrates over those material points, e.g. for the internal force
   * @f[
   *   \boldsymbol{f}^{\mathrm{int}}_A = \sum_{p \in \text{cell}} \boldsymbol{f}_{A}(\text{MP}_p)\,V_p ,
   * @f]
   * with the material point volume @f$ V_p @f$. The concrete kernel (which stress, which volume, which
   * configuration) is defined by the implementing cell, see e.g. Marmot::Discretizations::Cells::DisplacementCell.
   *
   * **Dof layout.** The cell vectors (Q, fInt, fExt) have getNDofPerCell() entries, the matrices are
   * dense with getNDofPerCell() @f$ \times @f$ getNDofPerCell() entries. getNodeFields() lists for each node
   * the names of the fields it carries; getDofIndicesPermutationPattern() is the permutation that maps
   * the node-major dof order of the host framework to the cell-internal order (for the cells in Marmot:
   * blocked, i.e. all dofs of the first field, then all dofs of the second field; within a field
   * node-major, component fastest).
   *
   * **Call sequence of the host framework** (e.g. EdelweissMeshfree), per time increment:
   *  1. Marmot::Discretizations::MaterialPoints::MaterialPoint::prepareYourself() on all material points,
   *  2. isCoordinateInCell() / getBoundingBox() to find the cells hit by the material points, then
   *     assignMaterialPoints() on every active cell (connectivity update),
   *  3. in every Newton iteration: Marmot::Discretizations::MaterialPoints::MaterialPoint::prepareYourself() (resets
   * the accumulated increment), interpolateFieldsToMaterialPoints() on all active cells,
   *     Marmot::Discretizations::MaterialPoints::MaterialPoint::computeYourself() on all material points, then
   *     computeMaterialPointKernels(), computeBodyLoad(), computeDistributedLoad() (and the inertia
   *     methods for dynamic analyses) on all active cells,
   *  4. after convergence, Marmot::Discretizations::MaterialPoints::MaterialPoint::acceptStateAndPosition() on all
   * material points.
   *
   * Instances are created by name via Marmot::Registry::CellFactory.
   */
  class Cell {

  public:
    /// Virtual destructor; cells are owned through Cell pointers.
    virtual ~Cell() = default;

    /**
     * @brief Names of the fields carried by each node.
     * @return One entry per node, each a list of field names (e.g. @c "displacement").
     */
    virtual const std::vector< std::vector< std::string > >& getNodeFields() const = 0;

    /**
     * @brief Permutation from the host's node-major dof order to the cell-internal dof order.
     * @return Vector of length getNDofPerCell(); entry @c i is the node-major index of cell dof @c i.
     */
    virtual const std::vector< int >& getDofIndicesPermutationPattern() const = 0;

    /**
     * @brief Number of cell nodes.
     * @return @f$ n_\mathrm{nodes} @f$.
     */
    virtual int getNNodes() const = 0;

    /**
     * @brief Total number of dofs of the cell.
     * @return Size of the cell vectors (and number of rows/columns of the cell matrices).
     */
    virtual int getNDofPerCell() const = 0;

    /**
     * @brief Shape of the cell for output.
     * @return Ensight Gold shape name, e.g. @c "quad4" or @c "hexa8".
     */
    virtual std::string getCellShape() const = 0;

    /**
     * @brief Checks whether a point lies inside the cell.
     * @param[in] coordinates Point coordinates (nDim values).
     * @return @c true if the point is inside the cell.
     */
    virtual bool isCoordinateInCell( const double* coordinates ) const = 0;

    /**
     * @brief Axis-aligned bounding box of the cell.
     * @param[out] boundingBoxMin Lower corner (nDim values).
     * @param[out] boundingBoxMax Upper corner (nDim values).
     */
    virtual void getBoundingBox( double* boundingBoxMin, double* boundingBoxMax ) const = 0;

    /**
     * @brief Assigns the material points that currently lie in the cell.
     * @details Replaces any previous assignment. Implementations typically evaluate and store the
     * interpolation (and its gradient) at each material point here, so this must be called again whenever
     * the material points have moved (connectivity update).
     * @param[in] materialPoints Material points located in the cell (not owned by the cell).
     * @throws std::invalid_argument (in the Marmot cells) if a material point is of the wrong type.
     */
    virtual void assignMaterialPoints(
      const std::vector< Marmot::Discretizations::MaterialPoints::MaterialPoint* >& materialPoints ) = 0;

    /**
     * @brief Optional: sets the point at which a cell-level projection of the cell (e.g. the volume change of an
     * F-bar cell) is evaluated, instead of the cell's own choice, until the next assignMaterialPoints().
     * @details Meant for cells that share one interpolation, e.g. the cells of an aggregate in the aggregated material
     * point method, whose boundary cells carry the extension of the field of an interior cell: evaluating the
     * projection of all of them at one common point imposes one constraint per aggregate instead of one per cell. The
     * point may lie outside the cell (the cell's basis is then extended). A cell without such a projection ignores
     * the point and returns false (the default).
     * @param[in] coordinates The point, in physical coordinates (nDim values).
     * @return Whether the cell uses the point.
     */
    virtual bool setProjectionPoint( [[maybe_unused]] const double* coordinates ) { return false; }

    /**
     * @name Optional: a projection over a patch of cells
     * @brief A scalar that the cells of a patch share, e.g. the volume change of an F-bar-patch formulation, computed
     * by the host in two passes.
     *
     * @details The patch value is the weighted mean of the contributions of the material points of all cells of the
     * patch,
     * @f[
     *   \bar{v} = \frac{\sum_c v_c(\boldsymbol{Q}_c)}{\sum_c w_c}, \qquad v_c = \sum_{p \in c} w_p\, v_p, \quad
     *   w_c = \sum_{p \in c} w_p ,
     * @f]
     * which couples the degrees of freedom of all cells of the patch. Per evaluation, the host
     *  1. calls computeProjectionContribution() on every cell of the patch (before
     *     interpolateFieldsToMaterialPoints(); it changes no material point), forms @f$ \bar{v} @f$,
     *  2. passes it to every cell with setProjectedValue(), then interpolates, computes the material points and the
     *     kernels as usual: computeMaterialPointKernels() (and the loads) return the residual of the cell and its
     *     tangent at fixed @f$ \bar{v} @f$,
     *  3. assembles the coupling of the patch, the outer products
     *     @f$ \partial \boldsymbol{R}_c/\partial \bar{v} \otimes \partial v_d/\partial \boldsymbol{Q}_d / \sum w @f$
     *     for all cells @f$ c, d @f$ of the patch, with getProjectionSensitivity().
     *
     * A cell on its own is a patch of one cell. The default implementations decline (hasPatchProjection() is false)
     * or throw.
     * @{
     */

    /**
     * @brief Whether the cell supports a projection over a patch of cells.
     * @return false (the default).
     */
    virtual bool hasPatchProjection() const { return false; }

    /**
     * @brief The contribution of the cell's material points to the patch value at the dof vector @p Q, and its
     * derivative; changes no material point.
     * @param[in]  Q         Cell dof vector (getNDofPerCell() values, the same as for
     *                       interpolateFieldsToMaterialPoints()).
     * @param[out] value     @f$ v_c = \sum_p w_p v_p @f$.
     * @param[out] weight    @f$ w_c = \sum_p w_p @f$ (independent of @p Q).
     * @param[out] dValue_dQ @f$ \partial v_c/\partial \boldsymbol{Q} @f$ (getNDofPerCell() values), overwritten.
     * @throws std::invalid_argument if not supported (default).
     */
    virtual void computeProjectionContribution( [[maybe_unused]] const double* Q,
                                                [[maybe_unused]] double*       value,
                                                [[maybe_unused]] double*       weight,
                                                [[maybe_unused]] double*       dValue_dQ ) const
    {
      throw std::invalid_argument( Marmot::Utilities::MakeString() << __PRETTY_FUNCTION__ << ": not supported" );
    }

    /**
     * @brief Sets the patch value @f$ \bar{v} @f$ the following evaluations use, and switches the cell to the patch
     * projection until the next assignMaterialPoints().
     * @param[in] value The patch value.
     * @throws std::invalid_argument if not supported (default).
     */
    virtual void setProjectedValue( [[maybe_unused]] double value )
    {
      throw std::invalid_argument( Marmot::Utilities::MakeString() << __PRETTY_FUNCTION__ << ": not supported" );
    }

    /**
     * @brief The derivative of the residual contributions of the cell (computeMaterialPointKernels() and
     * computeDistributedLoad(), as evaluated since the last interpolateFieldsToMaterialPoints()) w.r.t. the patch
     * value, in the sign convention of the tangents.
     * @param[out] dR_dValue getNDofPerCell() values, overwritten.
     * @throws std::invalid_argument if not supported (default).
     */
    virtual void getProjectionSensitivity( [[maybe_unused]] double* dR_dValue ) const
    {
      throw std::invalid_argument( Marmot::Utilities::MakeString() << __PRETTY_FUNCTION__ << ": not supported" );
    }
    /** @} */

    /**
     * @brief Integrates the internal force vector and its tangent over the assigned material points.
     * @details Uses the material point responses computed in the preceding
     * Marmot::Discretizations::MaterialPoints::MaterialPoint::computeYourself(). Results are added to the output
     * arrays.
     * @param[in]     Q        Cell dof vector (getNDofPerCell() values; the Marmot cells expect the
     *                         increment of the current step).
     * @param[in,out] fInt     Internal force vector (getNDofPerCell() values), accumulated.
     * @param[in,out] dFInt_dQ Tangent @f$ \partial \boldsymbol{f}^{\mathrm{int}}/\partial \boldsymbol{Q} @f$,
     *                         getNDofPerCell() @f$ \times @f$ getNDofPerCell(), accumulated.
     * @param[in]     timeNew  Time at the end of the increment.
     * @param[in]     dT       Time increment.
     */
    virtual void computeMaterialPointKernels( const double* Q,
                                              double*       fInt,
                                              double*       dFInt_dQ,
                                              double        timeNew,
                                              double        dT ) const = 0;

    /**
     * @brief Computes the lumped (diagonal) inertia of the cell.
     * @param[out] I Lumped mass vector (getNDofPerCell() values).
     * @throws std::invalid_argument if not implemented by the concrete cell (default).
     */
    virtual void computeLumpedInertia( [[maybe_unused]] double* I )
    {
      throw std::invalid_argument( Marmot::Utilities::MakeString() << __PRETTY_FUNCTION__ << "not yet implemented" );
    }

    /**
     * @brief Computes the consistent inertia (mass) matrix of the cell.
     * @param[out] I Mass matrix, getNDofPerCell() @f$ \times @f$ getNDofPerCell().
     * @throws std::invalid_argument if not implemented by the concrete cell (default).
     */
    virtual void computeConsistentInertia( [[maybe_unused]] double* I )
    {
      throw std::invalid_argument( Marmot::Utilities::MakeString() << __PRETTY_FUNCTION__ << "not yet implemented" );
    }

    /**
     * @brief Integrates a body load over the assigned material points.
     * @details The Marmot cells add the load with a negative sign, i.e. as its contribution to the
     * residual.
     * @param[in]     type    Load type id, from getSupportedBodyLoadTypes().
     * @param[in]     load    Load vector (e.g. force per unit volume, nDim values).
     * @param[in,out] fExt    External load vector (getNDofPerCell() values), accumulated.
     * @param[in,out] dExt_dQ Tangent of @p fExt w.r.t. the cell dofs, accumulated.
     * @param[in]     timeNew Time at the end of the increment.
     * @param[in]     dT      Time increment.
     * @throws std::invalid_argument (in the Marmot cells) for an unsupported load type.
     */
    virtual void computeBodyLoad( int           type,
                                  const double* load,
                                  double*       fExt,
                                  double*       dExt_dQ,
                                  double        timeNew,
                                  double        dT ) const = 0;

    /**
     * @brief Applies a distributed (surface) load carried by one assigned material point.
     * @details In MPM, boundary loads are attached to material points: only the material point with
     * number @p materialPointNumber contributes. The Marmot cells add the load with a negative sign, i.e.
     * as its contribution to the residual.
     * @param[in]     type                Load type id, from getSupportedDistributedLoadTypes().
     * @param[in]     surfaceID           Surface id of the load.
     * @param[in]     materialPointNumber Number of the material point that carries the load.
     * @param[in]     load                Load vector (nDim values).
     * @param[in,out] fExt                External load vector (getNDofPerCell() values), accumulated.
     * @param[in,out] dExt_dQ             Tangent of @p fExt w.r.t. the cell dofs, accumulated.
     * @param[in]     timeNew             Time at the end of the increment.
     * @param[in]     dT                  Time increment.
     * @throws std::invalid_argument (in the Marmot cells) for an unsupported load type.
     */
    virtual void computeDistributedLoad( int           type,
                                         int           surfaceID,
                                         int           materialPointNumber,
                                         const double* load,
                                         double*       fExt,
                                         double*       dExt_dQ,
                                         double        timeNew,
                                         double        dT ) const = 0;

    /**
     * @brief Evaluates the cell interpolation @f$ N_A @f$ at a point.
     * @param[out] vec         Shape function values (getNNodes() values).
     * @param[in]  coordinates Point coordinates (nDim values), inside the cell.
     */
    virtual void getInterpolationVector( double* vec, const double* coordinates ) const = 0;

    /**
     * @brief Evaluates the gradients @f$ \partial N_A/\partial X_i @f$ of the cell interpolation at a point, e.g. for
     * the normal derivatives on the faces of a ghost penalty stabilization. Optional: cells deriving from
     * CellWithGeometry implement it from their geometry (also outside of the cell, the gradients of the polynomial
     * extension of getInterpolationVector()); the default throws.
     * @param[out] dNdX        Gradients, nDim x getNNodes(), column-major (the nDim derivatives of node 0 first).
     * @param[in]  coordinates Point coordinates (nDim values).
     * @throws std::runtime_error if the cell does not provide its interpolation gradients.
     */
    virtual void getInterpolationGradient( [[maybe_unused]] double*       dNdX,
                                           [[maybe_unused]] const double* coordinates ) const
    {
      throw std::runtime_error( "Cell::getInterpolationGradient is not implemented for cell shape " + getCellShape() );
    }

    /**
     * @brief Body load types supported by the cell.
     * @return Map from upper-case load type name (e.g. @c "BODYFORCE") to the id passed to computeBodyLoad().
     */
    virtual const std::unordered_map< std::string, int >& getSupportedBodyLoadTypes() const = 0;

    /**
     * @brief Distributed load types supported by the cell.
     * @return Map from upper-case load type name (e.g. @c "PRESSURE") to the id passed to
     * computeDistributedLoad().
     */
    virtual const std::unordered_map< std::string, int >& getSupportedDistributedLoadTypes() const = 0;

    /**
     * @brief Interpolates the nodal fields to the assigned material points (grid-to-particle map).
     * @details The Marmot cells interpolate the dof increment and its gradient w.r.t. the reference
     * position and add them to the kinematic state of each material point; this accumulated increment is
     * reset in Marmot::Discretizations::MaterialPoints::MaterialPoint::prepareYourself().
     * @param[in] Q Cell dof vector (getNDofPerCell() values; the Marmot cells expect the increment of the
     *              current step).
     */
    virtual void interpolateFieldsToMaterialPoints( const double* Q ) const = 0;
  };

} // namespace Marmot::Discretizations::Cells
