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
#include "Eigen/Sparse"
#include "Marmot/MarmotElement.h"
#include "Marmot/MarmotElementProperty.h"
#include <functional>
#include <memory>

/**
 * @class MarmotElementSpatialWrapper
 * @brief Wrapper that embeds a lower-dimensional child element (e.g.\ a truss) into a
 *        higher-dimensional ambient space (2-D or 3-D).
 *
 * The projection transformation is constructed automatically from the supplied
 * nodal coordinates.  The child element is created via a user-provided factory
 * functor so that the wrapper remains independent of the concrete child type.
 */
class MarmotElementSpatialWrapper : public MarmotElement {

public:
  const int nDim;                                ///< Number of spatial dimensions of the ambient space.
  const int nDimChild;                           ///< Number of spatial dimensions of the child element.
  const int nNodes;                              ///< Number of nodes shared by parent and child element.
  const int nRhsChild;                           ///< Size of the child element's right-hand-side vector.
  const Eigen::Map< const Eigen::VectorXi >
            rhsIndicesToBeProjected;             ///< Indices in the child RHS vector that need projection.
  const int projectedSize;                       ///< Number of projected DOFs (child-element dimension).
  const int unprojectedSize;                     ///< Number of unprojected DOFs (ambient-space dimension).

  std::unique_ptr< MarmotElement > childElement; ///< Owned child element instance.
  Eigen::MatrixXd                  T;            ///< Coordinate transformation matrix from child to parent space.
  Eigen::MatrixXd                  P;            ///< Projection matrix mapping parent DOFs to child DOFs.
  Eigen::MatrixXd                  projectedCoordinates; ///< Nodal coordinates expressed in the child (local) frame.

  /**
   * @brief Construct the spatial wrapper.
   * @param[in] nDim                      Number of spatial dimensions of the ambient space.
   * @param[in] nChildDim                 Number of spatial dimensions of the child element.
   * @param[in] nNodes                    Number of nodes.
   * @param[in] sizeRhsChild              Size of the child element's right-hand-side vector.
   * @param[in] rhsIndicesToBeWrapped_    Array of child-RHS indices to be projected.
   * @param[in] nRhsIndicesToBeWrapped    Length of @p rhsIndicesToBeWrapped_.
   * @param[in] childElement              Owning pointer to the child element instance.
   */
  MarmotElementSpatialWrapper( int                              nDim,
                               int                              nChildDim,
                               int                              nNodes,
                               int                              sizeRhsChild,
                               const int                        rhsIndicesToBeWrapped_[],
                               int                              nRhsIndicesToBeWrapped,
                               std::unique_ptr< MarmotElement > childElement );

  /// @copydoc MarmotElement::getNumberOfRequiredStateVars
  int getNumberOfRequiredStateVars();

  /// @copydoc MarmotElement::getNodeFields
  std::vector< std::vector< std::string > > getNodeFields();

  /// @copydoc MarmotElement::getDofIndicesPermutationPattern
  std::vector< int > getDofIndicesPermutationPattern();

  /// @copydoc MarmotElement::getNNodes
  int getNNodes();

  /// @copydoc MarmotElement::getNSpatialDimensions
  int getNSpatialDimensions();

  /// @copydoc MarmotElement::getNDofPerElement
  int getNDofPerElement();

  /// @copydoc MarmotElement::getElementShape
  std::string getElementShape();

  /// @copydoc MarmotElement::assignStateVars
  void assignStateVars( double* stateVars, int nStateVars );

  /// @copydoc MarmotElement::assignProperty(const ElementProperties&)
  void assignProperty( const ElementProperties& property );

  /// @copydoc MarmotElement::assignProperty(const MarmotMaterialSection&)
  void assignProperty( const MarmotMaterialSection& property );

  /// @copydoc MarmotElement::assignNodeCoordinates
  void assignNodeCoordinates( const double* coordinates );

  /// @copydoc MarmotElement::initializeYourself
  void initializeYourself();

  /**
   * @brief Perform element computations with coordinate transformation.
   * @param[in]  QTotal  Total dof vector in the ambient (parent) space.
   * @param[in]  dQ      Incremental dof vector in the ambient space.
   * @param[out] Pe      Internal force vector in the ambient space.
   * @param[out] Ke      Stiffness matrix in the ambient space.
   * @param[in]  time    Current time.
   * @param[in]  dT      Time step size.
   */
  void computeKernels( const double* QTotal, const double* dQ, double* Pe, double* Ke, double time, double dT );

  /// @copydoc MarmotElement::setInitialConditions
  void setInitialConditions( StateTypes state, const double* values );

  /**
   * @brief Compute contribution from distributed surface loads with coordinate transformation.
   * @param[in]  loadType    Type of distributed load.
   * @param[out] P           External load vector in the ambient space.
   * @param[out] K           Load stiffness matrix in the ambient space.
   * @param[in]  elementFace Index of element face.
   * @param[in]  load        Applied load values.
   * @param[in]  QTotal      Total dof vector.
   * @param[in]  time        Current time.
   * @param[in]  dT          Time step size.
   */
  void computeDistributedLoad( DistributedLoadTypes loadType,
                               double*              P,
                               double*              K,
                               int                  elementFace,
                               const double*        load,
                               const double*        QTotal,
                               double               time,
                               double               dT );

  /**
   * @brief Compute body force contribution with coordinate transformation.
   * @param[out] P      External load vector in the ambient space.
   * @param[out] K      Load stiffness matrix in the ambient space.
   * @param[in]  load   Applied body force values.
   * @param[in]  QTotal Total dof vector.
   * @param[in]  time   Current time.
   * @param[in]  dT     Time step size.
   */
  void computeBodyForce( double* P, double* K, const double* load, const double* QTotal, double time, double dT );

  /// @copydoc MarmotElement::getStateView
  StateView getStateView( const std::string& stateName, int quadraturePoint );

  /// @copydoc MarmotElement::getCoordinatesAtCenter
  std::vector< double > getCoordinatesAtCenter();

  /// @copydoc MarmotElement::getCoordinatesAtQuadraturePoints
  std::vector< std::vector< double > > getCoordinatesAtQuadraturePoints();

  /// @copydoc MarmotElement::getNumberOfQuadraturePoints
  int getNumberOfQuadraturePoints();
};
