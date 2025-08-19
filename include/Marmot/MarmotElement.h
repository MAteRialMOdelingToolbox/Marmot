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
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotUtils.h"
#include <stdexcept>
#include <string>
#include <vector>

/**
 * @class MarmotElement
 * @brief Abstract base class for finite elements in the Marmot framework.
 *
 * This class defines the generic interface for finite elements, including
 * methods for state variable handling, geometry, degrees of freedom,
 * initialization, loading, and numerical integration. Concrete element
 * implementations must override the pure virtual functions.
 */
class MarmotElement {

public:
  /** @brief Types of element state variables used in initialization and output. */
  enum StateTypes {
    Sigma11,                     ///< Stress component \f$\sigma_{11}\f$
    Sigma22,                     ///< Stress component \f$\sigma_{22}\f$
    Sigma33,                     ///< Stress component \f$\sigma_{33}\f$
    HydrostaticStress,           ///< Hydrostatic stress state
    GeostaticStress,             ///< Geostatic stress state
    MarmotMaterialStateVars,     ///< Generic material state variables
    MarmotMaterialInitialization ///< Initialization state for materials
  };

  /** @brief Types of distributed loads applicable to element boundaries. */
  enum DistributedLoadTypes {
    Pressure,        ///< Pressure load (normal to surface)
    SurfaceTorsion,  ///< Surface torsional load
    SurfaceTraction, ///< Surface traction vector
  };

  /** @brief Virtual destructor for safe polymorphic cleanup. */
  virtual ~MarmotElement();

  /** @return Number of state variables required by the element. */
  virtual int getNumberOfRequiredStateVars() = 0;

  /**
   * @brief Get the nodal field names (e.g. displacement, rotation).
   * @return A 2D vector of strings representing the fields per node.
   */
  virtual std::vector< std::vector< std::string > > getNodeFields() = 0;

  /**
   * @brief Get permutation pattern for degrees of freedom.
   * @return Vector of indices describing the permutation.
   */
  virtual std::vector< int > getDofIndicesPermutationPattern() = 0;

  /** @return Number of nodes in the element. */
  virtual int getNNodes() = 0;

  /** @return Number of spatial dimensions (2D/3D). */
  virtual int getNSpatialDimensions() = 0;

  /** @return Number of degrees of freedom per element. */
  virtual int getNDofPerElement() = 0;

  /** @return String describing the element shape in Ensight Gold notation (e.g. "quad4", "hexa8"). */
  virtual std::string getElementShape() = 0;

  /**
   * @brief Assign state variable array to element.
   * @param stateVars Pointer to state variable array.
   * @param nStateVars Number of state variables.
   */
  virtual void assignStateVars( double* stateVars, int nStateVars ) = 0;

  /**
   * @brief Assign element property set.
   * @param property Element property object containing material, geometry, etc.
   */
  virtual void assignProperty( const ElementProperties& property );

  /**
   * @brief Assign material section property.
   * @param property Material section definition (e.g. cross-sectional data).
   */
  virtual void assignProperty( const MarmotMaterialSection& property );

  /**
   * @brief Assign nodal coordinates to element.
   * @param coordinates Pointer to array of nodal coordinates.
   */
  virtual void assignNodeCoordinates( const double* coordinates ) = 0;

  /** @brief Initialize element state and internal variables. */
  virtual void initializeYourself() = 0;

  /**
   * @brief Apply initial conditions to the element.
   * @param state State type to be set.
   * @param values Array of initial values.
   */
  virtual void setInitialConditions( StateTypes state, const double* values ) = 0;

  /**
   * @brief Perform element computations (stiffness, residual, etc.).
   * @param QTotal Total dof vector.
   * @param dQ Incremental dof vector.
   * @param Pint Internal force vector (output).
   * @param K Stiffness matrix (output).
   * @param time Current time.
   * @param dT Time step size.
   * @param pNewdT Suggested new time step size (output).
   */
  virtual void computeYourself( const double* QTotal,
                                const double* dQ,
                                double*       Pint,
                                double*       K,
                                const double* time,
                                double        dT,
                                double&       pNewdT ) = 0;

  /**
   * @brief Compute contribution from distributed surface loads.
   * @param loadType Type of load.
   * @param Pext External load vector (output).
   * @param K Stiffness matrix (output).
   * @param elementFace Index of element face.
   * @param load Applied load values.
   * @param QTotal Total displacement vector.
   * @param time Current time.
   * @param dT Time step size.
   */
  virtual void computeDistributedLoad( DistributedLoadTypes loadType,
                                       double*              Pext,
                                       double*              K,
                                       int                  elementFace,
                                       const double*        load,
                                       const double*        QTotal,
                                       const double*        time,
                                       double               dT ) = 0;

  /**
   * @brief Compute contribution from body forces.
   * @param Pext External load vector (output).
   * @param K Stiffness matrix (output).
   * @param load Body force vector.
   * @param QTotal Total displacement vector.
   * @param time Current time.
   * @param dT Time step size.
   */
  virtual void computeBodyForce( double*       Pext,
                                 double*       K,
                                 const double* load,
                                 const double* QTotal,
                                 const double* time,
                                 double        dT ) = 0;

  /**
   * @brief Compute lumped inertia matrix.
   * @param I Inertia matrix (output).
   * @note Default implementation throws an exception.
   */
  virtual void computeLumpedInertia( double* I )
  {
    throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << " not yet implemented" );
  };

  /**
   * @brief Compute consistent inertia matrix.
   * @param I Inertia matrix (output).
   * @note Default implementation throws an exception.
   */
  virtual void computeConsistentInertia( double* I )
  {
    throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << " not yet implemented" );
  };

  /**
   * @brief Access element state at a quadrature point.
   * @param stateName Name of the state variable.
   * @param quadraturePoint Index of quadrature point.
   * @return View into state variable.
   */
  virtual StateView getStateView( const std::string& stateName, int quadraturePoint ) = 0;

  /**
   * @brief Get coordinates of element center.
   * @return Vector of coordinates at element centroid.
   */
  virtual std::vector< double > getCoordinatesAtCenter() = 0;

  /**
   * @brief Get coordinates of quadrature points.
   * @return 2D vector of coordinates at quadrature points.
   */
  virtual std::vector< std::vector< double > > getCoordinatesAtQuadraturePoints() = 0;

  /** @return Number of quadrature points used by the element. */
  virtual int getNumberOfQuadraturePoints() = 0;
};
