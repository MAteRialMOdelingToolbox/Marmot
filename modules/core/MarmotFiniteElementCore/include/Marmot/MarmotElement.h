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

    Sigma11,                      ///< Normal stress component σ₁₁.
    Sigma22,                      ///< Normal stress component σ₂₂.
    Sigma33,                      ///< Normal stress component σ₃₃.
    HydrostaticStress,            ///< Hydrostatic (mean) stress.
    GeostaticStress,              ///< Geostatic in-situ stress state.
    MarmotMaterialStateVars,      ///< Internal material state variables.
    MarmotMaterialInitialization, ///< Trigger for material model initialization.
    HasEigenDeformation,          ///< Flag indicating presence of eigen (initial) deformation.
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
   * @param[in,out] stateVars Pointer to state variable array.
   * @param[in] nStateVars Number of state variables.
   */
  virtual void assignStateVars( double* stateVars, int nStateVars ) = 0;

  /**
   * @brief Assign element property set.
   * @param[in] property Element property object containing material, geometry, etc.
   */
  virtual void assignProperty( const ElementProperties& property );

  /**
   * @brief Assign material section property.
   * @param[in] property Material section definition (e.g. cross-sectional data).
   */
  virtual void assignProperty( const MarmotMaterialSection& property );

  /**
   * @brief Assign a single property of the element by name.
   * @param[in] propertyName Name of the property.
   * @param[in] properties Pointer to the array of property values.
   * @param[in] nProperties Number of values behind that pointer.
   * @note Default implementation throws an exception, as an element not overriding this
   * interface does not support any named properties.
   * @note The count is part of the interface rather than implied by the name: the caller is
   * typically a scripting layer handing over a user-supplied list, and an implementation reading
   * a fixed number of values from an unchecked pointer would read past the end of it silently.
   */
  virtual void assignProperty( const std::string& propertyName, const double* properties, int nProperties )
  {
    throw std::invalid_argument( MakeString()
                                 << __PRETTY_FUNCTION__ << ": unsupported named property '" << propertyName << "'" );
  };

  /**
   * @brief Get the names of all the valid properties of the element.
   * @return Vector of strings containing the property names.
   * @note Default implementation returns an empty vector, as an element not overriding this
   * interface does not expose any named properties.
   */
  virtual std::vector< std::string > getPropertyNames() const { return {}; };

  /**
   * @brief Assign nodal coordinates to element.
   * @param[in] coordinates Pointer to array of nodal coordinates.
   */
  virtual void assignNodeCoordinates( const double* coordinates ) = 0;

  /** @brief Initialize element state and internal variables. */
  virtual void initializeYourself() = 0;

  /**
   * @brief Apply initial conditions to the element.
   * @param[in] state State type to be set.
   * @param[in] values Array of initial values.
   */
  virtual void setInitialConditions( StateTypes state, const double* values ) = 0;

  /**
   * @brief Perform element computations (stiffness, residual, etc.).
   * @param[in] QTotal Total dof vector.
   * @param[in] dQ Incremental dof vector.
   * @param[out] Pint Internal force vector.
   * @param[out] K Stiffness matrix.
   * @param[in] time Current time.
   * @param[in] dT Time step size.
   */
  virtual void computeKernels( const double* QTotal,
                               const double* dQ,
                               double*       Pint,
                               double*       K,
                               double        time,
                               double        dT ) = 0;

  /**
   * @brief Perform element computations for explicit time integration.
   * @param[in] QTotal Total dof vector.
   * @param[in] dQ Incremental dof vector.
   * @param[out] Pint Internal force vector.
   * @param[in] time Current time.
   * @param[in] dT Time step size.
   *
   * @note Default implementation throws an exception.
   */
  virtual void computeKernelsExplicit( const double* QTotal, const double* dQ, double* Pint, double time, double dT )
  {
    throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << " not yet implemented" );
  };
  /**
   * @brief Compute contribution from distributed surface loads.
   * @param[in] loadType Type of load.
   * @param[out] Pext External load vector.
   * @param[out] K Stiffness matrix.
   * @param[in] elementFace Index of element face.
   * @param[in] load Applied load values.
   * @param[in] QTotal Total dof vector.
   * @param[in] time Current time.
   * @param[in] dT Time step size.
   */
  virtual void computeDistributedLoad( DistributedLoadTypes loadType,
                                       double*              Pext,
                                       double*              K,
                                       int                  elementFace,
                                       const double*        load,
                                       const double*        QTotal,
                                       double               time,
                                       double               dT ) = 0;

  /**
   * @brief Compute contribution from body forces.
   * @param[out] Pext External load vector.
   * @param[out] K Stiffness matrix.
   * @param[in] load Body force vector.
   * @param[in] QTotal Total displacement vector.
   * @param[in] time Current time.
   * @param[in] dT Time step size.
   */
  virtual void computeBodyForce( double*       Pext,
                                 double*       K,
                                 const double* load,
                                 const double* QTotal,
                                 double        time,
                                 double        dT ) = 0;

  /**
   * @brief Compute the lumped (diagonal) inertia of the element, over every field it carries.
   * @param[out] I Diagonal of the lumped inertia, in the element's dof order.
   * @details The coefficient of each field's SECOND time derivative: mass on the displacement
   * block, and, on a non-local block whose field has been given the named property "nonlocal micro
   * inertia", that micro-inertia. Zero on any non-local block that has not been -- carrying none is
   * what keeps that field first order in time; see computeLumpedDamping() for what integrates it in
   * that case.
   * @note Default implementation throws an exception.
   */
  virtual void computeLumpedInertia( double* I )
  {
    throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << " not yet implemented" );
  };

  /**
   * @brief Compute consistent inertia matrix.
   * @param[out] I Inertia matrix.
   * @note Default implementation throws an exception.
   */
  virtual void computeConsistentInertia( double* I )
  {
    throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << " not yet implemented" );
  };

  /**
   * @brief Compute critical time step for explicit dynamics.
   * @param[out] criticalTimeStep Suggested critical time step size.
   * @param[in] QTotal Total dof vector.
   * @note Default implementation throws an exception.
   */
  virtual void computeCriticalTimeStepForExplicitDynamics( double& criticalTimeStep, const double* QTotal )
  {
    throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << " not yet implemented" );
  };

  /**
   * @brief Compute internal energy of the element.
   * @param[out] internalEnergy Computed internal energy.
   * @note Default implementation throws an exception.
   */
  virtual void computeInternalEnergy( double& internalEnergy )
  {
    throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << " not yet implemented" );
  };

  /**
   * @brief Access element state at a quadrature point.
   * @param[in] stateName Name of the state variable.
   * @param[in] quadraturePoint Index of quadrature point.
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

  /**
   * @brief Compute the lumped (diagonal) damping of the element, over every field it carries.
   * @param[out] C Diagonal of the lumped damping, in the element's dof order.
   * @details The coefficient of each field's FIRST time derivative: zero on the displacement
   * block, where no device reports through this path, and the non-local viscosity on a non-local
   * block -- always, whether or not that field has been given a micro-inertia (see
   * computeLumpedInertia()). A first-order non-local field is integrated by this term alone; a
   * second-order one is damped by it. Unlike most of its siblings here this default does NOT
   * throw: carrying no damping is the ordinary answer for the overwhelming majority of elements,
   * not an unimplemented case, and the caller assembles over every element in the model. The
   * buffer is supplied zero-initialised, so a default that leaves it untouched reports exactly
   * that.
   *
   * @note Declared LAST, away from computeLumpedInertia() where it belongs by subject, because
   * adding a virtual function in the middle of this class shifts the vtable slot of every virtual
   * declared after it. Everything that links against Marmot has to be rebuilt when this header
   * changes either way -- but a stale library mixed with a fresh consumer then dispatches existing
   * calls to the wrong function, silently, instead of failing on the one call that is actually
   * new. Appending keeps that failure mode confined to callers of this method.
   */
  virtual void computeLumpedDamping( double* C ) { static_cast< void >( C ); };
};
