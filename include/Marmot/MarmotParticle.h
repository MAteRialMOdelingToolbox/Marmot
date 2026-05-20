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

#include "Marmot/MarmotMeshfreeKernelFunction.h"
#include "Marmot/MarmotUtils.h"
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace Marmot::Meshfree {

  /**
   * @brief The MarmotParticle class is the base class for all particles in the meshfree framework.
   * It provides the interface for the interaction with the meshfree framework.
   *
   * The MarmotParticle class is the base class for all particles in the meshfree framework.
   * It interacts with degrees of freedom (living on the nodes) via shape functions (MarmotMeshfreeShapeFunction).
   *
   * Similar to finite elements, the MarmotParticle class is responsible for the computation of the internal forces,
   * the body loads and the distributed loads.
   * Also, it may have a distinct shape (commonly a point).
   *
   * For the computation of the residual vector and the stiffness matrix, unlike elements and cells, no (field-)blocked
   * layout is possible. This is due to the fact that the number of nodes per particle is not fixed and may vary, and
   * accordingly, the dofIndicesPermutationPattern used to designate the structure of a blocked storage may vary.
   * Accordingly, for perfomance reasons, no dofIndicesPermutationPattern is provided, and the residual vector and the
   * stiffness matrix are computed in a non-blocked, node-wise layout. Example: [node_1_displacement,
   * node_1_temperature, node_2_displacement, node_2_temperature, ...], in which nodes designate the nodes of the kernel
   * functions attached to the particle.
   *
   * For the stiffness matrix, a column-major layout is used.
   */

  class MarmotParticle {

  public:
    MarmotParticle() = default;

    virtual ~MarmotParticle() = default;

    /// Assign all the properties of the particle. The properties are stored in a double array, which is passed to the
    /// function. The number of properties is passed as an argument. The properties need to be follow the
    /// order of the property names returned by the getPropertyNames() function.
    virtual void setProperties( const double* properties, int nProperties ) = 0;

    /// Assign a single property of the particle. The property is stored in a double array, which is passed to the
    /// function.
    virtual void setProperty( const std::string& propertyName, const double* property ) = 0;

    /// Get the names of all the valid properties of the particle.
    /// Also, the order of the properties is important, as it is used to assign the properties in the
    /// assignProperties function.
    virtual std::vector< std::string > getPropertyNames() const = 0;

    virtual void assignMeshfreeKernelFunctions(
      const std::vector< const MarmotMeshfreeKernelFunction* >& kernelFunctions ) = 0;

    virtual void assignStateVars( double* stateVars, int nStateVars ) = 0;

    virtual void initializeYourself() = 0;

    virtual void acceptStateAndPosition() = 0;

    virtual StateView getStateView( const std::string& stateName, int qp ) const = 0;

    virtual int getNumberOfRequiredStateVars() const = 0;

    virtual int getDimension() const = 0;

    virtual int getNumberOfVertices() const = 0;

    virtual double getVolumeUndeformed() const = 0;

    /// Get the shape of the material point (e.g., point, line, triangle, tetrahedron) in Ensight Gold notation
    virtual std::string getParticleShape() const = 0;

    /// Get the coordinates of the vertices of the particle (e.g., the nodes of a tetrahedron) for visualization
    virtual void getVisualizationVertexCoordinates( double* coordinates ) const = 0;

    /// Get the coordinates of the vertices of the particle (e.g., the nodes of a tetrahedron) for testing the coverage
    /// of the particle by shape functions.
    virtual void getVertexCoordinates( double* coordinates ) const = 0;

    /// Get surface coordinates of a face of the particle (e.g., the face of a tetrahedron) for boundary loads
    virtual void getFaceCoordinates( int faceID, double* coordinates ) const = 0;

    /// Get the coordinates of the center of the particle (e.g., the center of a tetrahedron)
    virtual void getCenterCoordinates( double* coordinates ) const = 0;

    /// Get the number of dofs per attached node (e.g., 3 for displacement, 1 for temperature, = 4 in total)
    /// The actual number of dofs results from the number of attached nodes and the number of dofs per node.
    virtual int getNBaseDof() const = 0;

    /// Get the fields on which the particle lives (e.g., displacement, temperature, ...). For a particle, the
    /// contribution to each attached node is equal.
    virtual const std::vector< std::string >& getFields() const = 0;

    virtual void computePhysicsKernels( const double* dQ,
                                        double*       fInt,
                                        double*       dFInt_ddQ,
                                        double        timeNew,
                                        double        dT ) = 0;

    virtual void computeBodyLoad( int           type,
                                  const double* load,
                                  double*       fExt,
                                  double*       dExt_dQ,
                                  double        timeNew,
                                  double        dT ) const = 0;

    virtual void computeDistributedLoad( int           type,
                                         int           boundaryFaceID,
                                         const double* load,
                                         double*       fExt,
                                         double*       dExt_dQ,
                                         double        timeNew,
                                         double        dT ) const = 0;

    virtual void updatePhysicsExplicit( const double* dQ, double timeNew, double dT ){};

    virtual void computePhysicsKernelsExplicit( double* fInt ){};

    virtual void computeBodyLoadExplicit( int           type,
                                          const double* load,
                                          double*       fExt,
                                          double        timeNew,
                                          double        dT ) const {};

    virtual void computeDistributedExplicit( int           type,
                                             int           boundaryFaceID,
                                             const double* load,
                                             double*       fExt,
                                             double        timeNew,
                                             double        dT ) const {};

    virtual void computeLumpedInertia( double* mLumped ) const { throw std::runtime_error( "Not implemented yet!" ); };

    virtual void computeLumpedMomentum( double* mLumped ) const { throw std::runtime_error( "Not implemented yet!" ); };

    virtual void getInterpolationVector( double* vec, const double* coordinates ) const = 0;

    virtual const std::unordered_map< std::string, int >& getSupportedBodyLoadTypes() const = 0;

    virtual const std::unordered_map< std::string, int >& getSupportedDistributedLoadTypes() const = 0;

    /// Get the coordinates of the vertices of the particle (e.g., the nodes of a tetrahedron) for testing the coverage
    /// of the particle by shape functions.
    virtual void getEvaluationCoordinates( double* coordinates ) const = 0;

    virtual int getNumberOfEvaluationPoints() const = 0;

    /**
     * We also implement the Variationally Consistent Integration (VCI)
     * according to the work of Chen, Hillman and Rueter (2013). The VCI method is used to modify the test functions by
     * adding correction terms.
     */

    /// Get the number of VCI constraints
    virtual int vci_getNumberOfConstraints() = 0;

    /// Compute the boundary integral of the test function for all constraints
    virtual void vci_compute_Test_P_BoundaryIntegral( double*       R_AiC_RowMajor,
                                                      const double* boundarySurfaceVector,
                                                      int           boundaryFaceID ) = 0;

    /// Compute the volume integral of the gradient of the test function for all constraints
    virtual void vci_compute_TestGradient_P_Integral( double* R_AiC_RowMajow ) = 0;

    /// Compute the kernel localization integral for all constraints
    virtual void vci_compute_Test_PGradient_Integral( double* R_AiC_RowMajor ) = 0;

    /// Compute the kernel localization integral for all constraints
    virtual void vci_compute_MMatrix( double* M_ACD_RowMajor ) = 0;

    /// Assign the correction terms to the test (shape) functions for all constraints
    virtual void vci_assignTestFunctionCorrectionTerms( const double* eta_AiC_RowMajor ) = 0;

    virtual void setInitialCondition( const std::string& conditionName, const double* value ) = 0;
  };

} // namespace Marmot::Meshfree
