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
#include "Marmot/MarmotMeshfreeApproximation.h"
#include "Marmot/MarmotMonomialBasisFunctions.h"
#include "Marmot/MarmotParticle.h"
#include <Eigen/Dense> // For Eigen::Matrix
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace Marmot::Meshfree {

  /**
   * @class Marmot::Meshfree::GenericParticle
   * @brief Common base of the point particles (one evaluation point at the particle center), independent of the
   * physics.
   *
   * GenericParticle implements the parts of the MarmotParticle interface that do not depend on the physics:
   *  - the geometry of a point particle: a single vertex at the center @f$ \boldsymbol{Y}_c @f$ in the reference
   *    intermediate configuration (the last accepted configuration), shape "point";
   *  - the shape functions: assignMeshfreeKernelFunctions() evaluates the trial functions @f$ N_A @f$ and their
   *    gradients @f$ \partial N_A / \partial Y_i @f$ at the center (direct nodal integration); the test functions
   *    @f$ T_A @f$, @f$ \partial T_A / \partial Y_i @f$ start as copies of them (Bubnov-Galerkin) and are changed
   *    only by the VCI correction (Petrov-Galerkin). Derived classes (e.g., the SQCNI particles) may replace the
   *    gradients by smoothed ones;
   *  - the property "VCI order" and the particle contributions to the variationally consistent integration (VCI)
   *    of Chen, Hillman and Rüter (2013).
   *
   * **VCI.** For the VCI order @f$ k @f$, the basis @f$ \boldsymbol{P}(\boldsymbol{Y}) @f$ is the complete monomial
   * basis of degree @f$ \le k @f$ in @f$ n_{\text{dim}} @f$ variables (Math::computeMonomialBasis()), with
   * @f$ n_C = \binom{k + n_{\text{dim}}}{n_{\text{dim}}} @f$ entries, evaluated at the particle center. With the
   * particle volume @f$ V_Y @f$ in the reference intermediate configuration as the integration weight, the particle
   * contributes
   * @f[
   *   \int \partial_i T_A\, P_C \,dV \approx \partial_i T_A\, P_C\, V_Y, \qquad
   *   \int T_A\, \partial_i P_C \,dV \approx T_A\, \partial_i P_C\, V_Y, \qquad
   *   M_{ACD} \approx \chi_A\, P_C\, P_D\, V_Y,
   * @f]
   * where @f$ \chi_A \in \{0, 1\} @f$ indicates whether the particle center is in the support of kernel function
   * @f$ A @f$. The correction assigned by vci_assignTestFunctionCorrectionTerms() is
   * @f$ \partial_i T_A \mathrel{+}= \chi_A \sum_C \eta_{AiC} P_C @f$. The boundary integral depends on the
   * kinematics of the physics and is left to the derived classes.
   *
   * Derived classes must keep @c _centerReferenceIntermediate, @c _volReferenceIntermediate, @c _P and
   * @c _P_Gradient up to date in acceptStateAndPosition().
   *
   * @tparam nDim The number of dimensions (2 or 3).
   */
  template < int nDim >
  class GenericParticle : public Marmot::Meshfree::MarmotParticle {

  protected:
    /// Center coordinates in the undeformed configuration.
    Eigen::Matrix< double, nDim, 1 > _centerCoordinatesUndeformed;
    /// Center coordinates @f$ \boldsymbol{Y}_c @f$ in the reference intermediate configuration.
    Eigen::Matrix< double, nDim, 1 > _centerReferenceIntermediate;
    /// Volume @f$ V_Y @f$ in the reference intermediate configuration (the integration weight of the VCI integrals).
    double _volReferenceIntermediate;

    /// The meshfree approximation (shape functions).
    const MarmotMeshfreeApproximation& _meshfreeApproximation;
    /// The assigned kernel functions (nodes), in dof order.
    std::vector< const MarmotMeshfreeKernelFunction* > _assignedKernelFunctions;

    /// The number of currently assigned nodes (= meshfree kernel functions)
    int _nNodes;

    int             _vciOrder;        ///< Order @f$ k @f$ of the VCI polynomial basis.
    int             _nVCIConstraints; ///< Number of VCI constraints @f$ n_C @f$ (size of the basis).
    Eigen::VectorXd _P;               ///< VCI basis @f$ P_C @f$ at the particle center, size @f$ n_C @f$.
    Eigen::MatrixXd _P_Gradient;      ///< Gradient @f$ \partial P_C / \partial Y_i @f$ at the particle center,
                                      ///< @f$ n_C \times n_{\text{dim}} @f$.

    /// The vector of trial shape functions @f$ N_A @f$ (@f$ 1 \times n_{\text{nodes}} @f$)
    Eigen::MatrixXd _N;
    /// The matrix of trial shape functions gradients @f$ \partial N_A / \partial Y_i @f$
    /// (@f$ n_{\text{dim}} \times n_{\text{nodes}} @f$)
    Eigen::MatrixXd _dN_dY;

    /// The vector of test shape functions @f$ T_A @f$ (@f$ 1 \times n_{\text{nodes}} @f$)
    Eigen::MatrixXd _T;
    /// The matrix of test shape functions gradients @f$ \partial T_A / \partial Y_i @f$
    /// (@f$ n_{\text{dim}} \times n_{\text{nodes}} @f$), including the VCI correction
    Eigen::MatrixXd _dT_dY;

    /// static vector of valid properties
    inline static const std::vector< std::string > _validProperties = {
      "VCI order",
    };

    /**
     * @brief Set the VCI order, resize the VCI basis and evaluate it at the current particle center.
     * @details The number of constraints is the number of monomials of degree @f$ \le @f$ @p order in nDim
     * variables (Math::computeSizeOfMonomialBasisVector()). The basis is evaluated right away, since VCI may run
     * before the first accepted increment updates it.
     * @param[in] order The VCI order @f$ k @f$.
     */
    void setVCIOrder( int order )
    {
      _vciOrder = order;
      // number of monomials of degree <= order in nDim variables
      _nVCIConstraints = Math::computeSizeOfMonomialBasisVector( order, nDim );
      _P.resize( _nVCIConstraints );
      _P_Gradient.resize( _nVCIConstraints, nDim );
      // evaluate the basis right away: VCI may run before the first accepted increment updates it
      Math::computeMonomialBasis( order, _centerReferenceIntermediate, _P );
      Math::computeMonomialBasisGradient( order, _centerReferenceIntermediate, _P_Gradient );
    };

  public:
    /**
     * @brief Construct a point particle.
     * @param[in] elementID The number of the particle (unused here).
     * @param[in] centerCoordinates0 Center coordinates in the undeformed configuration.
     * @param[in] nCenterCoordinates0 Number of values in @p centerCoordinates0; must be nDim.
     * @param[in] volume The undeformed volume, used as the initial reference intermediate volume.
     * @param[in] approximation The meshfree approximation; it must outlive the particle.
     * @throws std::invalid_argument if @p nCenterCoordinates0 != nDim.
     */
    GenericParticle( int                                elementID,
                     const double*                      centerCoordinates0,
                     int                                nCenterCoordinates0,
                     double                             volume,
                     const MarmotMeshfreeApproximation& approximation );

    /**
     * @brief Set the properties of GenericParticle, in the order of its own property list ("VCI order").
     * @param[in] properties Property values.
     * @param[in] nProperties Number of values; must be 1.
     * @throws std::runtime_error if the number of properties does not match.
     */
    virtual void setProperties( const double* properties, int nProperties ) override
    {
      if ( nProperties != static_cast< int >( _validProperties.size() ) ) {
        std::ostringstream oss;
        oss << "Error in " << __PRETTY_FUNCTION__ << ": ";
        oss << "Expected " << _validProperties.size() << " properties, but got " << nProperties << ". ";
        oss << "Valid properties are: ";
        for ( const auto& prop : _validProperties ) {
          oss << prop << ", ";
        }
        throw std::runtime_error( oss.str() );
      }

      for ( int i = 0; i < nProperties; i++ ) {
        setProperty( _validProperties[i], &properties[i] );
      }
    };

    /**
     * @brief Set a property by name; GenericParticle knows only "VCI order" (see setVCIOrder()).
     * @param[in] propertyName Name of the property.
     * @param[in] property Pointer to the value (the VCI order is truncated to an integer).
     * @throws std::runtime_error for an unknown property name.
     */
    virtual void setProperty( const std::string& propertyName, const double* property ) override
    {
      if ( propertyName == "VCI order" ) {
        _vciOrder = static_cast< int >( property[0] );
        this->setVCIOrder( _vciOrder );
      }
      else {
        std::ostringstream oss;
        oss << "Property " << propertyName << " not supported by GenericParticle! Valid properties are: ";
        for ( const auto& prop : _validProperties ) {
          oss << prop << ", ";
        }
        throw std::runtime_error( oss.str() );
      }
    };

    /// @brief Get the property names of GenericParticle.
    /// @return {"VCI order"}.
    virtual std::vector< std::string > getPropertyNames() const override { return _validProperties; };

    /**
     * @brief Store the kernel functions and evaluate @f$ N_A @f$ and @f$ \partial N_A / \partial Y_i @f$ at the
     * particle center (getCenterCoordinates()); the test functions are set equal to the trial functions, which
     * resets any previous VCI correction.
     * @param[in] kernelFunctions The kernel functions (nodes) covering the particle.
     */
    virtual void assignMeshfreeKernelFunctions(
      const std::vector< const MarmotMeshfreeKernelFunction* >& kernelFunctions ) override
    {
      _assignedKernelFunctions = kernelFunctions;

      _nNodes = _assignedKernelFunctions.size();

      Eigen::Matrix< double, nDim, 1 > coords;
      getCenterCoordinates( coords.data() );

      _N     = Eigen::MatrixXd::Zero( 1, _nNodes );
      _dN_dY = Eigen::MatrixXd::Zero( nDim, _nNodes );

      _meshfreeApproximation.computeShapeFunctionsAndGradients( coords.data(),
                                                                _assignedKernelFunctions,
                                                                _N.data(),
                                                                _dN_dY.data() );

      _T     = _N;
      _dT_dY = _dN_dY;
    }

    /// @brief Get the only vertex of a point particle, its center in the reference intermediate configuration.
    /// @param[out] coordinates Array of nDim values.
    virtual void getVertexCoordinates( double* coordinates ) const override
    {
      // Default implementation: particle center is its only vertex
      Eigen::Map< Eigen::Matrix< double, nDim, 1 > > coordinatesMap( coordinates );
      coordinatesMap = _centerReferenceIntermediate;
    }

    /// @brief A point particle has no faces.
    /// @param[in] faceID Unused.
    /// @param[out] coordinates Unused.
    /// @throws std::runtime_error always.
    virtual void getFaceCoordinates( int faceID, double* coordinates ) const override
    {
      throw std::runtime_error( "Error: GenericParticle::getFaceCoordinates not implemented." );
    }

    /// @brief Get the center in the reference intermediate configuration (same as getVertexCoordinates()).
    /// @param[out] coordinates Array of nDim values.
    virtual void getCenterCoordinates( double* coordinates ) const override { getVertexCoordinates( coordinates ); }

    /// @brief Get the visualization vertex, the particle center.
    /// @param[out] coordinates Array of nDim values.
    virtual void getVisualizationVertexCoordinates( double* coordinates ) const override
    {
      getVertexCoordinates( coordinates );
    };

    /// @brief A point particle has one vertex.
    /// @return 1.
    virtual int getNumberOfVertices() const override { return 1; };

    /// @brief Get the Ensight shape of a point particle.
    /// @return "point".
    virtual std::string getParticleShape() const override { return "point"; }

    /// @brief Get the spatial dimension.
    /// @return nDim.
    virtual int getDimension() const override { return nDim; };

    /// @brief Evaluate the shape functions of the assigned kernel functions at @p coordinates.
    /// @param[out] vec Array of @f$ n_{\text{nodes}} @f$ values.
    /// @param[in] coordinates Evaluation point (nDim values).
    virtual void getInterpolationVector( double* vec, const double* coordinates ) const override
    {
      _meshfreeApproximation.computeShapeFunctions( coordinates, _assignedKernelFunctions, vec );
    };

    /// @brief Get the only evaluation point, the particle center.
    /// @param[out] coordinates Array of nDim values.
    virtual void getEvaluationCoordinates( double* coordinates ) const override { getVertexCoordinates( coordinates ); }

    /// @brief A point particle has one evaluation point.
    /// @return 1.
    virtual int getNumberOfEvaluationPoints() const override
    {
      return 1; // only one evaluation point at the center of the particle
    };

    // VCI:
    /// @brief Get the number of VCI constraints.
    /// @return @f$ n_C @f$.
    virtual int vci_getNumberOfConstraints() override { return _nVCIConstraints; }

    /**
     * @brief Not available: the boundary integral needs the physics-specific kinematics (@f$ dY/dX @f$).
     * @param[in,out] R_AiC_RowMajor Unused.
     * @param[in] boundarySurfaceVector Unused.
     * @param[in] boundaryFaceID Unused.
     * @throws std::runtime_error always; derived classes must override it.
     */
    virtual void vci_compute_Test_P_BoundaryIntegral( double*       R_AiC_RowMajor,
                                                      const double* boundarySurfaceVector,
                                                      int           boundaryFaceID ) override
    {
      // This method depends on dY_dX, which is physics-specific.
      // It must be implemented in derived classes.
      throw std::runtime_error( "Error: GenericParticle::vci_compute_Test_P_BoundaryIntegral not implemented. "
                                "Must be implemented in a physics-specific derived class." );
    };

    /**
     * @brief Accumulate @f$ R_{AiC} \mathrel{+}= \partial_i T_A\, P_C\, V_Y @f$.
     * @param[in,out] R_AiC_RowMajor Row-major array @f$ n_{\text{nodes}} \times n_{\text{dim}} \times n_C @f$.
     */
    virtual void vci_compute_TestGradient_P_Integral( double* R_AiC_RowMajor ) override
    {
      // dimensions of R_AiC_RowMajor: _nNodes x nDim x _nVCIConstraints
      //
      for ( int A = 0; A < _nNodes; A++ )
        for ( int i = 0; i < nDim; i++ )
          for ( int C = 0; C < _nVCIConstraints; C++ )
            R_AiC_RowMajor[A * ( nDim * _nVCIConstraints ) + i * _nVCIConstraints + C] += _dT_dY( i, A ) * _P( C ) *
                                                                                          _volReferenceIntermediate;
    };

    /**
     * @brief Accumulate @f$ R_{AiC} \mathrel{+}= T_A\, \partial_i P_C\, V_Y @f$.
     * @param[in,out] R_AiC_RowMajor Row-major array @f$ n_{\text{nodes}} \times n_{\text{dim}} \times n_C @f$.
     */
    virtual void vci_compute_Test_PGradient_Integral( double* R_AiC_RowMajor ) override
    {
      // dimensions of R_AiC_RowMajor: _nNodes x nDim x _nVCIConstraints
      for ( int A = 0; A < _nNodes; A++ )
        for ( int i = 0; i < nDim; i++ )
          for ( int C = 0; C < _nVCIConstraints; C++ )
            R_AiC_RowMajor[A * ( nDim * _nVCIConstraints ) + i * _nVCIConstraints + C] += _T( A ) *
                                                                                          _P_Gradient( C, i ) *
                                                                                          _volReferenceIntermediate;
    };

    /**
     * @brief Accumulate @f$ M_{ACD} \mathrel{+}= \chi_A P_C P_D V_Y @f$, with @f$ \chi_A = 1 @f$ if the particle
     * center is in the support of kernel function @f$ A @f$, else 0.
     * @param[in,out] mMatrix_ACD_RowMajor Row-major array @f$ n_{\text{nodes}} \times n_C \times n_C @f$.
     */
    virtual void vci_compute_MMatrix( double* mMatrix_ACD_RowMajor ) override
    {
      // dimensions of R_AiC_RowMajor: _nNodes x nDim x _nVCIConstraints

      for ( int A = 0; A < _nNodes; A++ ) {
        const double R_A = _assignedKernelFunctions[A]->isInSupport( _centerReferenceIntermediate.data() ) ? 1.0 : 0.0;

        for ( int C = 0; C < _nVCIConstraints; C++ )
          for ( int D = 0; D < _nVCIConstraints; D++ )
            mMatrix_ACD_RowMajor[A * ( _nVCIConstraints * _nVCIConstraints ) + C * _nVCIConstraints +
                                 D] += R_A * _P( C ) * _P( D ) * _volReferenceIntermediate;
      }
    };

    /**
     * @brief Add the VCI correction to the test function gradients,
     * @f$ \partial_i T_A \mathrel{+}= \chi_A \sum_C \eta_{AiC} P_C @f$.
     * @details The correction is added to the current test function gradients, so it is meant to be applied once
     * after assignMeshfreeKernelFunctions() (which resets them).
     * @param[in] eta_AiC_RowMajor Row-major array @f$ n_{\text{nodes}} \times n_{\text{dim}} \times n_C @f$.
     */
    virtual void vci_assignTestFunctionCorrectionTerms( const double* eta_AiC_RowMajor ) override
    {

      for ( int A = 0; A < _nNodes; A++ ) {
        const double R_A = _assignedKernelFunctions[A]->isInSupport( _centerReferenceIntermediate.data() ) ? 1.0 : 0.0;
        for ( int i = 0; i < nDim; i++ ) {
          for ( int C = 0; C < _nVCIConstraints; C++ ) {
            _dT_dY( i, A ) += eta_AiC_RowMajor[A * ( nDim * _nVCIConstraints ) + i * _nVCIConstraints + C] * R_A *
                              _P( C );
          }
        }
      }
    };

    /**
     * @brief GenericParticle has no states of its own; derived classes provide them.
     * @param[in] stateName Name of the state.
     * @param[in] qp Unused.
     * @return Never returns.
     * @throws std::runtime_error always.
     */
    virtual StateView getStateView( const std::string& stateName, int qp ) const override
    {
      throw std::runtime_error( MakeString() << __PRETTY_FUNCTION__ << ": State " << stateName
                                             << " not supported by GenericParticle." );
    }
  };

  template < int nDim >
  GenericParticle< nDim >::GenericParticle( int                                                  elementID,
                                            const double*                                        centerCoordinates0,
                                            int                                                  nCenterCoordinates0,
                                            double                                               volume,
                                            const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
    : _centerCoordinatesUndeformed( Eigen::Map< const Eigen::Matrix< double, nDim, 1 > >( centerCoordinates0 ) ),
      _centerReferenceIntermediate( _centerCoordinatesUndeformed ),
      _volReferenceIntermediate( volume ), // Initial volume is undeformed volume
      _meshfreeApproximation( approximation ),
      _nNodes( 0 ),
      _vciOrder( 0 ), // Default VCI order
      _nVCIConstraints( 0 )
  {
    if ( nCenterCoordinates0 != nDim ) {
      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": size of center coordinates must be "
                                                << nDim << ", but got " << nCenterCoordinates0 );
    }
    this->setVCIOrder( _vciOrder ); // Initialize VCI members
  }

} // namespace Marmot::Meshfree
