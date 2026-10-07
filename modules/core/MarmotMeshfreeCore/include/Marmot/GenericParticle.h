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

#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMeshfreeApproximation.h"
#include "Marmot/MarmotParticle.h"
#include <Eigen/Dense> // For Eigen::Matrix
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace Marmot::Meshfree {

  template < int nDim >
  class GenericParticle : public Marmot::Meshfree::MarmotParticle {

  protected:
    Eigen::Matrix< double, nDim, 1 > _centerCoordinatesUndeformed;
    Eigen::Matrix< double, nDim, 1 > _centerReferenceIntermediate;
    double                           _volReferenceIntermediate;

    const MarmotMeshfreeApproximation&                 _meshfreeApproximation;
    std::vector< const MarmotMeshfreeKernelFunction* > _assignedKernelFunctions;

    /// The number of currently assigned nodes (=meshfree kernel functions}
    int _nNodes;

    int             _vciOrder; // order of the VCI polynomial basis
    int             _nVCIConstraints;
    Eigen::VectorXd _P;
    Eigen::MatrixXd _P_Gradient;

    /// The vector of trial shape functions
    Eigen::MatrixXd _N;
    /// The matrix of trial shape functions gradients
    Eigen::MatrixXd _dN_dY;

    /// The vector of test shape functions
    Eigen::MatrixXd _T;
    /// The matrix of test shape functions gradients
    Eigen::MatrixXd _dT_dY;

    /// static vector of valid properties
    inline static const std::vector< std::string > _validProperties = {
      "VCI order",
    };

    // Helper to set VCI order and resize related members
    void setVCIOrder( int order )
    {
      _vciOrder        = order;
      _nVCIConstraints = ( order + 1 ) * ( order + 2 ) / 2; // number of VCI constraints for polynomial basis of order
      _P.resize( _nVCIConstraints );
      _P_Gradient.resize( _nVCIConstraints, nDim );
    };

  public:
    GenericParticle( int                                elementID,
                     const double*                      centerCoordinates0,
                     int                                nCenterCoordinates0,
                     double                             volume,
                     const MarmotMeshfreeApproximation& approximation );

    // MarmotParticle interface overrides
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

    virtual std::vector< std::string > getPropertyNames() const override { return _validProperties; };

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

    virtual void getVertexCoordinates( double* coordinates ) const override
    {
      // Default implementation: particle center is its only vertex
      Eigen::Map< Eigen::Matrix< double, nDim, 1 > > coordinatesMap( coordinates );
      coordinatesMap = _centerReferenceIntermediate;
    }

    virtual void getFaceCoordinates( int faceID, double* coordinates ) const override
    {
      throw std::runtime_error( "Error: GenericParticle::getFaceCoordinates not implemented." );
    }

    virtual void getCenterCoordinates( double* coordinates ) const override { getVertexCoordinates( coordinates ); }

    virtual void getVisualizationVertexCoordinates( double* coordinates ) const override
    {
      getVertexCoordinates( coordinates );
    };

    virtual int getNumberOfVertices() const override { return 1; };

    virtual std::string getParticleShape() const override { return "point"; }

    virtual int getDimension() const override { return nDim; };

    virtual void getInterpolationVector( double* vec, const double* coordinates ) const override
    {
      _meshfreeApproximation.computeShapeFunctions( coordinates, _assignedKernelFunctions, vec );
    };

    virtual void getEvaluationCoordinates( double* coordinates ) const override { getVertexCoordinates( coordinates ); }

    virtual int getNumberOfEvaluationPoints() const override
    {
      return 1; // only one evaluation point at the center of the particle
    };

    // VCI:
    virtual int vci_getNumberOfConstraints() override { return _nVCIConstraints; }

    virtual void vci_compute_Test_P_BoundaryIntegral( double*       R_AiC_RowMajor,
                                                      const double* boundarySurfaceVector,
                                                      int           boundaryFaceID ) override
    {
      // This method depends on dY_dX, which is physics-specific.
      // It must be implemented in derived classes.
      throw std::runtime_error( "Error: GenericParticle::vci_compute_Test_P_BoundaryIntegral not implemented. "
                                "Must be implemented in a physics-specific derived class." );
    };

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
