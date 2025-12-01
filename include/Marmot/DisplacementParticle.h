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

#include "Marmot/DisplacementMaterialPoint.h"
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include "Marmot/MarmotMeshfreeApproximation.h"
#include "Marmot/MarmotMonomialBasisFunctions.h"
#include "Marmot/MarmotParticle.h"
#include "Marmot/MarmotUtils.h"
#include "Marmot/NewmarkBetaIntegrator.h"
#include <vector>

#include "Marmot/MarmotMeshfreeQuadHexCell.h"

namespace Marmot::Meshfree {

  template < int nDim >
  class DisplacementParticle : public Marmot::Meshfree::MarmotParticle {

  protected:
    Eigen::Matrix< double, nDim, 1 > _centerCoordinatesUndeformed;
    Eigen::Matrix< double, nDim, 1 > _centerReferenceIntermediate;
    double                           _volReferenceIntermediate;

    MaterialPoints::DisplacementMaterialPoint< nDim >* __mp;
    MaterialPoints::DisplacementMaterialPoint< nDim >& _mp;
    const MarmotMeshfreeApproximation&                 _meshfreeApproximation;
    std::vector< const MarmotMeshfreeKernelFunction* > _assignedKernelFunctions;

    double _newmark_beta;
    double _newmark_gamma;

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
      "newmark-beta beta",
      "newmark-beta gamma",
      "VCI order",
    };

  public:
    enum BodyLoadTypes {
      BodyForce,
    };

    enum DistributedLoadTypes { Pressure };

    const std::unordered_map< std::string, int >& getSupportedBodyLoadTypes() const override
    {
      static const std::unordered_map< std::string, int > _supportedBodyLoadTypes = { { "BODYFORCE", BodyForce } };
      return _supportedBodyLoadTypes;
    };

    const std::unordered_map< std::string, int >& getSupportedDistributedLoadTypes() const override
    {
      static const std::unordered_map< std::string, int > _supportedDistributedLoadTypes = { { "PRESSURE", Pressure } };
      return _supportedDistributedLoadTypes;
    };

    static constexpr int nDofPerNodeU = nDim; // Displacement   field U

    using Material = MarmotMaterialFiniteStrain;

    using ForceSized = Eigen::Matrix< double, nDim, 1 >;

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

    virtual void setProperty( const std::string& propertyName, const double* property )
    {
      if ( propertyName == "newmark-beta beta" ) {
        _newmark_beta = property[0];
      }
      else if ( propertyName == "newmark-beta gamma" ) {
        _newmark_gamma = property[0];
      }
      else if ( propertyName == "VCI order" ) {
        _vciOrder = static_cast< int >( property[0] );
        this->setVCIOrder( _vciOrder );
      }
      else {
        std::ostringstream oss;
        oss << "Property " << propertyName << " not supported! Valid properties are: ";
        for ( const auto& prop : _validProperties ) {
          oss << prop << ", ";
        }
        throw std::runtime_error( oss.str() );
      }
    };

    virtual std::vector< std::string > getPropertyNames() const { return _validProperties; };

    virtual int getNumberOfRequiredStateVars() const override { return _mp.getNumberOfRequiredStateVars(); };

    void assignStateVars( double* stateVars, int nStateVars ) override { _mp.assignStateVars( stateVars, nStateVars ); }

    void assignMeshfreeKernelFunctions(
      const std::vector< const MarmotMeshfreeKernelFunction* >& kernelFunctions ) override
    {
      _assignedKernelFunctions = kernelFunctions;

      _nNodes = _assignedKernelFunctions.size();

      Eigen::Matrix< double, nDim, 1 > coords;
      _mp.getCoordinatesAtCenter( coords.data() );

      _N     = Eigen::MatrixXd::Zero( 1, _nNodes );
      _dN_dY = Eigen::MatrixXd::Zero( nDim, _nNodes );

      _meshfreeApproximation.computeShapeFunctionsAndGradients( coords.data(),
                                                                _assignedKernelFunctions,
                                                                _N.data(),
                                                                _dN_dY.data() );

      _T     = _N;
      _dT_dY = _dN_dY;
    }

    virtual int getNBaseDof() const { return nDofPerNodeU; }

    virtual const std::vector< std::string >& getFields() const override
    {
      static const std::vector< std::string > nodeFields = { "displacement" };
      return nodeFields;
    };

    virtual void getVertexCoordinates( double* coordinates ) const override { _mp.getVertexCoordinates( coordinates ); }

    virtual void getFaceCoordinates( int faceID, double* coordinates ) const override
    {
      throw std::runtime_error( "Error: DisplacementParticle::getFaceCoordinates not implemented." );
    }

    virtual void getCenterCoordinates( double* coordinates ) const override { getVertexCoordinates( coordinates ); }

    virtual void getVisualizationVertexCoordinates( double* coordinates ) const override
    {
      getVertexCoordinates( coordinates );
    };

    virtual int getNumberOfVertices() const override { return 1; };

    virtual std::string getParticleShape() const override { return "point"; }

    virtual int getDimension() const override { return nDim; };

    DisplacementParticle( int                                elementID,
                          const double*                      nodeCoordinates,
                          int                                nNodeCoordiantes,
                          double                             volume,
                          const std::string&                 materialName,
                          const double*                      materialProperties,
                          int                                nMaterialProperties,
                          const MarmotMeshfreeApproximation& approximation );

    void initializeYourself() override { _mp.initializeYourself(); };

    virtual void acceptStateAndPosition() override
    {

      _mp.acceptStateAndPosition();
      _mp.prepareYourself( 0, 0 );

      updateVolumeToReferenceIntermediate();
      updateParticlePositionToReferenceIntermediate();

      Math::computeMonomialBasis( _vciOrder, _centerReferenceIntermediate, _P );
      Math::computeMonomialBasisGradient( _vciOrder, _centerReferenceIntermediate, _P_Gradient );
    };

    virtual void computePhysicsKernels( const double* dQ,
                                        double*       fInt,
                                        double*       dFInt_ddQ,
                                        double        timeNew,
                                        double        dT ) override;

    virtual void computeBodyLoad( int           type,
                                  const double* load,
                                  double*       fExt,
                                  double*       dExt_dQ,
                                  double        timeNew,
                                  double        dT ) const override;

    virtual void computeDistributedLoad( int           type,
                                         int           surfaceID,
                                         const double* load,
                                         double*       fExt,
                                         double*       dExt_dQ,
                                         double        timeNew,
                                         double        dT ) const override;

    virtual StateView getStateView( const std::string& stateName, int qp ) const override;

    virtual void getInterpolationVector( double* vec, const double* coordinates ) const override
    {
      _meshfreeApproximation.computeShapeFunctions( coordinates, _assignedKernelFunctions, vec );
    };

    // VCI:

    virtual int vci_getNumberOfConstraints() override { return _nVCIConstraints; }

    virtual void vci_compute_Test_P_BoundaryIntegral( double*       R_AiC_RowMajor,
                                                      const double* boundarySurfaceVector,
                                                      int           boundaryFaceID ) override
    {

      using namespace Fastor;

      Tensor< double, nDim > n_dA0( boundarySurfaceVector ); // undeformed load vector p * N_I * dA_0

      // apply Nanson's formula
      const Tensor< double, nDim, nDim > FInv = inverse( _mp.dY_dX() );
      const double                       J    = determinant( _mp.dY_dX() );

      const Tensor< double, nDim > n_dAY = J * transpose( FInv ) % n_dA0;

      for ( int A = 0; A < _nNodes; A++ )
        for ( int i = 0; i < nDim; i++ )
          for ( int C = 0; C < _nVCIConstraints; C++ )
            R_AiC_RowMajor[A * ( nDim * _nVCIConstraints ) + i * _nVCIConstraints + C] += _T( A ) * _P( C ) * n_dAY[i];
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

    virtual double getVolumeUndeformed() const { return _mp.getVolumeUndeformed(); };

    virtual void getEvaluationCoordinates( double* coordinates ) const { getVertexCoordinates( coordinates ); }

    virtual int getNumberOfEvaluationPoints() const
    {
      return 1; // only one evaluation point at the center of the particle
    };

    virtual void setInitialCondition( const std::string& conditionName, const double* value ) override
    {
      if ( conditionName == "geostaticstress" ) {
        _mp.setInitialCondition( conditionName, value );
      }
      else {
        throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": invalid initial condition" );
      }
    };

  private:
    virtual void updateParticlePositionToReferenceIntermediate()
    {
      _mp.getVertexCoordinates( _centerReferenceIntermediate.data() );
    };

    virtual void updateVolumeToReferenceIntermediate()
    {
      _volReferenceIntermediate = _mp.getVolumeUndeformed() * determinant( _mp.dY_dX() );
    };

    void setVCIOrder( int order )
    {
      _vciOrder        = order;
      _nVCIConstraints = ( order + 1 ) * ( order + 2 ) / 2; // number of VCI constraints for polynomial basis of order
      _P.resize( _nVCIConstraints );
      _P_Gradient.resize( _nVCIConstraints, nDim );
    };
  };

  template < int nDim >
  StateView DisplacementParticle< nDim >::getStateView( const std::string& stateName, int qp ) const
  {
    return _mp.getStateView( stateName );
  }

  template < int nDim >
  DisplacementParticle< nDim >::DisplacementParticle(
    int                                                  elementID,
    const double*                                        centerCoordinates0,
    int                                                  sizeCenterCoordinates0,
    double                                               volume,
    const std::string&                                   materialName,
    const double*                                        materialProperties,
    int                                                  nMaterialProperties,
    const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
    : _centerCoordinatesUndeformed( Eigen::Map< const Eigen::Matrix< double, nDim, 1 > >( centerCoordinates0 ) ),
      _centerReferenceIntermediate( _centerCoordinatesUndeformed ),
      __mp( new MaterialPoints::DisplacementMaterialPoint2D( elementID,
                                                             _centerCoordinatesUndeformed.data(),
                                                             _centerCoordinatesUndeformed.size(),
                                                             volume ) ),
      _mp( *__mp ),
      _meshfreeApproximation( approximation ),
      _newmark_beta( 0. ),
      _newmark_gamma( 0. ),
      _vciOrder( 0 )
  {
    if ( sizeCenterCoordinates0 != nDim ) {
      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": size of center coordinates must be "
                                                << nDim << ", but got " << sizeCenterCoordinates0 );
    }

    int                   materialCode = MarmotLibrary::MarmotMaterialFactory::getMaterialCodeFromName( materialName );
    MarmotMaterialSection section( materialCode, materialProperties, nMaterialProperties );

    _mp.assignMaterial( section );

    this->setVCIOrder( _vciOrder );
  }

  template < int nDim >
  void DisplacementParticle< nDim >::computePhysicsKernels( const double* dQ,
                                                            double*       fInt,
                                                            double*       dFInt_ddQ,
                                                            double        timeNew,
                                                            double        dT )
  {

    using namespace Marmot::FastorIndices;
    using namespace Fastor;
    using to_jk = Fastor::OIndex< j_, k_ >;

    const static Tensor< double, nDim, nDim > I(
      ( Eigen::Matrix< double, nDim, nDim >() << Eigen::Matrix< double, nDim, nDim >::Identity() ).finished().data() );

    constexpr int nodeBlockSize = nDim;

    Tensor< double, nDim > du( 0.0 );

    Tensor< double, nDim, nDim > du_dY( 0.0 );

    for ( int B = 0; B < _nNodes; B++ ) {

      const int idxB_u = nodeBlockSize * B;

      const double N_B     = _N( B );
      const auto   dN_B_dY = Tensor< double, nDim >( _dN_dY.col( B ).data() ); // works because ColumnMajor of Eigen

      const auto dQU = Tensor< double, nDim >( dQ + idxB_u );

      du += N_B * dQU;

      du_dY += einsum< i, j >( dQU, dN_B_dY );
    }

    _mp.prepareYourself( timeNew, dT );
    _mp.incrementDeformation( du, du_dY );
    _mp.computeYourself( timeNew, dT );

    const double density0 = _mp.getDensityUndeformed();

    auto v = _mp.getVelocity();
    auto a = _mp.getAcceleration();

    Tensor< double, nDim, nDim > da_ddu( 0.0 );
    Marmot::TimeIntegration::newmarkBetaIntegration< nDim >( du.data(),
                                                             v.data(),
                                                             a.data(),
                                                             dT,
                                                             this->_newmark_beta,
                                                             this->_newmark_gamma,
                                                             da_ddu.data() );
    _mp.setVelocity( v );
    _mp.setAcceleration( a );

    Tensor< double, nDim > r_U( 0.0 );

    Tensor< double, nDim, nDim > k_UU( 0.0 );

    const auto& S = _mp.response.S;

    const double V0 = getVolumeUndeformed();

    const auto& t = _mp.tangents;

    Eigen::Map< Eigen::VectorXd > P( fInt, _nNodes * nodeBlockSize );
    Eigen::Map< Eigen::MatrixXd > K( dFInt_ddQ, _nNodes * nodeBlockSize, _nNodes * nodeBlockSize );

    // clang-format off
    for ( int A = 0; A < _nNodes; A++ ) {

      const double T_A = _T( A );
      const auto                   dT_A_dY = TensorMap< const double, nDim >( _dT_dY.col( A ).data() );
      const Tensor< double, nDim > dT_A_dx = einsum< ji, j >( inv( _mp.dx_dY() ), dT_A_dY );
      const Tensor< double, nDim > dT_A_dX = einsum< ji, j >( _mp.dY_dX(), dT_A_dY );


        const int idxA_u = nodeBlockSize * A;

        r_U = ( +einsum< i, ij >( dT_A_dx, S ) ) * V0;

        // add inertia
        r_U += density0 * a * T_A * V0;

        {
          using namespace Eigen;
          P.template segment< nDim >( idxA_u ) += Map< Matrix< double, nDim, 1 > >( r_U.data() );
        }

        for ( int B = 0; B < _nNodes; B++ ) {

          const int idxB_u = nodeBlockSize * B;

          const double                 N_B     = _N( B );
          const auto dN_B_dY = TensorMap< const double, nDim >( _dN_dY.col(B).data() );
          const auto dN_B_dx = evaluate( einsum< ji, j >( inv( _mp.dx_dY() ), dN_B_dY ) );
          const auto dN_B_dX = evaluate( einsum< ji, j >( _mp.dY_dX(), dN_B_dY ) ); // no dependence on current deformations!

          // aux stiffness tensors
          const auto dS_dqU_B = evaluate ( + einsum < ijkl, l > ( t.dS_dDeltaF, dN_B_dY )                                            );

          k_UU  = ( + einsum< i, ijk        > ( dT_A_dx, dS_dqU_B )                                                       ) * V0;

          k_UU += ( - einsum< k, ij, i, to_jk >( dT_A_dx, S, dN_B_dx ) ) * V0;

          k_UU += density0 * da_ddu * T_A * N_B * V0;

          {
              using namespace Eigen;
              // TODO: check if we can use transpose instead of torowmajor:
              K.template block< nDim, nDim >( idxA_u, idxB_u ) += Map< Matrix< double, nDim, nDim > >( torowmajor( k_UU ).data() );
          }
      }
    }
    // clang-format on
  }

  template < int nDim >
  void DisplacementParticle< nDim >::computeDistributedLoad( int           type,
                                                             int           surfaceID,
                                                             const double* load,
                                                             double*       fExt,
                                                             double*       dExt_dQ,
                                                             double        timeNew,
                                                             double        dT ) const
  {
  }

  template < int nDim >
  void DisplacementParticle< nDim >::computeBodyLoad( int           type,
                                                      const double* load,
                                                      double*       fExt,
                                                      double*       dExt_dQ,
                                                      double        timeNew,
                                                      double        dT ) const
  {
  }

} // namespace Marmot::Meshfree
