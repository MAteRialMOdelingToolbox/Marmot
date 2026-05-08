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
#include "Marmot/GenericParticle.h" // New base class
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include "Marmot/MarmotMeshfreeApproximation.h"
#include "Marmot/MarmotMonomialBasisFunctions.h"
#include "Marmot/MarmotParticle.h"
#include "Marmot/MarmotUtils.h"
#include "Marmot/NewmarkBetaIntegrator.h"
#include <vector>

template < int nDim >
using MaterialPointType = std::conditional_t< nDim == 2,
                                              Marmot::MaterialPoints::DisplacementMaterialPoint2D,
                                              Marmot::MaterialPoints::DisplacementMaterialPoint3D >;

namespace Marmot::Meshfree {

  template < int nDim >
  class DisplacementParticle : public Marmot::Meshfree::GenericParticle< nDim > {

    using TensorD  = Fastor::Tensor< double, nDim >;
    using TensorDD = Fastor::Tensor< double, nDim, nDim >;

  protected:
    std::unique_ptr< MaterialPointType< nDim > > _mp;

    double _newmark_beta;
    double _newmark_gamma;

    /// static vector of valid properties
    inline static const std::vector< std::string > _validProperties = {
      "newmark-beta beta",
      "newmark-beta gamma",
    };

    virtual TensorDD dY_dX() const { return ( this->_mp->dY_dX() ); }

    virtual TensorDD dx_dY() const { return ( this->_mp->dx_dY() ); }

    virtual TensorD getDisplacementAtCenter() const
    {
      TensorD u( 0.0 );
      this->_mp->getCenterDisplacement( u.data() );
      return u;
    }

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
      // Combine property names from base and derived
      std::vector< std::string > allPropertyNames = Marmot::Meshfree::GenericParticle< nDim >::getPropertyNames();
      allPropertyNames.insert( allPropertyNames.end(), _validProperties.begin(), _validProperties.end() );

      if ( nProperties != static_cast< int >( allPropertyNames.size() ) ) {
        std::ostringstream oss;
        oss << "Error in " << __PRETTY_FUNCTION__ << ": ";
        oss << "Expected " << allPropertyNames.size() << " properties, but got " << nProperties << ". ";
        oss << "Valid properties are: ";
        for ( const auto& prop : allPropertyNames ) {
          oss << prop << ", ";
        }
        throw std::runtime_error( oss.str() );
      }

      for ( int i = 0; i < nProperties; i++ ) {
        setProperty( allPropertyNames[i], &properties[i] );
      }
    };

    virtual void setProperty( const std::string& propertyName, const double* property ) override
    {
      if ( propertyName == "newmark-beta beta" ) {
        _newmark_beta = property[0];
      }
      else if ( propertyName == "newmark-beta gamma" ) {
        _newmark_gamma = property[0];
      }
      else {
        // If not a DisplacementParticle specific property, try the base class
        Marmot::Meshfree::GenericParticle< nDim >::setProperty( propertyName, property );
      }
    };

    virtual std::vector< std::string > getPropertyNames() const override
    {
      std::vector< std::string > names = Marmot::Meshfree::GenericParticle< nDim >::getPropertyNames();
      names.insert( names.end(), _validProperties.begin(), _validProperties.end() );
      return names;
    };

    virtual int getNumberOfRequiredStateVars() const override { return _mp->getNumberOfRequiredStateVars(); };

    void assignStateVars( double* stateVars, int nStateVars ) override
    {
      _mp->assignStateVars( stateVars, nStateVars );
    }

    virtual int getNBaseDof() const { return nDofPerNodeU; }

    virtual const std::vector< std::string >& getFields() const override
    {
      static const std::vector< std::string > nodeFields = { "displacement" };
      return nodeFields;
    };

    DisplacementParticle( int                                elementID,
                          const double*                      nodeCoordinates,
                          int                                nNodeCoordiantes,
                          double                             volume,
                          const std::string&                 materialName,
                          const double*                      materialProperties,
                          int                                nMaterialProperties,
                          const MarmotMeshfreeApproximation& approximation );

    void initializeYourself() override { _mp->initializeYourself(); };

    virtual void acceptStateAndPosition() override
    {
      _mp->acceptStateAndPosition();
      _mp->prepareYourself( 0, 0 );

      this->updateVolumeToReferenceIntermediate();
      this->updateParticlePositionToReferenceIntermediate();

      // Use base class members for VCI
      Math::computeMonomialBasis( this->_vciOrder, this->_centerReferenceIntermediate, this->_P );
      Math::computeMonomialBasisGradient( this->_vciOrder, this->_centerReferenceIntermediate, this->_P_Gradient );
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

    virtual void computeLumpedInertia( double* mLumped ) const override;

    virtual void computeLumpedMomentum( double* mLumped ) const override;

    virtual StateView getStateView( const std::string& stateName, int qp ) const override;

    // VCI methods are now in GenericParticle, but vci_compute_Test_P_BoundaryIntegral needs override
    // because it depends on dY_dX() which is physics-specific.
    virtual void vci_compute_Test_P_BoundaryIntegral( double*       R_AiC_RowMajor,
                                                      const double* boundarySurfaceVector,
                                                      int           boundaryFaceID ) override
    {
      using namespace Fastor;

      Tensor< double, nDim > n_dA0( boundarySurfaceVector ); // undeformed load vector p * N_I * dA_0

      // apply Nanson's formula
      const Tensor< double, nDim, nDim > FInv = inverse( dY_dX() );
      const double                       J    = determinant( dY_dX() );

      const Tensor< double, nDim > n_dAY = J * transpose( FInv ) % n_dA0;

      for ( int A = 0; A < this->_nNodes; A++ )                            // Use base class _nNodes
        for ( int i = 0; i < nDim; i++ )
          for ( int C = 0; C < this->_nVCIConstraints; C++ )               // Use base class _nVCIConstraints
            R_AiC_RowMajor[A * ( nDim * this->_nVCIConstraints ) + i * this->_nVCIConstraints +
                           C] += this->_T( A ) * this->_P( C ) * n_dAY[i]; // Use base class _T, _P
    };

    virtual double getVolumeUndeformed() const { return _mp->getVolumeUndeformed(); };

    virtual void setInitialCondition( const std::string& conditionName, const double* value ) override
    {
      if ( conditionName == "geostaticstress" ) {
        _mp->setInitialCondition( conditionName, value );
      }
      else {
        throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": invalid initial condition" );
      }
    };

  private:
    virtual void updateParticlePositionToReferenceIntermediate()
    {
      _mp->getVertexCoordinates(
        this->_centerReferenceIntermediate.data() ); // Use base class _centerReferenceIntermediate
    };

    virtual void updateVolumeToReferenceIntermediate()
    {
      this->_volReferenceIntermediate = getVolumeUndeformed() *
                                        determinant( dY_dX() ); // Use base class _volReferenceIntermediate
    };
  };

  template < int nDim >
  StateView DisplacementParticle< nDim >::getStateView( const std::string& stateName, int qp ) const
  {
    return _mp->getStateView( stateName );
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
    : Marmot::Meshfree::GenericParticle< nDim >( elementID,
                                                 centerCoordinates0,
                                                 sizeCenterCoordinates0,
                                                 volume,
                                                 approximation ), // Call virtual base constructor
      _mp( std::make_unique< MaterialPointType< nDim > >( elementID,
                                                          Eigen::Map< const Eigen::Matrix< double, nDim, 1 > >(
                                                            centerCoordinates0 )
                                                            .data(),
                                                          sizeCenterCoordinates0,
                                                          volume ) ),
      _newmark_beta( 0. ),
      _newmark_gamma( 0. )
  {
    MarmotMaterialSection section( materialName, materialProperties, nMaterialProperties );

    _mp->assignMaterial( section );
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

    for ( int B = 0; B < this->_nNodes; B++ ) { // Use base class _nNodes

      const int idxB_u = nodeBlockSize * B;

      const double N_B     = this->_N( B );                                          // Use base class _N
      const auto   dN_B_dY = Tensor< double, nDim >( this->_dN_dY.col( B ).data() ); // Use base class _dN_dY

      const auto dQU = Tensor< double, nDim >( dQ + idxB_u );

      du += N_B * dQU;

      du_dY += einsum< i, j >( dQU, dN_B_dY );
    }

    _mp->prepareYourself( timeNew, dT );
    _mp->incrementDeformation( du, du_dY );
    _mp->computeYourself( timeNew, dT );

    const double density0 = _mp->getDensityUndeformed();

    auto v = _mp->getVelocity();
    auto a = _mp->getAcceleration();

    Tensor< double, nDim, nDim > da_ddu( 0.0 );
    Marmot::TimeIntegration::newmarkBetaIntegration< nDim >( du.data(),
                                                             v.data(),
                                                             a.data(),
                                                             dT,
                                                             this->_newmark_beta,
                                                             this->_newmark_gamma,
                                                             da_ddu.data() );
    _mp->setVelocity( v );
    _mp->setAcceleration( a );

    Tensor< double, nDim > r_U( 0.0 );

    Tensor< double, nDim, nDim > k_UU( 0.0 );

    const auto& S = _mp->response.S;

    const double V0 = getVolumeUndeformed();

    const auto& t = _mp->tangents;

    Eigen::Map< Eigen::VectorXd > P( fInt, this->_nNodes * nodeBlockSize ); // Use base class _nNodes
    Eigen::Map< Eigen::MatrixXd > K( dFInt_ddQ,
                                     this->_nNodes * nodeBlockSize,
                                     this->_nNodes * nodeBlockSize ); // Use base class _nNodes

    // clang-format off
    for ( int A = 0; A < this->_nNodes; A++ ) { // Use base class _nNodes

      const double T_A = this->_T( A ); // Use base class _T
      const auto                   dT_A_dY = TensorMap< const double, nDim >( this->_dT_dY.col( A ).data() ); // Use base class _dT_dY
      const Tensor< double, nDim > dT_A_dx = einsum< ji, j >( inv( _mp->dx_dY() ), dT_A_dY );

        const int idxA_u = nodeBlockSize * A;

        r_U = ( +einsum< i, ij >( dT_A_dx, S ) ) * V0;

        // add inertia
        r_U += density0 * a * T_A * V0;

        {
          using namespace Eigen;
          P.template segment< nDim >( idxA_u ) += Map< Matrix< double, nDim, 1 > >( r_U.data() );
        }

        for ( int B = 0; B < this->_nNodes; B++ ) { // Use base class _nNodes

          const int idxB_u = nodeBlockSize * B;

          const double                 N_B     = this->_N( B ); // Use base class _N
          const auto dN_B_dY = TensorMap< const double, nDim >( this->_dN_dY.col(B).data() ); // Use base class _dN_dY
          const auto dN_B_dx = evaluate( einsum< ji, j >( inv( _mp->dx_dY() ), dN_B_dY ) );

          // aux stiffness tensors
          const auto dS_dqU_B = evaluate ( + einsum < ijkl, l > ( t.dS_dDeltaF, dN_B_dY )                                            );

          k_UU  = ( + einsum< i, ijk        >  ( dT_A_dx, dS_dqU_B )   ) * V0;
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

  template < int nDim >
  void DisplacementParticle< nDim >::computeLumpedInertia( double* mLumped ) const
  {
    const double density0 = _mp->getDensityUndeformed();
    const double V0       = getVolumeUndeformed();

    for ( int A = 0; A < this->_nNodes; A++ ) { // Use base class _nNodes
      const double T_A    = this->_T( A );      // Use base class _T
      const int    idxA_u = nDofPerNodeU * A;

      for ( int i = 0; i < nDofPerNodeU; i++ ) {
        mLumped[idxA_u + i] += density0 * T_A * V0;
      }
    }
  }

  template < int nDim >
  void DisplacementParticle< nDim >::computeLumpedMomentum( double* mLumped ) const
  {
    const double density0 = _mp->getDensityUndeformed();
    const double V0       = getVolumeUndeformed();
    const auto   v        = _mp->getVelocity();

    for ( int A = 0; A < this->_nNodes; A++ ) { // Use base class _nNodes
      const double T_A    = this->_T( A );      // Use base class _T
      const int    idxA_u = nDofPerNodeU * A;

      for ( int i = 0; i < nDofPerNodeU; i++ ) {
        mLumped[idxA_u + i] += density0 * T_A * V0 * v[i];
      }
    }
  }

} // namespace Marmot::Meshfree
