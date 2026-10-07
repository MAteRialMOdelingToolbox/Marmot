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

#include "Marmot/GradientEnhancedFiniteStrainParticleSQCNI.h"
#include "Marmot/MarmotMeshfreeApproximation.h"
#include "Marmot/NewmarkBetaIntegrator.h"
#include <Fastor/Fastor.h>

namespace Marmot::Meshfree {

  /**
   * @brief SQCNI with NATURALLY STABILIZED nodal integration (NSNI), non-micropolar.
   *
   * The non-micropolar port of @ref GradientEnhancedMicropolarParticleSQCNIxNSNI.  On top of
   * the smoothed gradients of @ref GradientEnhancedFiniteStrainParticleSQCNI it adds the NSNI
   * stabilization term, built from the SECOND derivatives of the shape functions (also obtained
   * by boundary integration over the smoothing domain) contracted with the second moments of
   * the particle domain about its centroid.  That term is what removes the spurious zero-energy
   * modes of nodal integration -- which matters precisely on a fine mesh in softening, where an
   * unstabilised nodal scheme can produce oscillations that look like localisation.
   *
   * The micropolar version's `stabilize angular momentum` option and every couple-stress /
   * Levi-Civita term are gone: there is no micro-rotation field here to stabilise.
   */
  template < int nDim, int nVertices >
  class GradientEnhancedFiniteStrainParticleSQCNIxNSNI
    : public GradientEnhancedFiniteStrainParticleSQCNI< nDim, nVertices > {

    using ParentPointParticle = GradientEnhancedFiniteStrainParticle< nDim >;
    using ParentSQCNIParticle = GradientEnhancedFiniteStrainParticleSQCNI< nDim, nVertices >;
    using LagrangeCellType    = ParentSQCNIParticle::LagrangeCellType;

    using TensorD    = Fastor::Tensor< double, nDim >;
    using TensorDD   = Fastor::Tensor< double, nDim, nDim >;
    using TensorDDD  = Fastor::Tensor< double, nDim, nDim, nDim >;
    using TensorDDDD = Fastor::Tensor< double, nDim, nDim, nDim, nDim >;

    TensorDD                            _momentsOfInertia_Undeformed;
    TensorDD                            _momentsOfInertia_IntermediateReference;
    std::array< Eigen::MatrixXd, nDim > _d2N_dYdY;

  public:
    GradientEnhancedFiniteStrainParticleSQCNIxNSNI(
      int                                            elementID,
      const double*                                  nodeCoordinates,
      int                                            nNodeCoordiantes,
      double                                         volume,
      const std::string&                             materialName,
      const double*                                  materialProperties,
      int                                            sizeMaterialProperties,
      const MarmotMeshfreeApproximation&             approximation,
      ParentSQCNIParticle::SmoothingDomainUpdateType smoothingVolumeUpdateType );

    void assignMeshfreeKernelFunctions(
      const std::vector< const MarmotMeshfreeKernelFunction* >& kernelFunctions ) override
    {
      ParentPointParticle::_assignedKernelFunctions = kernelFunctions;

      ParentPointParticle::_nNodes = GradientEnhancedFiniteStrainParticle< nDim >::_assignedKernelFunctions.size();

      Eigen::Matrix< double, nDim, 1 > coords;
      ParentPointParticle::_mp.getCoordinatesAtCenter( coords.data() );

      Eigen::Matrix< double, nDim, nVertices > vertexCoordinates;
      this->getVertexCoordinates( vertexCoordinates.data() );

      ParentPointParticle::_N     = Eigen::MatrixXd::Zero( 1, ParentPointParticle::_nNodes );
      ParentPointParticle::_dN_dY = Eigen::MatrixXd::Zero( nDim, ParentPointParticle::_nNodes );
      for ( int i = 0; i < nDim; i++ )
        _d2N_dYdY[i] = Eigen::MatrixXd::Zero( nDim, ParentPointParticle::_nNodes );

      ParentPointParticle::_meshfreeApproximation.computeShapeFunctions( coords.data(),
                                                                         ParentPointParticle::_assignedKernelFunctions,
                                                                         ParentPointParticle::_N.data() );

      Eigen::MatrixXd NBoundary( 1, ParentPointParticle::_nNodes );
      Eigen::MatrixXd dN_dY_Boundary = Eigen::MatrixXd::Zero( nDim, ParentPointParticle::_nNodes );

      const LagrangeCellType cell( vertexCoordinates.data(), nDim * nVertices );

      for ( int i = 0; i < cell.getNumberOfFaces(); i++ ) {

        const Eigen::Matrix< double, nDim, 1 > n_dA       = cell.boundarySurfaceVector( i + 1 );
        const Eigen::Matrix< double, nDim, 1 > faceCenter = cell.getFaceCenterCoordinates( i + 1 );

        ParentPointParticle::_meshfreeApproximation
          .computeShapeFunctionsAndGradients( faceCenter.data(),
                                              ParentPointParticle::_assignedKernelFunctions,
                                              NBoundary.data(),
                                              dN_dY_Boundary.data() );

        ParentPointParticle::_dN_dY += n_dA * NBoundary;

        for ( int j = 0; j < nDim; j++ )
          _d2N_dYdY[j] += n_dA( j ) * dN_dY_Boundary;
      }

      const double VSmoothing = cell.volume();
      ParentPointParticle::_dN_dY /= VSmoothing;
      for ( int i = 0; i < nDim; i++ ) {
        _d2N_dYdY[i] /= VSmoothing;
      }

      ParentPointParticle::_T     = ParentPointParticle::_N;
      ParentPointParticle::_dT_dY = ParentPointParticle::_dN_dY;
    }

    void computePhysicsKernels( const double* dQ, double* fInt, double* dFInt_ddQ, double timeNew, double dT ) override;

    /// \brief Extract the second derivative of the shape function for a given node
    /// \param d2N_dYdY The second derivative of the shape function
    /// \param node The node for which the second derivative is extracted
    /// \return The second derivative of the shape function for the given node
    /// \details The second derivative is averaged over the two indices to ensure symmetry.
    inline TensorDD extract_d2N_dYdY_for_node( const std::array< Eigen::MatrixXd, nDim >& d2N_dYdY, int node ) const
    {
      TensorDD d2Nnode_dYdY( 0.0 );
      for ( int J = 0; J < nDim; J++ ) {
        for ( int I = 0; I < nDim; I++ ) {
          d2Nnode_dYdY( I, J ) += 0.5 * ( d2N_dYdY[J]( I, node ) + d2N_dYdY[I]( J, node ) );
        }
      }
      return d2Nnode_dYdY;
    }

    virtual void acceptStateAndPosition() override
    {
      ParentSQCNIParticle::acceptStateAndPosition();

      const auto   F    = this->_mp.dY_dX();
      const double detJ = det( F );

      // I2_ij = FiI * I2_IJ * F_jJ^T * detJ
      _momentsOfInertia_IntermediateReference = F % _momentsOfInertia_Undeformed % transpose( F ) * detJ;
    };
  };

  template < int nDim, int nVertices >
  GradientEnhancedFiniteStrainParticleSQCNIxNSNI< nDim, nVertices >::GradientEnhancedFiniteStrainParticleSQCNIxNSNI(
    int                                                  elementID,
    const double*                                        vertexCoordinates,
    int                                                  nVertexCoordinates,
    double                                               volume,
    const std::string&                                   materialName,
    const double*                                        materialProperties,
    int                                                  sizeMaterialProperties,
    const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation,
    ParentSQCNIParticle::SmoothingDomainUpdateType       smoothingVolumeUpdateType )
    : GradientEnhancedFiniteStrainParticleSQCNI< nDim, nVertices >( elementID,
                                                                    vertexCoordinates,
                                                                    nVertexCoordinates,
                                                                    volume,
                                                                    materialName,
                                                                    materialProperties,
                                                                    sizeMaterialProperties,
                                                                    approximation,
                                                                    smoothingVolumeUpdateType )
  {
    // second moments of the undeformed particle domain about its centroid
    // (for an affinely mapped cell this is identical to the former
    // detJ * 16/12 * J J^T quad formula)
    const Eigen::Matrix< double, nDim, nDim > secondMoments = this->_makeUndeformedCell().secondMoments();

    for ( int i = 0; i < nDim; i++ )
      for ( int j = 0; j < nDim; j++ )
        _momentsOfInertia_Undeformed( i, j ) = secondMoments( i, j );

    // the intermediate reference coincides with the undeformed one until the first accepted increment;
    // without this, the stabilization of the first increment is scaled by uninitialized memory
    _momentsOfInertia_IntermediateReference = _momentsOfInertia_Undeformed;
  }

  template < int nDim, int nVertices >
  void GradientEnhancedFiniteStrainParticleSQCNIxNSNI< nDim, nVertices >::computePhysicsKernels( const double* dQ,
                                                                                                 double*       fInt,
                                                                                                 double* dFInt_ddQ,
                                                                                                 double  timeNew,
                                                                                                 double  dT )
  {
    using namespace Marmot::FastorIndices;
    using namespace Fastor;
    using ink   = Fastor::Index< i_, n_, k_ >;
    using to_jk = Fastor::OIndex< j_, k_ >;
    using ijmM  = Index< i_, j_, m_, M_ >;
    using mMK   = Index< m_, M_, K_ >;
    using ijK   = Fastor::Index< i_, j_, K_ >;
    using ijkM  = Fastor::Index< i_, j_, k_, M_ >;
    using jK    = Fastor::Index< j_, K_ >;

    const auto& _nNodes = this->_nNodes;
    const auto& _N      = this->_N;
    const auto& _dN_dY  = this->_dN_dY;
    const auto& _T      = this->_T;
    const auto& _dT_dY  = this->_dT_dY;
    auto&       _mp     = this->_mp;

    const static TensorDD I(
      ( Eigen::Matrix< double, nDim, nDim >() << Eigen::Matrix< double, nDim, nDim >::Identity() ).finished().data() );

    constexpr int nodeBlockSize = nDim + 1;

    TensorD  du( 0.0 );
    TensorDD du_dY( 0.0 );

    double    dn = 0.0;
    TensorD   dn_dY( 0.0 );
    TensorDDD d2x_dYdY( 0.0 );

    for ( int B = 0; B < _nNodes; B++ ) {

      const int idxB_u = nodeBlockSize * B;
      const int idxB_n = nodeBlockSize * B + nDim;

      const double N_B        = _N( B );
      const auto   dN_B_dY    = TensorD( _dN_dY.col( B ).data() ); // works because ColumnMajor of Eigen
      const auto   d2N_B_dYdY = extract_d2N_dYdY_for_node( _d2N_dYdY, B );

      const auto dQU_B = TensorD( dQ + idxB_u );
      const auto dQN_B = dQ[idxB_n];

      du += N_B * dQU_B;
      dn += N_B * dQN_B;

      du_dY += einsum< i, j >( dQU_B, dN_B_dY );
      dn_dY += ( dQN_B * dN_B_dY );

      d2x_dYdY += einsum< i, jk >( dQU_B, d2N_B_dYdY );
    }

    _mp.prepareYourself( timeNew, dT );
    _mp.incrementDeformation( du, du_dY, dn );
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

    TensorD r_U( 0.0 );
    double  r_N( 0.0 );

    TensorDD k_UU( 0.0 );
    TensorD  k_UN( 0.0 );
    TensorD  k_NU( 0.0 );
    double   k_NN( 0.0 );

    const auto&  S           = _mp.response.S;
    const auto&  dLocalField = _mp.response.dL;
    const double c           = _mp.response.nonLocalRadius * _mp.response.nonLocalRadius;

    const double V0 = this->getVolumeUndeformed();

    const auto& t = _mp.tangents;

    Eigen::Map< Eigen::VectorXd > P( fInt, _nNodes * nodeBlockSize );
    Eigen::Map< Eigen::MatrixXd > K( dFInt_ddQ, _nNodes * nodeBlockSize, _nNodes * nodeBlockSize );

    const auto   dY_dx            = evaluate( inv( _mp.dx_dY() ) );
    const double detJIntermediate = determinant( _mp.dY_dX() );

    // Gradient of the Kirchhoff stress along Y, which is what the NSNI stabilization charges.
    // Only the deformation-gradient term survives here: the micropolar reference also carried
    // dS_dW and dS_ddWdY, and like it we do NOT include the dS_dN contribution of the nonlocal
    // field.  Note that the tangent of the stabilization below is APPROXIMATE: it differentiates
    // d2x_dYdY but not dS_dDeltaF itself, i.e. it omits d2tau/dF2 (not exposed by the material
    // interface) and the dependence of dS_dDeltaF on the nonlocal field (the damage). The nonlocal
    // rows are exact; see the module test for the measured size of the omitted terms.
    const auto dS_dY = evaluate( einsum< ijmM, mMK >( t.dS_dDeltaF, d2x_dYdY ) );

    // clang-format off
    for ( int A = 0; A < _nNodes; A++ ) {

      const double T_A = _T( A );
      const auto  dT_A_dY = TensorMap< const double, nDim >( _dT_dY.col( A ).data() );
      const TensorD dT_A_dx = einsum< ji, j >( dY_dx , dT_A_dY );
      const TensorD dT_A_dX = einsum< ji, j >( _mp.dY_dX(), dT_A_dY );

      const TensorD dn_dX = einsum< ji, j >( _mp.dY_dX(), dn_dY );

      const int idxA_u = nodeBlockSize * A;
      const int idxA_n = nodeBlockSize * A + nDim;

      r_U = ( +einsum< i, ij >( dT_A_dx, S ) ) * V0;
      r_N = evaluate( ( T_A * dn + c * einsum< i, i >( dT_A_dX, dn_dX ) - T_A * dLocalField ) * V0 ) .toscalar();

      // add inertia
      r_U += density0 * a * T_A * V0;

      const auto d2NA_dYdY = extract_d2N_dYdY_for_node( _d2N_dYdY, A );

      const TensorDD d2NA_dYdY_x_MOIScaled = einsum< ij, jk >( d2NA_dYdY, _momentsOfInertia_IntermediateReference ) / detJIntermediate;
      const TensorDD d2NA_dxdY_x_MOIScaled = einsum< ji, jk >( dY_dx, d2NA_dYdY_x_MOIScaled ) ;
      const TensorD dTA_dY_x_MOIScaled = einsum< j, jk >( dT_A_dY, _momentsOfInertia_IntermediateReference ) / detJIntermediate;

      TensorD rU_Stab = einsum< iK, ijK >( d2NA_dxdY_x_MOIScaled, dS_dY );

      r_U += rU_Stab;

      {
        using namespace Eigen;
        P.template segment< nDim >( idxA_u ) += Map< Matrix< double, nDim, 1 > >( r_U.data() );
        P( idxA_n ) += r_N;
      }

      for ( int B = 0; B < _nNodes; B++ ) {

        const int idxB_u = nodeBlockSize * B;
        const int idxB_n = nodeBlockSize * B + nDim;

        const double                 N_B     = _N( B );
        const auto dN_B_dY = TensorMap< const double, nDim >( _dN_dY.col(B).data() );
        const auto dN_B_dx = evaluate( einsum< ji, j >( dY_dx, dN_B_dY ) ); //
        const auto dN_B_dX = evaluate( einsum< ji, j >( _mp.dY_dX(), dN_B_dY ) ); // no dependence on current deformations!
        const auto d2NB_dYdY = extract_d2N_dYdY_for_node( _d2N_dYdY, B );

        // aux stiffness tensors
        const auto dS_dqU_B = evaluate ( + einsum < ijkl, l > ( t.dS_dDeltaF, dN_B_dY )                                            );
        const auto dS_dqN_B = evaluate (                      ( t.dS_dN*       N_B    )                                            );
        const auto dL_dqU_B = evaluate ( + einsum <   kl, l > ( t.dL_dDeltaF, dN_B_dY )                                             );

        k_UU  = ( + einsum< i, ijk        > ( dT_A_dx, dS_dqU_B )                                                       ) * V0;
        k_UN  = ( + einsum< i,  ij        > ( dT_A_dx, dS_dqN_B )                                                       ) * V0;

        k_NU  = (                                                 - ( T_A * dL_dqU_B )                                  ) * V0;
        k_NN  = ( + T_A * N_B + inner( dT_A_dX, dN_B_dX ) *  c                                                          ) * V0;

        k_UU += ( - einsum< k, ij, i, to_jk >( dT_A_dx, S, dN_B_dx ) ) * V0;

        k_UU += density0 * da_ddu * T_A * N_B * V0;

        // clang-format on
        const TensorDDDD d2S_dqu_dY = einsum< ijkM, MK >( t.dS_dDeltaF, d2NB_dYdY );

        const TensorDD dNB_dx_times_dS_dY = einsum< ijK, i >( dS_dY, dN_B_dx );

        const TensorDD dRU_Stab_dU_1   = einsum< iK, ijkK, to_jk >( d2NA_dxdY_x_MOIScaled, d2S_dqu_dY );
        const TensorDD dRU_Stab_x_dU_2 = -einsum< kK, jK, to_jk >( d2NA_dxdY_x_MOIScaled, dNB_dx_times_dS_dY );

        k_UU += dRU_Stab_dU_1 + dRU_Stab_x_dU_2;

        // clang-format off
        {
            using namespace Eigen;
            // TODO: check if we can use transpose instead of torowmajor:
            K.template block< nDim, nDim >( idxA_u, idxB_u ) += Map< Matrix< double, nDim, nDim > >( torowmajor( k_UU ).data() );
            K.template block< nDim,    1 >( idxA_u, idxB_n ) += Map< Matrix< double, nDim,    1 > >( torowmajor( k_UN ).data() );
            K.template block<    1, nDim >( idxA_n, idxB_u ) += Map< Matrix< double,    1, nDim > >( torowmajor( k_NU ).data() );
            K                             ( idxA_n, idxB_n ) +=                                                  k_NN           ;
        }
      }
    }
    // clang-format on
  }

} // namespace Marmot::Meshfree
