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

#include "Marmot/DisplacementFiniteStrainULElement.h"
#include "Marmot/MarmotDeformationMeasures.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include <Fastor/expressions/linalg_ops/unary_trans_op.h>
#include <Fastor/tensor/Tensor.h>
#include <Fastor/tensor_algebra/einsum_explicit.h>
#include <Fastor/tensor_algebra/indicial.h>

namespace Marmot::Elements {

  /**
   * @class Marmot::Elements::DisplacementFiniteStrainFBarElement
   * @tparam nDim   Number of spatial dimensions. Only @c nDim=3 (Solid section) is currently
   * implemented; other values throw at runtime.
   * @tparam nNodes Number of nodes of the element. Intended for the low-order shapes (Tetra4,
   * Hexa8) that suffer most from volumetric locking.
   *
   * @brief F-bar displacement element (de Souza Neto, Perić, Dutko, Owen, 1996) for locking-free
   * analysis of nearly incompressible finite-strain materials on low-order shapes.
   *
   * @details Rather than evaluating the material at the quadrature point's own deformation
   * gradient \f$\boldsymbol F\f$, this element evaluates it at the F-bar modified deformation
   * gradient
   * \f[
   *   \bar{\boldsymbol F} = \left(\frac{J_0}{J}\right)^{1/3} \boldsymbol F,
   * \f]
   * where \f$J=\det\boldsymbol F\f$ is the quadrature point's own volume ratio and \f$J_0\f$ is
   * the volume ratio sampled once at the element centroid
   * (DeformationMeasures::FBarDeformationGradient()). This decouples the volumetric response
   * (effectively sampled only once per element) from the deviatoric response (still sampled at
   * every quadrature point), which alleviates the volumetric locking that a standard low-order
   * displacement element suffers from with nearly incompressible materials - without adding any
   * degrees of freedom.
   *
   * Unlike Marmot::Elements::DisplacementPressureFiniteStrainElement, this is a locking
   * *mitigation* for a material with genuine (if large) volumetric stiffness - not a mechanism for
   * *exact* incompressibility. It must be paired with a material that has a volumetric energy term
   * \f$U(J)\f$ (e.g. Marmot::Materials::CompressibleNeoHooke with a large bulk modulus \f$K\f$
   * relative to the shear modulus \f$G\f$); it is not compatible with the purely isochoric
   * materials in Marmot::Materials::Ogden, Marmot::Materials::IncompressibleNeoHooke or
   * Marmot::Materials::IncompressibleMooneyRivlin, which have no volumetric stiffness for F-bar to
   * redistribute in the first place.
   *
   * Because \f$J_0\f$ depends on every node's degrees of freedom (not just those of the shape
   * functions local to a quadrature point), the consistent tangent has an additional correction
   * term beyond the standard material and geometric stiffness - see computeKernels() for the
   * derivation. Only computeKernels() and computeKernelsExplicit() are overridden; everything else
   * (geometry, quadrature, DOF layout, loads, inertia) is inherited unchanged from
   * DisplacementFiniteStrainULElement, since F-bar does not introduce any new fields or degrees of
   * freedom.
   */
  template < int nDim, int nNodes >
  class DisplacementFiniteStrainFBarElement : public DisplacementFiniteStrainULElement< nDim, nNodes > {

  public:
    using Parent = DisplacementFiniteStrainULElement< nDim, nNodes >;
    using Parent::Parent;

    using SectionType   = typename Parent::SectionType;
    using XiSized       = typename Parent::XiSized;
    using JacobianSized = typename Parent::JacobianSized;
    using dNdXiSized    = typename Parent::dNdXiSized;
    using KSizedMatrix  = typename Parent::KSizedMatrix;
    using Material      = typename Parent::Material;

    /**
     * @brief Compute the negative element residual vector and the consistent F-bar tangent.
     *
     * @details At each quadrature point, the deformation gradient \f$\boldsymbol F\f$ and its
     * volume ratio \f$J\f$ are computed as usual, together with the deformation gradient
     * \f$\boldsymbol F_0\f$ and volume ratio \f$J_0\f$ at the element centroid (computed once per
     * call, shared by all quadrature points). The material is evaluated at
     * \f$\bar{\boldsymbol F}=(J_0/J)^{1/3}\boldsymbol F\f$, giving \f$\bar{\boldsymbol\tau}\f$ and
     * \f$\boldsymbol a=\partial\bar{\boldsymbol\tau}/\partial\bar{\boldsymbol F}\f$. The residual
     * uses \f$\bar{\boldsymbol\tau}\f$ together with the spatial gradient operator built from the
     * *actual* (unmodified) \f$\boldsymbol F\f$:
     * \f[
     *   \boldsymbol R_u = \sum_{qp} (\nabla_x\boldsymbol N)^{\!\top}\,\bar{\boldsymbol\tau}\,J_0
     *   w_{qp}.
     * \f]
     * Writing \f$\alpha=(J_0/J)^{1/3}\f$, the tangent follows from
     * \f$\partial\bar{\boldsymbol F}/\partial\boldsymbol q_B = \alpha\,\partial\boldsymbol
     * F/\partial\boldsymbol q_B + \boldsymbol F\otimes\partial\alpha/\partial\boldsymbol q_B\f$
     * with \f$\partial\alpha/\partial\boldsymbol q_B = \tfrac13\alpha\,(\nabla_x^0\boldsymbol
     * N_B-\nabla_x\boldsymbol N_B)\f$ (using \f$\partial J/\partial\boldsymbol q_B=J\,\nabla_x
     * \boldsymbol N_B\f$ at both the quadrature point and the centroid). This gives three
     * contributions to \f$K_{uu}\f$ per quadrature point: the standard geometric term (using
     * \f$\bar{\boldsymbol\tau}\f$), the standard material term scaled by \f$\alpha\f$ (using
     * \f$\boldsymbol a\f$ contracted with the quadrature point's own shape function gradients),
     * and an additional F-bar correction term
     * \f[
     *   \tfrac{\alpha}{3}\, \big(\boldsymbol a:\boldsymbol F\big)\cdot(\nabla_x\boldsymbol N_A) \otimes
     *   \big(\nabla_x^0\boldsymbol N_B - \nabla_x\boldsymbol N_B\big)
     * \f]
     * coupling every quadrature point's stiffness to the centroid's sensitivity to every node's
     * degrees of freedom - the term a naive "just substitute \f$\bar{\boldsymbol F}\f$ for
     * \f$\boldsymbol F\f$" implementation would miss.
     */
    void computeKernels( const double* QTotal, const double* dQ, double* Pe, double* Ke, double time, double dT );

    /**
     * @brief Compute the negative element residual vector only (no tangent), for explicit time
     * integration, using the same F-bar substitution as computeKernels().
     */
    void computeKernelsExplicit( const double* QTotal, const double* dQ, double* Pe, double time, double dT );
  };

  template < int nDim, int nNodes >
  void DisplacementFiniteStrainFBarElement< nDim, nNodes >::computeKernels( const double* qTotal,
                                                                            const double* dQ,
                                                                            double*       rightHandSide,
                                                                            double*       stiffnessMatrix,
                                                                            double        time,
                                                                            double        dT )
  {
    if constexpr ( nDim != 3 ) {
      throw std::runtime_error( "DisplacementFiniteStrainFBarElement currently only supports nDim=3 (Solid section)" );
    }
    else {
      using namespace Fastor;
      using namespace Marmot;
      using namespace Marmot::FastorIndices;
      using namespace Marmot::FastorStandardTensors;

      const auto& I = Spatial3D::I;

      const auto qU_np = TensorMap< const double, nNodes, nDim >( qTotal );

      TensorMap< double, nNodes, nDim >            r_U( rightHandSide );
      Tensor< double, nDim, nNodes, nDim, nNodes > k_UU( 0.0 );

      // centroid quantities (F0, J0, and the corresponding pushed-forward shape function
      // gradients), shared by every quadrature point and recomputed fresh each call - cheap
      // relative to the O(nNodes^2) stiffness assembly, and avoids touching the base class's
      // per-quadrature-point data layout.
      const auto          dNdXi_0 = this->dNdXi( XiSized::Zero() );
      const JacobianSized Jac_0   = this->Jacobian( dNdXi_0 );
      const auto          dNdX_0_ = this->dNdX( dNdXi_0, Jac_0.inverse() );
      const auto          dNdX_0  = Tensor< double, nDim, nNodes >( dNdX_0_.data(), ColumnMajor );

      const auto   F_0    = evaluate( einsum< Ai, jA >( qU_np, dNdX_0 ) + I );
      const double J_0    = determinant( F_0 );
      const auto   dNdx_0 = evaluate( einsum< ji, jA >( inv( F_0 ), dNdX_0 ) );

      for ( auto& qp : this->qps ) {

        const auto& dNdX_ = qp.dNdX;
        const auto  dNdX  = Tensor< double, nDim, nNodes >( dNdX_.data(), ColumnMajor );

        const auto   F_qp  = evaluate( einsum< Ai, jA >( qU_np, dNdX ) + I );
        const double J_qp  = determinant( F_qp );
        const double alpha = pow( J_0 / J_qp, 1. / 3. );

        const auto Fbar = ContinuumMechanics::DeformationMeasures::FBarDeformationGradient( F_qp, J_qp, J_0 );

        const typename Material::template Deformation< nDim > deformation{ Fbar };
        const typename Material::TimeIncrement                timeIncrement{ time, dT };
        typename Material::template ConstitutiveResponse< nDim >
          response( Tensor< double, nDim, nDim >( qp.managedStateVars->stress.data(), ColumnMajor ),
                    qp.managedStateVars->elasticEnergy / qp.J0xW,
                    qp.managedStateVars->dissipation / qp.J0xW,
                    qp.managedStateVars->materialStateVars.data() );
        typename Material::template AlgorithmicModuli< nDim > tangents;

        qp.material->computeStress( response, tangents, deformation, timeIncrement );

        qp.managedStateVars->stress            = Marmot::mapEigenToFastor( response.tau ).reshaped();
        qp.managedStateVars->F                 = Marmot::mapEigenToFastor( F_qp ).reshaped();
        qp.managedStateVars->elasticEnergy     = response.elasticEnergyDensity * qp.J0xW;
        qp.managedStateVars->dissipation       = response.dissipation * qp.J0xW;
        qp.managedStateVars->totalStrainEnergy = ( response.elasticEnergyDensity + response.dissipation ) * qp.J0xW;

        // spatial gradient operator from the *actual* F, not Fbar
        const auto dNdx = evaluate( einsum< ji, jA >( inv( F_qp ), dNdX ) );

        const double& J0xW   = qp.J0xW;
        const auto&   taubar = response.tau;
        const auto&   a      = tangents.dTau_dF; // dTaubar/dFbar, from the material

        r_U += ( +einsum< iA, ij >( dNdx, taubar ) ) * J0xW;

        // standard geometric term (using taubar and the real dNdx - identical structure to the
        // plain UL element)
        k_UU += ( -einsum< kA, ij, iB, to_jAkB >( dNdx, taubar, dNdx ) ) * J0xW;

        // standard material term, scaled by alpha (since dFbar/dF|_(F0 fixed) = alpha * dF/dF)
        const auto dTau_dqU = evaluate( alpha * einsum< ijkl, lB >( a, dNdX ) );
        k_UU += ( +einsum< iA, ijkB, to_jAkB >( dNdx, dTau_dqU ) ) * J0xW;

        // F-bar correction term, from the dependence of alpha on both F (local) and F0
        // (centroid, i.e. every node's DOFs): with G_ij = a_ijkl F_kl and H_jA = dNdx_iA G_ij,
        //   dtaubar_ij/dqU_Bn += (alpha/3) * G_ij * (dNdx0_nB - dNdx_nB)
        // contracted with dNdx_iA on the residual gives the correction to k_UU below.
        Tensor< double, nDim, nDim > G( 0.0 );
        for ( int i = 0; i < nDim; ++i )
          for ( int j = 0; j < nDim; ++j )
            for ( int k = 0; k < nDim; ++k )
              for ( int l = 0; l < nDim; ++l )
                G( i, j ) += a( i, j, k, l ) * F_qp( k, l );

        Tensor< double, nDim, nNodes > H( 0.0 );
        for ( int j = 0; j < nDim; ++j )
          for ( int A = 0; A < nNodes; ++A )
            for ( int i = 0; i < nDim; ++i )
              H( j, A ) += dNdx( i, A ) * G( i, j );

        for ( int j = 0; j < nDim; ++j )
          for ( int A = 0; A < nNodes; ++A )
            for ( int n = 0; n < nDim; ++n )
              for ( int B = 0; B < nNodes; ++B )
                k_UU( j, A, n, B ) += ( alpha / 3. ) * H( j, A ) * ( dNdx_0( n, B ) - dNdx( n, B ) ) * J0xW;
      }

      using namespace Eigen;

      Map< KSizedMatrix > K( stiffnessMatrix );
      K.template block< Parent::bsU, Parent::bsU >( Parent::idxU, Parent::idxU ) += Map<
        Matrix< double, Parent::bsU, Parent::bsU > >( torowmajor( k_UU ).data() );
    }
  }

  template < int nDim, int nNodes >
  void DisplacementFiniteStrainFBarElement< nDim, nNodes >::computeKernelsExplicit( const double* qTotal,
                                                                                    const double* dQ,
                                                                                    double*       rightHandSide,
                                                                                    double        time,
                                                                                    double        dT )
  {
    if constexpr ( nDim != 3 ) {
      throw std::runtime_error( "DisplacementFiniteStrainFBarElement currently only supports nDim=3 (Solid section)" );
    }
    else {
      using namespace Fastor;
      using namespace Marmot;
      using namespace Marmot::FastorIndices;
      using namespace Marmot::FastorStandardTensors;

      const auto& I = Spatial3D::I;

      const auto                        qU_np = TensorMap< const double, nNodes, nDim >( qTotal );
      TensorMap< double, nNodes, nDim > r_U( rightHandSide );

      const auto          dNdXi_0 = this->dNdXi( XiSized::Zero() );
      const JacobianSized Jac_0   = this->Jacobian( dNdXi_0 );
      const auto          dNdX_0_ = this->dNdX( dNdXi_0, Jac_0.inverse() );
      const auto          dNdX_0  = Tensor< double, nDim, nNodes >( dNdX_0_.data(), ColumnMajor );

      const auto   F_0 = evaluate( einsum< Ai, jA >( qU_np, dNdX_0 ) + I );
      const double J_0 = determinant( F_0 );

      for ( auto& qp : this->qps ) {

        const auto& dNdX_ = qp.dNdX;
        const auto  dNdX  = Tensor< double, nDim, nNodes >( dNdX_.data(), ColumnMajor );

        const auto   F_qp = evaluate( einsum< Ai, jA >( qU_np, dNdX ) + I );
        const double J_qp = determinant( F_qp );

        const auto Fbar = ContinuumMechanics::DeformationMeasures::FBarDeformationGradient( F_qp, J_qp, J_0 );

        const typename Material::template Deformation< nDim > deformation{ Fbar };
        const typename Material::TimeIncrement                timeIncrement{ time, dT };
        typename Material::template ConstitutiveResponse< nDim >
          response( Tensor< double, nDim, nDim >( qp.managedStateVars->stress.data(), ColumnMajor ),
                    qp.managedStateVars->elasticEnergy / qp.J0xW,
                    qp.managedStateVars->dissipation / qp.J0xW,
                    qp.managedStateVars->materialStateVars.data() );

        qp.material->computeStressExplicit( response, deformation, timeIncrement );

        qp.managedStateVars->stress            = Marmot::mapEigenToFastor( response.tau ).reshaped();
        qp.managedStateVars->F                 = Marmot::mapEigenToFastor( F_qp ).reshaped();
        qp.managedStateVars->elasticEnergy     = response.elasticEnergyDensity * qp.J0xW;
        qp.managedStateVars->dissipation       = response.dissipation * qp.J0xW;
        qp.managedStateVars->totalStrainEnergy = ( response.elasticEnergyDensity + response.dissipation ) * qp.J0xW;

        const auto dNdx = evaluate( einsum< ji, jA >( inv( F_qp ), dNdX ) );

        r_U += ( +einsum< iA, ij >( dNdx, response.tau ) ) * qp.J0xW;
      }
    }
  }

} // namespace Marmot::Elements
