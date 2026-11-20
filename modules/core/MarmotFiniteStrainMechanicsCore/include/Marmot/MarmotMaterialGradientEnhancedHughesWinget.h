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

#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMaterialGradientEnhancedHypoElasticWrapper.h"
#include "Marmot/MarmotMaterialHughesWinget.h"
#include "Marmot/MarmotNumericalDifferentiation.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"
#include <cmath>
#include <cstring>
#include <limits>
#include <memory>
#include <type_traits>
#include <vector>

namespace Marmot::Materials {

  /**
   * @class GradientEnhancedHughesWingetWrapper
   * @brief Decorator presenting a small-strain gradient-enhanced material as a finite-strain one.
   *
   * @tparam BaseMaterialType A concrete class derived from Marmot::Materials::MaterialGradientEnhancedHypoElastic<n>;
   *                          the wrapper has the same n nonlocal variables.
   * @tparam tangentMode      How the algorithmic tangent is evaluated, see Marmot::Materials::HughesWingetTangent.
   *
   * The gradient-enhanced sibling of HughesWingetWrapper, and it works the same way: the wrapped material
   * is driven by the objective strain increment of the Hughes-Winget algorithm, evaluated on the mid-step
   * configuration,
   * @f[
   *   \Delta\boldsymbol{l} = \left(\boldsymbol{F}^{(n+1)} - \boldsymbol{F}^{(n)}\right)
   *                          \boldsymbol{F}_{\text{mid}}^{-1}, \qquad
   *   \boldsymbol{F}_{\text{mid}} = \tfrac{1}{2}\left(\boldsymbol{F}^{(n)} + \boldsymbol{F}^{(n+1)}\right),
   * @f]
   * split into @f$ \Delta\boldsymbol{\varepsilon} @f$ and @f$ \Delta\boldsymbol{\Omega} @f$, with the stress
   * of the previous increment rotated forward by the Cayley transform of @f$ \Delta\boldsymbol{\Omega} @f$
   * and the resulting Cauchy stress pushed to the Kirchhoff stress @f$ \boldsymbol{\tau} = J\,\boldsymbol{\sigma} @f$.
   *
   * On top of that it carries the nonlocal fields across the two interfaces, per nonlocal variable K:
   *
   * | this interface reports | from the wrapped material |
   * |---|---|
   * | @c response.L | @c KLocal, totals, which is what both sides mean |
   * | @c response.nonLocalRadius | @f$ \sqrt{c_K} @f$ |
   * | @c tangents.dTau_dN[K] | @f$ J\,\partial\boldsymbol{\sigma}/\partial\bar{N}_K @f$ |
   * | @c tangents.dL_dF[K] | @f$ \partial K^{\text{local}}_K/\partial\Delta\boldsymbol{\varepsilon} @f$ |
   * | @c tangents.dL_dN | @c dKLocalddK |
   *
   * @c dL_dF is chained through the Hughes-Winget kinematics. Each row is kept inside the 120-column
   * limit deliberately: `ReflowComments` would otherwise wrap a longer one onto the next line and
   * split the table in two.
   *
   * @note The wrapped material wants the nonlocal fields **and their increment**, while this interface hands
   * over only the totals. The increment is therefore formed against the value of the last call carried in this
   * wrapper's own state (@c HughesWinget_N_n), see GradientEnhancedHypoElasticWrapper.
   *
   * @note **The wrapped material must be incremental in stress.** It has to update the Cauchy stress it
   * is handed, not recompute one from a stored strain: the wrapper's whole mechanism is to hand over the
   * forward-rotated stress of the last increment, and a material that ignores that argument discards the
   * rotation along with it. `AT2PhaseField` is a counter-example: it carries a six-component `strain`
   * state (`AT2PhaseField.h:60`) and returns @f$ g(\varphi)\,\mathbb{C}:\boldsymbol{\varepsilon} @f$ from
   * it (`AT2PhaseField.cpp:71`), so its true
   * @f$ \partial\boldsymbol{\sigma}^{(n+1)}/\partial\boldsymbol{\sigma}_{\text{rot}} @f$ is zero
   * rather than the identity. Measured: wrapped, a rigid rotation on top of an anisotropic stretch leaves
   * its Kirchhoff stress bit-for-bit **unrotated**. It is therefore deliberately not registered with this
   * wrapper, and a total-strain-based model must not be either.
   *
   * @note **Only the stress is rotated.** Tensor-valued internal variables of the wrapped material are
   * passed through untouched and are therefore *not* objective under large incremental rotations. A model
   * whose internal state is entirely scalar is unaffected.
   *
   * @note The interaction @f$ c @f$ is reported through its square root, and its derivative
   * @f$ \partial c/\partial\bar{N} @f$ **cannot be forwarded at all**: the finite-strain interface has no
   * slot for it. Harmless for a material whose @f$ c @f$ is constant (`AT2PhaseField` is); a material
   * built on MarmotDecreasingInteractions would keep a correct residual but lose that term of its
   * consistent tangent.
   *
   * @note The nonlocal balance is formulated in the material configuration by the consumers of this
   * interface, so @f$ c @f$ is handed over as the material constant the wrapped model reports, neither
   * pushed forward nor weighted by @f$ J @f$.
   */
  template < typename BaseMaterialType, HughesWingetTangent tangentMode = HughesWingetTangent::Analytic >
  class GradientEnhancedHughesWingetWrapper : public GradientEnhancedHypoElasticWrapper< BaseMaterialType > {

    using Wrapper = GradientEnhancedHypoElasticWrapper< BaseMaterialType >;
    using typename Wrapper::BaseMaterial;
    using Wrapper::baseMaterial;
    using Wrapper::baseMaterialSlot;
    using Wrapper::fromEigen;
    using Wrapper::nNonlocal;
    using Wrapper::toEigen;

  public:
    using typename Wrapper::NonlocalVector;
    /// @brief Constitutive response, see Marmot::Materials::MaterialGradientEnhancedFiniteStrain::ConstitutiveResponse.
    template < int nDim >
    using ConstitutiveResponse = typename Wrapper::Interface::template ConstitutiveResponse< nDim >;
    /// @brief Algorithmic moduli, see Marmot::Materials::MaterialGradientEnhancedFiniteStrain::AlgorithmicModuli.
    template < int nDim >
    using AlgorithmicModuli = typename Wrapper::Interface::template AlgorithmicModuli< nDim >;
    /// @brief Deformation, see Marmot::Materials::MaterialGradientEnhancedFiniteStrain::Deformation.
    template < int nDim >
    using Deformation = typename Wrapper::Interface::template Deformation< nDim >;
    /// @brief Time increment, see Marmot::Materials::MaterialGradientEnhancedFiniteStrain::TimeIncrement.
    using TimeIncrement = typename Wrapper::Interface::TimeIncrement;

  protected:
    /// Name of the state layout slot holding the deformation gradient of the last accepted increment.
    static constexpr const char* deformationGradientSlot = "HughesWinget_F_n";
    /// Name of the state layout slot holding the Cauchy stress of the last accepted increment.
    static constexpr const char* stressSlot = "HughesWinget_sigma_n";
    /// Name of the state layout slot holding the nonlocal fields of the last accepted increment.
    static constexpr const char* nonlocalFieldSlot = "HughesWinget_N_n";

  public:
    /**
     * @brief Construct the wrapper and the wrapped material.
     * @param[in] matProperties_       Material property values; forwarded **unchanged** to the wrapped material.
     * @param[in] nMaterialProperties_ Number of material property values.
     * @param[in] materialNumber_      Unique identifier for this material instance.
     */
    GradientEnhancedHughesWingetWrapper( const double* matProperties_, int nMaterialProperties_, int materialNumber_ )
      : Wrapper( matProperties_, nMaterialProperties_, materialNumber_ )
    {
      initializeStateLayout();
    }

    virtual ~GradientEnhancedHughesWingetWrapper() = default;

    // The four-argument override below would otherwise hide the base class' five-argument
    // eigen-deformation overload for a caller holding a GradientEnhancedHughesWingetWrapper<...>
    // rather than a Marmot::Materials::MaterialGradientEnhancedFiniteStrain*.
    using Wrapper::Interface::computeStress;

    /**
     * @brief Register the state layout: the carried-over deformation gradient, Cauchy stress and nonlocal
     *        field, followed by the state variables of the wrapped material.
     */
    void initializeStateLayout()
    {
      this->stateLayout.add( deformationGradientSlot, 9 );
      this->stateLayout.add( stressSlot, 6 );
      this->stateLayout.add( nonlocalFieldSlot, nNonlocal );
      this->stateLayout.add( baseMaterialSlot, baseMaterial->getNumberOfRequiredStateVars() );
      this->stateLayout.finalize();
    }

    /**
     * @brief Initialise the state: @f$ \boldsymbol{F}^{(n)} = \boldsymbol{I} @f$,
     *        @f$ \boldsymbol{\sigma}^{(n)} = \boldsymbol{0} @f$, @f$ \bar{\boldsymbol{N}}^{(n)} = \boldsymbol{0} @f$,
     *        plus the wrapped material's own state.
     * @param[in,out] stateVars  Pointer to the state variable array.
     * @param[in]     nStateVars Number of state variables.
     */
    void initializeYourself( double* stateVars, int nStateVars ) override
    {
      using namespace Marmot::Math::Tensors::FastorTensors;

      TensorMap33d Fn = this->stateLayout.template getAs< TensorMap33d >( stateVars, deformationGradientSlot );
      std::memcpy( Fn.data(), Spatial3D::I.data(), 9 * sizeof( double ) );

      double* sigmaN = this->stateLayout.getPtr( stateVars, stressSlot );
      std::fill( sigmaN, sigmaN + 6, 0.0 );

      double* Nn = this->stateLayout.getPtr( stateVars, nonlocalFieldSlot );
      std::fill( Nn, Nn + nNonlocal, 0.0 );

      this->initializeBaseMaterial( stateVars );
    }

    /**
     * @brief Objective stress update of the wrapped small-strain gradient-enhanced material.
     * @param[in,out] response      Constitutive response; carries the Kirchhoff stress, the local driving
     *                              force, the nonlocal radius and the state variables.
     * @param[out]    tangents      Algorithmic moduli.
     * @param[in]     deformation   Current deformation gradient and nonlocal field.
     * @param[in]     timeIncrement Current (pseudo-)time and (pseudo-)time increment.
     */
    void computeStress( ConstitutiveResponse< 3 >& response,
                        AlgorithmicModuli< 3 >&    tangents,
                        const Deformation< 3 >&    deformation,
                        const TimeIncrement&       timeIncrement ) const override
    {
      if constexpr ( tangentMode == HughesWingetTangent::Numerical ) {
        computeNumericalTangent( response, tangents, deformation, timeIncrement );
        // The unperturbed evaluation runs last, so the committed state is the true one.
        computeStressCore( response, tangents, deformation, timeIncrement, false );
      }
      else {
        computeStressCore( response, tangents, deformation, timeIncrement, true );
      }
    }

    /**
     * @brief Explicit variant of @ref computeStress; no algorithmic tangent is evaluated.
     * @param[in,out] response      Constitutive response.
     * @param[in]     deformation   Current deformation gradient and nonlocal field.
     * @param[in]     timeIncrement Current (pseudo-)time and (pseudo-)time increment.
     */
    void computeStressExplicit( ConstitutiveResponse< 3 >& response,
                                const Deformation< 3 >&    deformation,
                                const TimeIncrement&       timeIncrement ) const override
    {
      AlgorithmicModuli< 3 > unused;
      computeStressCore( response, unused, deformation, timeIncrement, false );
    }

  protected:
    /**
     * @brief Recover @f$ \boldsymbol{S} = \partial\boldsymbol{\sigma}^{(n+1)} /
     *        \partial\boldsymbol{\sigma}_{\text{rot}} @f$ by forward differences.
     * @param[in] sigmaRotVoigt Rotated stress that was handed to the wrapped material.
     * @param[in] baseStateOld  Wrapped material state *before* the unperturbed evaluation.
     * @param[in] inc           Increment handed to the wrapped material.
     * @param[in] elasticEnergyDensity Elastic energy density of the unperturbed evaluation.
     * @param[in] dissipation   Dissipation of the unperturbed evaluation.
     * @return The sensitivity in Voigt form, with its three shear **columns** halved so that the result
     *         may be contracted as a full fourth-order tensor without double counting.
     */
    Marmot::Math::Matrix6d computeStressSensitivity( const Marmot::Math::Vector6d&           sigmaRotVoigt,
                                                     const std::vector< double >&            baseStateOld,
                                                     const typename BaseMaterial::increment& inc,
                                                     double                                  elasticEnergyDensity,
                                                     double                                  dissipation ) const
    {
      const int             nBase = baseMaterial->getNumberOfRequiredStateVars();
      std::vector< double > scratch( nBase );

      // A fresh copy of the wrapped material's pre-update state for every evaluation, including the
      // one MarmotMathCore's forwardDifference() takes at sigmaRotVoigt itself, so each probe -- and
      // the baseline -- starts from the same point regardless of evaluation order.
      auto evaluate = [&]( const Eigen::VectorXd& sigmaRot ) -> Eigen::VectorXd {
        if ( nBase > 0 )
          std::memcpy( scratch.data(), baseStateOld.data(), nBase * sizeof( double ) );

        typename BaseMaterial::response perturbed;
        perturbed.stress = sigmaRot;
        // BaseMaterial::response has no default member initialisers; a material that reads or
        // accumulates into KLocal/c would otherwise see garbage.
        perturbed.KLocal.setZero();
        perturbed.c.setZero();
        perturbed.stateVars            = nBase > 0 ? scratch.data() : nullptr;
        perturbed.elasticEnergyDensity = elasticEnergyDensity;
        perturbed.dissipation          = dissipation;

        // The tangent this probe would return is discarded -- the probe IS the tangent. Asking for
        // it costs a full consistent tangent per column per quadrature point.
        baseMaterial->computeStressExplicit( perturbed, inc );
        return perturbed.stress;
      };

      Marmot::Math::Matrix6d S = Marmot::Math::Differentiation::Numerical::forwardDifference( evaluate, sigmaRotVoigt );

      // Halve the shear columns: contracting S as a fourth-order tensor sums over both the (m,n) and
      // (n,m) entries of a symmetric stress, which would otherwise count each off-diagonal twice.
      S.rightCols( 3 ) *= 0.5;

      return S;
    }

    /**
     * @brief The Hughes-Winget update proper, optionally including the analytic algorithmic tangent.
     * @param[in,out] response       Constitutive response.
     * @param[out]    tangents       Algorithmic moduli; only written when @p computeTangent is true.
     * @param[in]     deformation    Current deformation gradient and nonlocal field.
     * @param[in]     timeIncrement  Current (pseudo-)time and (pseudo-)time increment.
     * @param[in]     computeTangent Whether the analytic tangent is to be evaluated.
     */
    void computeStressCore( ConstitutiveResponse< 3 >& response,
                            AlgorithmicModuli< 3 >&    tangents,
                            const Deformation< 3 >&    deformation,
                            const TimeIncrement&       timeIncrement,
                            bool                       computeTangent ) const
    {
      using namespace Fastor;
      using namespace Marmot::Math::Tensors::FastorTensors;
      using namespace Marmot::Math::Tensors::FastorIndices;
      using namespace Marmot::ContinuumMechanics::Voigt;

      const Tensor33d& Ident = Spatial3D::I;

      TensorMap33d Fn_ref = this->stateLayout.template getAs< TensorMap33d >( response.stateVars,
                                                                              deformationGradientSlot );
      Eigen::Map< Marmot::Math::Vector6d >
        sigmaN_n = this->stateLayout.template getAs< Eigen::Map< Marmot::Math::Vector6d > >( response.stateVars,
                                                                                             stressSlot );
      Eigen::Map< NonlocalVector >
        Nn = this->stateLayout.template getAs< Eigen::Map< NonlocalVector > >( response.stateVars, nonlocalFieldSlot );
      double* baseState = this->baseMaterialState( response.stateVars );

      // A zeroed state vector means pristine, so an all-zero F_n is the identity here -- without this,
      // F_mid would be singular on the very first increment for a host that never calls initializeYourself().
      const Tensor33d Fn  = ( Fastor::norm( Fn_ref ) == 0.0 ) ? Tensor33d( Ident ) : Tensor33d( Fn_ref );
      const Tensor33d Fn1 = deformation.F;

      // --- Hughes-Winget kinematics, evaluated on the mid-step configuration ---
      const Tensor33d M    = inverse( Tensor33d( 0.5 * ( Fn + Fn1 ) ) );
      const Tensor33d dl   = Tensor33d( Fn1 - Fn ) % M;
      const Tensor33d dEps = 0.5 * ( dl + transpose( dl ) );
      const Tensor33d dOm  = 0.5 * ( dl - transpose( dl ) );

      const Tensor33d Ainv = inverse( Tensor33d( Ident - 0.5 * dOm ) );
      const Tensor33d dR   = Ainv % Tensor33d( Ident + 0.5 * dOm );

      // --- rotate the carried-over Cauchy stress forward ---
      const Marmot::Math::Vector6d sigmaNVoigt = sigmaN_n;
      const Tensor33d              sigmaN      = fromEigen( stressMatrixFromVoigt< 3 >( sigmaNVoigt ) );
      const Tensor33d              sigmaRot    = dR % sigmaN % transpose( dR );

      // --- assemble the increment for the wrapped material ---
      const int             nBase = baseMaterial->getNumberOfRequiredStateVars();
      std::vector< double > baseStateOld;
      if constexpr ( tangentMode == HughesWingetTangent::Exact ) {
        if ( computeTangent && nBase > 0 ) {
          baseStateOld.resize( nBase );
          std::memcpy( baseStateOld.data(), baseState, nBase * sizeof( double ) );
        }
      }

      const Marmot::Math::Vector6d sigmaRotVoigt = stressToVoigt< double >( toEigen( sigmaRot ) );
      const Marmot::Math::Vector6d dEpsVoigt     = voigtFromStrainMatrix< 3 >( toEigen( dEps ) );

      // The interface hands over the total nonlocal fields; the wrapped material wants the increment too.
      const typename BaseMaterial::increment inc = this->baseMaterialIncrement( deformation.N,
                                                                                Nn,
                                                                                dEpsVoigt,
                                                                                timeIncrement );

      /* Converted back out of the REFERENCE-volume convention, because the round trip has to close.
       * This interface documents both quantities as densities per unit reference volume and the
       * element stores and reloads them that way, while the wrapped small-strain material knows only
       * its own current volume -- which is why the outputs below are multiplied by J. Handing the
       * stored value straight back in skipped the inverse, so a material that ACCUMULATES into
       * either field would have had its history rescaled by J on every increment and compounded.
       *
       * Jn, not J: the value on entry was written at the end of the last increment as Jn * (the
       * wrapped material's own density), with Jn the determinant of the deformation gradient this
       * increment starts from. At the first increment Fn is the identity and this is exactly a
       * no-op. An inverted element is already meaningless here, so a non-positive Jn falls back to
       * passing the value through rather than dividing by it.
       */
      const double Jn                = determinant( Fn );
      const double seedEnergyDensity = Jn > 0.0 ? response.elasticEnergyDensity / Jn : response.elasticEnergyDensity;
      const double seedDissipation   = Jn > 0.0 ? response.dissipation / Jn : response.dissipation;

      typename BaseMaterial::response res;
      res.stress = sigmaRotVoigt;
      // BaseMaterial::response has no default member initialisers, so KLocal and c are otherwise
      // indeterminate on entry -- harmless for a material that only writes them, but a material that
      // reads or accumulates into them would observe garbage.
      res.KLocal.setZero();
      res.c.setZero();
      res.elasticEnergyDensity = seedEnergyDensity;
      res.dissipation          = seedDissipation;
      res.stateVars            = baseState;

      /* Every use of `tan` below sits inside `if ( computeTangent )`, so on the explicit path it is
       * computed and thrown away -- and for a return-mapping material like GCDP that is a full
       * consistent tangent per quadrature point per increment, on the very path this wrapper exists
       * to serve.
       */
      typename BaseMaterial::tangents tan;
      if ( computeTangent )
        baseMaterial->computeStress( res, tan, inc );
      else
        baseMaterial->computeStressExplicit( res, inc );

      const Marmot::Math::Vector6d sigmaNp1Voigt = res.stress;
      const Tensor33d              sigmaNp1      = fromEigen( stressMatrixFromVoigt< 3 >( sigmaNp1Voigt ) );

      // --- push the Cauchy stress forward to the Kirchhoff stress ---
      const double    J    = determinant( Fn1 );
      const Tensor33d Finv = inverse( Fn1 );

      response.tau = J * sigmaNp1;
      response.L   = res.KLocal;
      // The interface reports the radius, the wrapped material the interaction c = R^2.
      for ( int K = 0; K < nNonlocal; K++ )
        response.nonLocalRadius( K ) = std::sqrt( std::max( 0.0, res.c( K ) ) );
      // The wrapped material reports densities per unit current volume; the finite-strain consumers
      // integrate against the reference volume.
      response.elasticEnergyDensity = J * res.elasticEnergyDensity;
      response.dissipation          = J * res.dissipation;

      if ( computeTangent ) {
        const Tensor33d Pm = Ident - 0.5 * dl;
        const Tensor33d Mt = transpose( M );

        // d(dl)_ij/dF_kl = P_ik M_lj ; the transposed pattern gives d(dl)_ji/dF_kl
        const Tensor3333d dl_dF   = einsum< ik, jl, to_ijkl >( Pm, Mt );
        const Tensor3333d dlT_dF  = einsum< jk, il, to_ijkl >( Pm, Mt );
        const Tensor3333d dEps_dF = 0.5 * ( dl_dF + dlT_dF );
        const Tensor3333d dOm_dF  = 0.5 * ( dl_dF - dlT_dF );

        // d(sigmaRot)_ij/d(dOmega)_kl = 0.5 ( Ainv_ik G_lj + Ainv_jk G_li ),  G = (I + dR) sigmaN dR^T
        const Tensor33d   G           = Tensor33d( Ident + dR ) % sigmaN % transpose( dR );
        const Tensor33d   Gt          = transpose( G );
        const Tensor3333d dSigRot_dOm = 0.5 * ( einsum< ik, jl, to_ijkl >( Ainv, Gt ) +
                                                einsum< jk, il, to_ijkl >( Ainv, Gt ) );

        Tensor3333d dSigRot_dF = einsum< ijmn, mnKL, to_ijKL >( dSigRot_dOm, dOm_dF );

        if constexpr ( tangentMode == HughesWingetTangent::Exact ) {
          // Replace the implicit d(sigma^(n+1))/d(sigmaRot) = I assumption by the true operator.
          const Marmot::Math::Matrix6d S = computeStressSensitivity( sigmaRotVoigt,
                                                                     baseStateOld,
                                                                     inc,
                                                                     seedEnergyDensity,
                                                                     seedDissipation );
          dSigRot_dF                     = einsum< ijmn, mnKL, to_ijKL >( voigtToStiffness( S ), dSigRot_dF );
        }

        // voigtToStiffness scatters without a factor 1/2, which is exactly what a *tensorial*
        // strain derivative requires: summing over both (m,n) and (n,m) reproduces the engineering shear.
        const Tensor3333d C4      = voigtToStiffness( tan.dStressddStrain );
        const Tensor3333d dSig_dF = dSigRot_dF + einsum< ijmn, mnKL, to_ijKL >( C4, dEps_dF );

        const Tensor33d dJ_dF = J * transpose( Finv ); // dJ/dF_kl = J Finv_lk

        tangents.dTau_dF = einsum< ij, kl, to_ijkl >( sigmaNp1, dJ_dF ) + J * dSig_dF;

        for ( int K = 0; K < nNonlocal; K++ ) {
          // The nonlocal fields do not deform, so they reach the Kirchhoff stress through sigma alone.
          const Tensor33d dSig_dN = fromEigen(
            stressMatrixFromVoigt< 3 >( Marmot::Math::Vector6d( tan.dStressddK.col( K ) ) ) );
          tangents.dTau_dN[K] = J * dSig_dN;

          // dKLocal/d(dEps) is reported against the engineering-shear Voigt strain, so scattering it as a
          // stress -- shear entries placed unhalved on both off-diagonals -- is what makes the contraction
          // over all nine (i,j) reproduce the engineering shear exactly once.
          const Tensor33d dL_dEps = fromEigen(
            stressMatrixFromVoigt< 3 >( Marmot::Math::Vector6d( tan.dKLocalddStrain.row( K ).transpose() ) ) );
          tangents.dL_dF[K] = einsum< ij, ijKL, Fastor::OIndex< K_, L_ > >( dL_dEps, dEps_dF );
        }
        tangents.dL_dN = tan.dKLocalddK;
      }

      // --- commit ---
      sigmaN_n = sigmaNp1Voigt;
      Fn_ref   = Fn1;
      Nn       = deformation.N;
    }

    /**
     * @brief Forward-difference the complete update with respect to @f$ F @f$ and @f$ \bar{N} @f$.
     *
     * Every evaluation starts from a pristine copy of the incoming state and none of them commits: the
     * caller runs the unperturbed update afterwards. This is the verification oracle for the analytic
     * tangent, and it differentiates all four blocks.
     *
     * @param[in]  response      Constitutive response, read for its state pointer and energies only.
     * @param[out] tangents      Algorithmic moduli.
     * @param[in]  deformation   Current deformation gradient and nonlocal field.
     * @param[in]  timeIncrement Current (pseudo-)time and (pseudo-)time increment.
     */
    void computeNumericalTangent( const ConstitutiveResponse< 3 >& response,
                                  AlgorithmicModuli< 3 >&          tangents,
                                  const Deformation< 3 >&          deformation,
                                  const TimeIncrement&             timeIncrement ) const
    {
      using namespace Marmot::Math::Tensors::FastorTensors;

      const int             nTotal = this->getNumberOfRequiredStateVars();
      std::vector< double > stateOld( nTotal );
      std::memcpy( stateOld.data(), response.stateVars, nTotal * sizeof( double ) );

      std::vector< double >  scratch( nTotal );
      AlgorithmicModuli< 3 > unused;

      auto evaluate = [&]( const Tensor33d& F, const NonlocalVector& N ) -> std::pair< Tensor33d, NonlocalVector > {
        std::memcpy( scratch.data(), stateOld.data(), nTotal * sizeof( double ) );
        ConstitutiveResponse< 3 > perturbed;
        perturbed.elasticEnergyDensity = response.elasticEnergyDensity;
        perturbed.dissipation          = response.dissipation;
        perturbed.stateVars            = scratch.data();
        computeStressCore( perturbed, unused, Deformation< 3 >{ F, N }, timeIncrement, false );
        return { perturbed.tau, perturbed.L };
      };

      const auto [tau0, L0] = evaluate( deformation.F, deformation.N );

      for ( int k = 0; k < 3; ++k )
        for ( int l = 0; l < 3; ++l ) {
          Tensor33d    F = deformation.F;
          const double h = std::sqrt( std::numeric_limits< double >::epsilon() ) *
                           std::max( 1.0, std::abs( F( k, l ) ) );
          F( k, l ) += h;

          const auto [tauPerturbed, LPerturbed] = evaluate( F, deformation.N );

          for ( int i = 0; i < 3; ++i )
            for ( int j = 0; j < 3; ++j )
              tangents.dTau_dF( i, j, k, l ) = ( tauPerturbed( i, j ) - tau0( i, j ) ) / h;

          for ( int K = 0; K < nNonlocal; ++K )
            tangents.dL_dF[K]( k, l ) = ( LPerturbed( K ) - L0( K ) ) / h;
        }

      for ( int M = 0; M < nNonlocal; ++M ) {
        const double hN = std::sqrt( std::numeric_limits< double >::epsilon() ) *
                          std::max( 1.0, std::abs( deformation.N( M ) ) );

        NonlocalVector N = deformation.N;
        N( M ) += hN;

        const auto [tauPerturbed, LPerturbed] = evaluate( deformation.F, N );

        for ( int i = 0; i < 3; ++i )
          for ( int j = 0; j < 3; ++j )
            tangents.dTau_dN[M]( i, j ) = ( tauPerturbed( i, j ) - tau0( i, j ) ) / hN;

        for ( int K = 0; K < nNonlocal; ++K )
          tangents.dL_dN( K, M ) = ( LPerturbed( K ) - L0( K ) ) / hN;
      }
    }
  };

} // namespace Marmot::Materials
