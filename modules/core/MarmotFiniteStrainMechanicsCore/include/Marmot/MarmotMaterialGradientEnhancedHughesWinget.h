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
#include "Marmot/MarmotMaterialGeneralGradientEnhancedHypoElastic.h"
#include "Marmot/MarmotMaterialGradientEnhancedFiniteStrain.h"
#include "Marmot/MarmotMaterialHughesWinget.h"
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
   * @tparam BaseMaterialType A concrete class derived from MarmotMaterialGeneralGradientEnhancedHypoElastic<1>.
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
   * On top of that it carries the nonlocal field across the two interfaces:
   *
   * | this interface reports | from the wrapped material |
   * |---|---|
   * | @c response.L | @c KLocal(0), a total, which is what both sides mean |
   * | @c response.nonLocalRadius | @f$ \sqrt{c(0)} @f$ |
   * | @c tangents.dTau_dN | @f$ J\,\partial\boldsymbol{\sigma}/\partial\bar{N} @f$ |
   * | @c tangents.dL_dF | @f$ \partial K^{\text{local}}/\partial\Delta\boldsymbol{\varepsilon} @f$, chained through the
   * Hughes-Winget kinematics | | @c tangents.dL_dN | @c dKLocalddK(0,0) |
   *
   * @note The wrapped material wants the nonlocal field **and its increment**, while this interface hands
   * over only the total. The increment is therefore formed here, against a value of the last accepted
   * increment carried in this wrapper's own state (@c HughesWinget_N_n) -- one slot more than the local
   * sibling needs.
   *
   * @note **The wrapped material must be incremental in stress.** It has to update the Cauchy stress it
   * is handed, not recompute one from a stored strain: the wrapper's whole mechanism is to hand over the
   * forward-rotated stress of the last increment, and a material that ignores that argument discards the
   * rotation along with it. `GCDPModel` qualifies -- it maps `res.stress` and updates it in place
   * (`GCDP.cpp:57`), and its entire state is four scalars. **`AT2PhaseField` does not**: it carries a
   * six-component `strain` state (`AT2PhaseField.h:60`) and returns @f$ g(\varphi)\,\mathbb{C}:
   * \boldsymbol{\varepsilon} @f$ from it (`AT2PhaseField.cpp:71`), so its true
   * @f$ \partial\boldsymbol{\sigma}^{(n+1)}/\partial\boldsymbol{\sigma}_{\text{rot}} @f$ is zero
   * rather than the identity. Measured: wrapped, a rigid rotation on top of an anisotropic stretch leaves
   * its Kirchhoff stress bit-for-bit **unrotated**. It is therefore deliberately not registered with this
   * wrapper, and a total-strain-based model must not be either.
   *
   * @note **Only the stress is rotated.** Tensor-valued internal variables of the wrapped material are
   * passed through untouched and are therefore *not* objective under large incremental rotations. A model
   * whose internal state is entirely scalar -- GCDP's `alphaP`, `alphaD`, `omega`, `I1p` -- is unaffected.
   *
   * @note The interaction @f$ c @f$ is reported through its square root, and its derivative
   * @f$ \partial c/\partial\bar{N} @f$ **cannot be forwarded at all**: the finite-strain interface has no
   * slot for it. Harmless for a material whose @f$ c @f$ is constant (GCDP, AT2PhaseField both are); a
   * material built on MarmotDecreasingInteractions would keep a correct residual but lose that term of
   * its consistent tangent.
   *
   * @note The nonlocal balance is formulated in the material configuration by the consumers of this
   * interface, so @f$ c @f$ is handed over as the material constant the wrapped model reports, neither
   * pushed forward nor weighted by @f$ J @f$.
   */
  template < typename BaseMaterialType, HughesWingetTangent tangentMode = HughesWingetTangent::Analytic >
  class GradientEnhancedHughesWingetWrapper : public MarmotMaterialGradientEnhancedFiniteStrain {

    static_assert( std::is_base_of_v< MarmotMaterialGeneralGradientEnhancedHypoElastic< 1 >, BaseMaterialType >,
                   "GradientEnhancedHughesWingetWrapper can only wrap materials derived from "
                   "MarmotMaterialGeneralGradientEnhancedHypoElastic<1>." );

    using BaseMaterial = MarmotMaterialGeneralGradientEnhancedHypoElastic< 1 >;

  protected:
    /// The wrapped small-strain material, held as the base interface.
    std::unique_ptr< BaseMaterial > baseMaterial;

    /// Name of the state layout slot holding the deformation gradient of the last accepted increment.
    static constexpr const char* deformationGradientSlot = "HughesWinget_F_n";
    /// Name of the state layout slot holding the Cauchy stress of the last accepted increment.
    static constexpr const char* stressSlot = "HughesWinget_sigma_n";
    /// Name of the state layout slot holding the nonlocal field of the last accepted increment.
    static constexpr const char* nonlocalFieldSlot = "HughesWinget_N_n";
    /// Name of the state layout slot holding the state variables of the wrapped material.
    static constexpr const char* baseMaterialSlot = "materialstate";

  public:
    /**
     * @brief Construct the wrapper and the wrapped material.
     * @param[in] matProperties_       Material property values; forwarded **unchanged** to the wrapped material.
     * @param[in] nMaterialProperties_ Number of material property values.
     * @param[in] materialNumber_      Unique identifier for this material instance.
     */
    GradientEnhancedHughesWingetWrapper( const double* matProperties_, int nMaterialProperties_, int materialNumber_ )
      : MarmotMaterialGradientEnhancedFiniteStrain( matProperties_, nMaterialProperties_, materialNumber_ ),
        baseMaterial( std::make_unique< BaseMaterialType >( matProperties_, nMaterialProperties_, materialNumber_ ) )
    {
      initializeStateLayout();
    }

    virtual ~GradientEnhancedHughesWingetWrapper() = default;

    /**
     * @brief Register the state layout: the carried-over deformation gradient, Cauchy stress and nonlocal
     *        field, followed by the state variables of the wrapped material.
     */
    void initializeStateLayout()
    {
      this->stateLayout.add( deformationGradientSlot, 9 );
      this->stateLayout.add( stressSlot, 6 );
      this->stateLayout.add( nonlocalFieldSlot, 1 );
      this->stateLayout.add( baseMaterialSlot, baseMaterial->getNumberOfRequiredStateVars() );
      this->stateLayout.finalize();
    }

    /**
     * @brief Resolve a state variable by name, falling back to the wrapped material.
     * @param[in] stateName Name of the state variable.
     * @param[in] stateVars Pointer to the state variable array of this wrapper.
     * @return A view of the requested state variable.
     */
    StateView getStateView( const std::string& stateName, double* stateVars ) const override
    {
      if ( stateName == deformationGradientSlot || stateName == stressSlot || stateName == nonlocalFieldSlot ||
           stateName == baseMaterialSlot )
        return this->stateLayout.getStateView( stateVars, stateName );

      return baseMaterial->getStateView( stateName, this->stateLayout.getPtr( stateVars, baseMaterialSlot ) );
    }

    /**
     * @brief Mass density of the wrapped material.
     * @param[in] stateVars Pointer to the state variable array of this wrapper.
     * @return Mass density.
     */
    double getDensity( const double* stateVars ) const override
    {
      return baseMaterial->getDensity(
        this->stateLayout.getPtr( const_cast< double* >( stateVars ), baseMaterialSlot ) );
    }

    /**
     * @brief Nonlocal viscosity of the wrapped small-strain material.
     * @param[in] stateVars Pointer to the state variable array of this wrapper.
     * @return Nonlocal viscosity @f$ \eta @f$, the first entry (and only, since @c BaseMaterialType is
     *         templated on a single nonlocal variable) of the wrapped material's own
     *         @c getNonlocalViscosity().
     */
    double getNonlocalViscosity( const double* stateVars ) const override
    {
      return baseMaterial->getNonlocalViscosity(
        this->stateLayout.getPtr( const_cast< double* >( stateVars ), baseMaterialSlot ) )[0];
    }

    /**
     * @brief Nonlocal micro-inertia of the wrapped small-strain material.
     * @param[in] stateVars Pointer to the state variable array of this wrapper.
     * @return Nonlocal micro-inertia @f$ m_k @f$, the first entry of the wrapped material's own
     *         @c getNonlocalMicroInertia().
     */
    double getNonlocalMicroInertia( const double* stateVars ) const override
    {
      return baseMaterial->getNonlocalMicroInertia(
        this->stateLayout.getPtr( const_cast< double* >( stateVars ), baseMaterialSlot ) )[0];
    }

    /**
     * @brief Initialise the state: @f$ \boldsymbol{F}^{(n)} = \boldsymbol{I} @f$,
     *        @f$ \boldsymbol{\sigma}^{(n)} = \boldsymbol{0} @f$, @f$ \bar{N}^{(n)} = 0 @f$, plus the
     *        wrapped material's own state.
     * @param[in,out] stateVars  Pointer to the state variable array.
     * @param[in]     nStateVars Number of state variables.
     */
    void initializeYourself( double* stateVars, int nStateVars ) override
    {
      using namespace Marmot::FastorStandardTensors;

      TensorMap33d Fn = this->stateLayout.getAs< TensorMap33d >( stateVars, deformationGradientSlot );
      std::memcpy( Fn.data(), Spatial3D::I.data(), 9 * sizeof( double ) );

      double* sigmaN = this->stateLayout.getPtr( stateVars, stressSlot );
      std::fill( sigmaN, sigmaN + 6, 0.0 );

      *this->stateLayout.getPtr( stateVars, nonlocalFieldSlot ) = 0.0;

      baseMaterial->initializeYourself( this->stateLayout.getPtr( stateVars, baseMaterialSlot ),
                                        baseMaterial->getNumberOfRequiredStateVars() );
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
    /// @brief Convert an Eigen 3x3 matrix into a Fastor tensor.
    static Marmot::FastorStandardTensors::Tensor33d fromEigen( const Eigen::Matrix3d& m )
    {
      Marmot::FastorStandardTensors::Tensor33d t;
      Marmot::mapEigenToFastor( t ) = m;
      return t;
    }

    /// @brief Convert a Fastor 3x3 tensor into an Eigen matrix.
    static Eigen::Matrix3d toEigen( const Marmot::FastorStandardTensors::Tensor33d& t )
    {
      return Eigen::Matrix3d( Marmot::mapEigenToFastor( t ) );
    }

    /**
     * @brief Recover @f$ \boldsymbol{S} = \partial\boldsymbol{\sigma}^{(n+1)} /
     *        \partial\boldsymbol{\sigma}_{\text{rot}} @f$ by forward differences.
     * @param[in] sigmaRotVoigt Rotated stress that was handed to the wrapped material.
     * @param[in] sigmaNp1Voigt Unperturbed updated stress returned by the wrapped material.
     * @param[in] baseStateOld  Wrapped material state *before* the unperturbed evaluation.
     * @param[in] inc           Increment handed to the wrapped material.
     * @return The sensitivity in Voigt form, with its three shear **columns** halved so that the result
     *         may be contracted as a full fourth-order tensor without double counting.
     */
    Marmot::Matrix6d computeStressSensitivity( const Marmot::Vector6d&        sigmaRotVoigt,
                                               const Marmot::Vector6d&        sigmaNp1Voigt,
                                               const std::vector< double >&   baseStateOld,
                                               const BaseMaterial::increment& inc ) const
    {
      const int nBase = baseMaterial->getNumberOfRequiredStateVars();

      const double scale = std::max( 1.0, sigmaRotVoigt.cwiseAbs().maxCoeff() );
      const double h     = std::sqrt( std::numeric_limits< double >::epsilon() ) * scale;

      Marmot::Matrix6d S = Marmot::Matrix6d::Zero();

      std::vector< double > scratch( nBase );
      for ( int j = 0; j < 6; ++j ) {
        if ( nBase > 0 )
          std::memcpy( scratch.data(), baseStateOld.data(), nBase * sizeof( double ) );

        BaseMaterial::response perturbed;
        perturbed.stress = sigmaRotVoigt;
        perturbed.stress( j ) += h;
        perturbed.stateVars = nBase > 0 ? scratch.data() : nullptr;

        BaseMaterial::tangents dummyTangents;
        baseMaterial->computeStress( perturbed, dummyTangents, inc );

        S.col( j ) = ( perturbed.stress - sigmaNp1Voigt ) / h;
      }

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
      using namespace Marmot::FastorStandardTensors;
      using namespace Marmot::FastorIndices;
      using namespace Marmot::ContinuumMechanics::VoigtNotation;

      const Tensor33d& Ident = Spatial3D::I;

      TensorMap33d Fn_ref    = this->stateLayout.getAs< TensorMap33d >( response.stateVars, deformationGradientSlot );
      double*      sigmaNPtr = this->stateLayout.getPtr( response.stateVars, stressSlot );
      double*      nPtr      = this->stateLayout.getPtr( response.stateVars, nonlocalFieldSlot );
      double*      baseState = this->stateLayout.getPtr( response.stateVars, baseMaterialSlot );

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
      Eigen::Map< const Marmot::Vector6d > sigmaNVoigtMap( sigmaNPtr );
      const Marmot::Vector6d               sigmaNVoigt = sigmaNVoigtMap;
      const Tensor33d                      sigmaN      = fromEigen( stressMatrixFromVoigt< 3 >( sigmaNVoigt ) );
      const Tensor33d                      sigmaRot    = dR % sigmaN % transpose( dR );

      // --- assemble the increment for the wrapped material ---
      const int             nBase = baseMaterial->getNumberOfRequiredStateVars();
      std::vector< double > baseStateOld;
      if constexpr ( tangentMode == HughesWingetTangent::Exact ) {
        if ( computeTangent && nBase > 0 ) {
          baseStateOld.resize( nBase );
          std::memcpy( baseStateOld.data(), baseState, nBase * sizeof( double ) );
        }
      }

      const Marmot::Vector6d sigmaRotVoigt = stressToVoigt< double >( toEigen( sigmaRot ) );
      const Marmot::Vector6d dEpsVoigt     = voigtFromStrainMatrix< 3 >( toEigen( dEps ) );

      // The interface hands over the total nonlocal field; the wrapped material wants the increment too.
      const double nNew = deformation.N;
      const double dN   = nNew - *nPtr;

      BaseMaterial::increment inc;
      inc.dStrain = dEpsVoigt;
      inc.K( 0 )  = nNew;
      inc.dK( 0 ) = dN;
      inc.time    = timeIncrement.time;
      inc.dT      = timeIncrement.dT;

      BaseMaterial::response res;
      res.stress               = sigmaRotVoigt;
      res.elasticEnergyDensity = response.elasticEnergyDensity;
      res.dissipation          = response.dissipation;
      res.stateVars            = baseState;

      BaseMaterial::tangents tan;
      baseMaterial->computeStress( res, tan, inc );

      const Marmot::Vector6d sigmaNp1Voigt = res.stress;
      const Tensor33d        sigmaNp1      = fromEigen( stressMatrixFromVoigt< 3 >( sigmaNp1Voigt ) );

      // --- push the Cauchy stress forward to the Kirchhoff stress ---
      const double    J    = determinant( Fn1 );
      const Tensor33d Finv = inverse( Fn1 );

      response.tau = J * sigmaNp1;
      response.L   = res.KLocal( 0 );
      // The interface reports the radius, the wrapped material the interaction c = R^2.
      response.nonLocalRadius = std::sqrt( std::max( 0.0, res.c( 0 ) ) );
      // The wrapped material reports densities per unit current volume; the finite-strain consumers
      // integrate against the reference volume.
      response.elasticEnergyDensity = J * res.elasticEnergyDensity;
      response.dissipation          = J * res.dissipation;

      if ( computeTangent ) {
        const Tensor33d P  = Ident - 0.5 * dl;
        const Tensor33d Mt = transpose( M );

        // d(dl)_ij/dF_kl = P_ik M_lj ; the transposed pattern gives d(dl)_ji/dF_kl
        const Tensor3333d dl_dF   = einsum< ik, jl, to_ijkl >( P, Mt );
        const Tensor3333d dlT_dF  = einsum< jk, il, to_ijkl >( P, Mt );
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
          const Marmot::Matrix6d S = computeStressSensitivity( sigmaRotVoigt, sigmaNp1Voigt, baseStateOld, inc );
          dSigRot_dF               = einsum< ijmn, mnKL, to_ijKL >( voigtToStiffnessFastor( S ), dSigRot_dF );
        }

        // voigtToStiffnessFastor scatters without a factor 1/2, which is exactly what a *tensorial*
        // strain derivative requires: summing over both (m,n) and (n,m) reproduces the engineering shear.
        const Tensor3333d C4      = voigtToStiffnessFastor( tan.dStressddStrain );
        const Tensor3333d dSig_dF = dSigRot_dF + einsum< ijmn, mnKL, to_ijKL >( C4, dEps_dF );

        const Tensor33d dJ_dF = J * transpose( Finv ); // dJ/dF_kl = J Finv_lk

        tangents.dTau_dF = einsum< ij, kl, to_ijkl >( sigmaNp1, dJ_dF ) + J * dSig_dF;

        // The nonlocal field does not deform, so it reaches the Kirchhoff stress through sigma alone.
        const Tensor33d dSig_dN = fromEigen(
          stressMatrixFromVoigt< 3 >( Marmot::Vector6d( tan.dStressddK.col( 0 ) ) ) );
        tangents.dTau_dN = J * dSig_dN;

        // dKLocal/d(dEps) is reported against the engineering-shear Voigt strain, so scattering it as a
        // stress -- shear entries placed unhalved on both off-diagonals -- is what makes the contraction
        // over all nine (i,j) reproduce the engineering shear exactly once.
        const Tensor33d dL_dEps = fromEigen(
          stressMatrixFromVoigt< 3 >( Marmot::Vector6d( tan.dKLocalddStrain.row( 0 ).transpose() ) ) );
        tangents.dL_dF = einsum< ij, ijKL, Fastor::OIndex< K_, L_ > >( dL_dEps, dEps_dF );
        tangents.dL_dN = tan.dKLocalddK( 0, 0 );
      }

      // --- commit ---
      Eigen::Map< Marmot::Vector6d > sigmaNOut( sigmaNPtr );
      sigmaNOut = sigmaNp1Voigt;
      Fn_ref    = Fn1;
      *nPtr     = nNew;
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
      using namespace Marmot::FastorStandardTensors;

      const int             nTotal = this->getNumberOfRequiredStateVars();
      std::vector< double > stateOld( nTotal );
      std::memcpy( stateOld.data(), response.stateVars, nTotal * sizeof( double ) );

      std::vector< double >  scratch( nTotal );
      AlgorithmicModuli< 3 > unused;

      auto evaluate = [&]( const Tensor33d& F, double N ) -> std::pair< Tensor33d, double > {
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

          tangents.dL_dF( k, l ) = ( LPerturbed - L0 ) / h;
        }

      {
        const double hN = std::sqrt( std::numeric_limits< double >::epsilon() ) *
                          std::max( 1.0, std::abs( deformation.N ) );

        const auto [tauPerturbed, LPerturbed] = evaluate( deformation.F, deformation.N + hN );

        for ( int i = 0; i < 3; ++i )
          for ( int j = 0; j < 3; ++j )
            tangents.dTau_dN( i, j ) = ( tauPerturbed( i, j ) - tau0( i, j ) ) / hN;

        tangents.dL_dN = ( LPerturbed - L0 ) / hN;
      }
    }
  };

} // namespace Marmot::Materials
