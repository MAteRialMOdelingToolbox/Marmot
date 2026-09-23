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
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include "Marmot/MarmotMaterialHypoElastic.h"
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
   * @brief Selects how the algorithmic tangent @f$ \partial\boldsymbol{\tau}/\partial\boldsymbol{F} @f$
   *        of @ref HughesWingetWrapper is evaluated.
   */
  enum class HughesWingetTangent {
    /**
     * Closed-form linearisation of the Hughes-Winget update, assuming
     * @f$ \partial\boldsymbol{\sigma}^{(n+1)}/\partial\boldsymbol{\sigma}_{\text{rot}} = \boldsymbol{I} @f$.
     * Exact while the wrapped material responds elastically. One wrapped-material evaluation.
     */
    Analytic,
    /**
     * As HughesWingetTangent::Analytic, but the operator
     * @f$ \boldsymbol{S} = \partial\boldsymbol{\sigma}^{(n+1)}/\partial\boldsymbol{\sigma}_{\text{rot}} @f$
     * is recovered by forward-differencing the wrapped material with respect to its six incoming stress
     * components. Fully consistent. Seven wrapped-material evaluations.
     */
    Exact,
    /**
     * The complete update is forward-differenced with respect to the nine components of
     * @f$ \boldsymbol{F} @f$. Fully consistent and model-agnostic, but the most expensive option
     * (eleven wrapped-material evaluations). Intended for verification and debugging.
     */
    Numerical
  };

  /**
   * @class HughesWingetWrapper
   * @brief Decorator presenting a small-strain (hypoelastic) material as a finite-strain material.
   *
   * @tparam BaseMaterialType A concrete class derived from MarmotMaterialHypoElastic.
   * @tparam tangentMode      How the algorithmic tangent is evaluated, see Marmot::Materials::HughesWingetTangent.
   *
   * The wrapped material is driven by the objective strain increment of the Hughes-Winget algorithm,
   * evaluated on the mid-step configuration,
   * @f[
   *   \Delta\boldsymbol{l} = \left(\boldsymbol{F}^{(n+1)} - \boldsymbol{F}^{(n)}\right)
   *                          \boldsymbol{F}_{\text{mid}}^{-1}, \qquad
   *   \boldsymbol{F}_{\text{mid}} = \tfrac{1}{2}\left(\boldsymbol{F}^{(n)} + \boldsymbol{F}^{(n+1)}\right),
   * @f]
   * split into its symmetric and skew parts @f$ \Delta\boldsymbol{\varepsilon} @f$ and
   * @f$ \Delta\boldsymbol{\Omega} @f$. The incremental rotation follows from the Cayley transform
   * @f[
   *   \Delta\boldsymbol{R} = \left(\boldsymbol{I} - \tfrac{1}{2}\Delta\boldsymbol{\Omega}\right)^{-1}
   *                          \left(\boldsymbol{I} + \tfrac{1}{2}\Delta\boldsymbol{\Omega}\right),
   * @f]
   * which is orthogonal to machine precision. The stress carried over from the previous increment is
   * rotated forward, @f$ \boldsymbol{\sigma}_{\text{rot}} = \Delta\boldsymbol{R}\,\boldsymbol{\sigma}^{(n)}
   * \Delta\boldsymbol{R}^{T} @f$, the wrapped material is evaluated, and the resulting Cauchy stress is
   * pushed to the Kirchhoff stress @f$ \boldsymbol{\tau} = J\,\boldsymbol{\sigma}^{(n+1)} @f$ expected by
   * the finite-strain consumers.
   *
   * @note **Only the stress is rotated.** Tensor-valued internal variables of the wrapped material
   *       (kinematic hardening back stresses, plastic strain tensors, anisotropic damage tensors) are
   *       passed through untouched and are therefore *not* objective under large incremental rotations.
   *       Models whose internal state is scalar or isotropic-invariant are unaffected.
   *
   * @note The characteristic element length is **not** derived from any geometry here. It must be
   *       assigned explicitly via @ref setCharacteristicElementLength; until then it is NaN, so a
   *       wrapped material that relies on it fails loudly instead of silently regularising incorrectly.
   */
  template < typename BaseMaterialType, HughesWingetTangent tangentMode = HughesWingetTangent::Analytic >
  class HughesWingetWrapper : public MarmotMaterialFiniteStrain {

    static_assert( std::is_base_of_v< MarmotMaterialHypoElastic, BaseMaterialType >,
                   "HughesWingetWrapper can only wrap materials derived from MarmotMaterialHypoElastic." );

  protected:
    /// The wrapped small-strain material. Held as the base interface: concrete materials are free to
    /// re-declare their overrides as protected, which would make them unreachable via BaseMaterialType.
    std::unique_ptr< MarmotMaterialHypoElastic > baseMaterial;

    /// Name of the state layout slot holding the deformation gradient of the last accepted increment.
    static constexpr const char* deformationGradientSlot = "HughesWinget_F_n";
    /// Name of the state layout slot holding the Cauchy stress of the last accepted increment.
    static constexpr const char* stressSlot = "HughesWinget_sigma_n";
    /// Name of the state layout slot holding the state variables of the wrapped material.
    static constexpr const char* baseMaterialSlot = "materialstate";

  public:
    /**
     * @brief Construct the wrapper and the wrapped material.
     * @param[in] matProperties_       Material property values; forwarded **unchanged** to the wrapped material.
     * @param[in] nMaterialProperties_ Number of material property values.
     * @param[in] materialNumber_      Unique identifier for this material instance.
     */
    HughesWingetWrapper( const double* matProperties_, int nMaterialProperties_, int materialNumber_ )
      : MarmotMaterialFiniteStrain( matProperties_, nMaterialProperties_, materialNumber_ ),
        baseMaterial( std::make_unique< BaseMaterialType >( matProperties_, nMaterialProperties_, materialNumber_ ) )
    {
      // Deliberately poisoned: an unassigned characteristic length must not silently behave like a valid
      // one. Materials that never read it are unaffected.
      baseMaterial->setCharacteristicElementLength( std::numeric_limits< double >::quiet_NaN() );

      initializeStateLayout();
    }

    virtual ~HughesWingetWrapper() = default;

    /**
     * @brief Register the state layout: the carried-over deformation gradient and Cauchy stress,
     *        followed by the state variables of the wrapped material.
     */
    void initializeStateLayout()
    {
      this->stateLayout.add( deformationGradientSlot, 9 );
      this->stateLayout.add( stressSlot, 6 );
      this->stateLayout.add( baseMaterialSlot, baseMaterial->getNumberOfRequiredStateVars() );
      this->stateLayout.finalize();
    }

    /**
     * @brief Forward the characteristic element length to the wrapped material.
     * @param[in] length Characteristic element length at the considered evaluation point.
     */
    void setCharacteristicElementLength( double length ) override
    {
      baseMaterial->setCharacteristicElementLength( length );
    }

    /**
     * @brief Resolve a state variable by name, falling back to the wrapped material.
     *
     * Without this, the wrapped material's own state -- CDP's @c omega, a plasticity model's
     * hardening variable -- would be unreachable through the wrapper, since the wrapper's layout only
     * knows its own three slots.
     *
     * @param[in] stateName Name of the state variable.
     * @param[in] stateVars Pointer to the state variable array of this wrapper.
     * @return A view of the requested state variable.
     */
    StateView getStateView( const std::string& stateName, double* stateVars ) const override
    {
      if ( stateName == deformationGradientSlot || stateName == stressSlot || stateName == baseMaterialSlot )
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
     * @brief Initialise the state: @f$ \boldsymbol{F}^{(n)} = \boldsymbol{I} @f$,
     *        @f$ \boldsymbol{\sigma}^{(n)} = \boldsymbol{0} @f$, plus the wrapped material's own state.
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

      baseMaterial->initializeYourself( this->stateLayout.getPtr( stateVars, baseMaterialSlot ),
                                        baseMaterial->getNumberOfRequiredStateVars() );
    }

    /**
     * @brief Objective stress update of the wrapped small-strain material.
     * @param[in,out] response      Constitutive response; carries the Kirchhoff stress and state variables.
     * @param[out]    tangents      Algorithmic moduli @f$ \partial\tau_{ij}/\partial F_{kl} @f$.
     * @param[in]     deformation   Current deformation gradient.
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
     * @param[in]     deformation   Current deformation gradient.
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
     *
     * @param[in] sigmaRotVoigt  Rotated stress that was handed to the wrapped material.
     * @param[in] sigmaNp1Voigt  Unperturbed updated stress returned by the wrapped material.
     * @param[in] baseStateOld   Wrapped material state *before* the unperturbed evaluation.
     * @param[in] dStrain        Strain increment handed to the wrapped material.
     * @param[in] timeInfo       Time information handed to the wrapped material.
     * @return The sensitivity in Voigt form, with its three shear **columns** halved so that the
     *         result may be contracted as a full fourth-order tensor without double counting the
     *         off-diagonal stress components.
     */
    Marmot::Matrix6d computeStressSensitivity( const Marmot::Vector6d&                    sigmaRotVoigt,
                                               const Marmot::Vector6d&                    sigmaNp1Voigt,
                                               const std::vector< double >&               baseStateOld,
                                               const Marmot::Vector6d&                    dStrain,
                                               const MarmotMaterialHypoElastic::timeInfo& timeInfo ) const
    {
      const int nBase = baseMaterial->getNumberOfRequiredStateVars();

      const double scale = std::max( 1.0, sigmaRotVoigt.cwiseAbs().maxCoeff() );
      const double h     = std::sqrt( std::numeric_limits< double >::epsilon() ) * scale;

      Marmot::Matrix6d S = Marmot::Matrix6d::Zero();

      std::vector< double > scratch( nBase );
      for ( int j = 0; j < 6; ++j ) {
        if ( nBase > 0 )
          std::memcpy( scratch.data(), baseStateOld.data(), nBase * sizeof( double ) );

        MarmotMaterialHypoElastic::state3D perturbed;
        perturbed.stress = sigmaRotVoigt;
        perturbed.stress( j ) += h;
        perturbed.stateVars = nBase > 0 ? scratch.data() : nullptr;

        Marmot::Matrix6d dummyTangent = Marmot::Matrix6d::Zero();
        baseMaterial->computeStress( perturbed, dummyTangent, dStrain, timeInfo );

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
     * @param[in]     deformation    Current deformation gradient.
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
      double*      baseState = this->stateLayout.getPtr( response.stateVars, baseMaterialSlot );

      // Hosts are not obliged to call initializeYourself(): DisplacementFiniteStrainULElement only does
      // so on an explicit "initialize material" initial condition. Marmot's convention is that a zeroed
      // state vector means pristine, so an all-zero F_n is the identity here -- without this, F_mid would
      // be singular on the very first increment and the wrapped material would see a garbage increment.
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

      // --- evaluate the wrapped small-strain material ---
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

      MarmotMaterialHypoElastic::state3D state;
      state.stress               = sigmaRotVoigt;
      state.elasticEnergyDensity = response.elasticEnergyDensity;
      state.dissipation          = response.dissipation;
      state.stateVars            = baseState;

      const MarmotMaterialHypoElastic::timeInfo timeInfo{ timeIncrement.time, timeIncrement.dT };

      Marmot::Matrix6d C = Marmot::Matrix6d::Zero();
      baseMaterial->computeStress( state, C, dEpsVoigt, timeInfo );

      const Marmot::Vector6d sigmaNp1Voigt = state.stress;
      const Tensor33d        sigmaNp1      = fromEigen( stressMatrixFromVoigt< 3 >( sigmaNp1Voigt ) );

      // --- push the Cauchy stress forward to the Kirchhoff stress ---
      const double    J    = determinant( Fn1 );
      const Tensor33d Finv = inverse( Fn1 );

      response.tau = J * sigmaNp1;
      // The wrapped material reports densities per unit current volume; the finite-strain consumers
      // integrate against the reference volume.
      response.elasticEnergyDensity = J * state.elasticEnergyDensity;
      response.dissipation          = J * state.dissipation;

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
          const Marmot::Matrix6d S = computeStressSensitivity( sigmaRotVoigt,
                                                               sigmaNp1Voigt,
                                                               baseStateOld,
                                                               dEpsVoigt,
                                                               timeInfo );
          dSigRot_dF               = einsum< ijmn, mnKL, to_ijKL >( voigtToStiffnessFastor( S ), dSigRot_dF );
        }

        // voigtToStiffnessFastor scatters without a factor 1/2, which is exactly what a *tensorial*
        // strain derivative requires: summing over both (m,n) and (n,m) reproduces the engineering shear.
        const Tensor3333d C4      = voigtToStiffnessFastor( C );
        const Tensor3333d dSig_dF = dSigRot_dF + einsum< ijmn, mnKL, to_ijKL >( C4, dEps_dF );

        const Tensor33d dJ_dF = J * transpose( Finv ); // dJ/dF_kl = J Finv_lk

        tangents.dTau_dF = einsum< ij, kl, to_ijkl >( sigmaNp1, dJ_dF ) + J * dSig_dF;
      }

      // --- commit ---
      Eigen::Map< Marmot::Vector6d > sigmaNOut( sigmaNPtr );
      sigmaNOut = sigmaNp1Voigt;
      Fn_ref    = Fn1;
    }

    /**
     * @brief Forward-difference the complete update with respect to the nine components of @f$ F @f$.
     *
     * Every evaluation starts from a pristine copy of the incoming state, and none of them is allowed
     * to commit: the caller runs the unperturbed update afterwards.
     *
     * @param[in]  response      Constitutive response, read for its state pointer only.
     * @param[out] tangents      Algorithmic moduli.
     * @param[in]  deformation   Current deformation gradient.
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

      auto evaluate = [&]( const Tensor33d& F ) -> Tensor33d {
        std::memcpy( scratch.data(), stateOld.data(), nTotal * sizeof( double ) );
        ConstitutiveResponse< 3 > perturbed;
        // Seed from the incoming (accumulated) response rather than the default-constructed zero:
        // computeStressCore feeds these into the wrapped material's state as its starting energy /
        // dissipation, and a material whose response depends on that history (not merely writing to
        // it) would otherwise be perturbed from a fictitious zero baseline instead of its real one.
        perturbed.elasticEnergyDensity = response.elasticEnergyDensity;
        perturbed.dissipation          = response.dissipation;
        perturbed.stateVars            = scratch.data();
        computeStressCore( perturbed, unused, Deformation< 3 >{ F }, timeIncrement, false );
        return perturbed.tau;
      };

      const Tensor33d tau0 = evaluate( deformation.F );

      for ( int k = 0; k < 3; ++k )
        for ( int l = 0; l < 3; ++l ) {
          Tensor33d    F = deformation.F;
          const double h = std::sqrt( std::numeric_limits< double >::epsilon() ) *
                           std::max( 1.0, std::abs( F( k, l ) ) );
          F( k, l ) += h;

          const Tensor33d tauPerturbed = evaluate( F );

          for ( int i = 0; i < 3; ++i )
            for ( int j = 0; j < 3; ++j )
              tangents.dTau_dF( i, j, k, l ) = ( tauPerturbed( i, j ) - tau0( i, j ) ) / h;
        }
    }
  };

} // namespace Marmot::Materials
