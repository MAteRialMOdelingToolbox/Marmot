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
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"
#include <Eigen/Eigenvalues>
#include <cmath>
#include <functional>

namespace Marmot::Materials {

  /**
   * @class GradientEnhancedCorotationalWrapper
   * @brief Decorator presenting a small-strain gradient-enhanced material as a finite-strain one, in the frame of the
   * polar rotation.
   *
   * @tparam BaseMaterialType A concrete class derived from Marmot::Materials::MaterialGradientEnhancedHypoElastic<n>;
   *                          the wrapper has the same n nonlocal variables.
   *
   * The deformation gradient is decomposed into the polar rotation and the right stretch,
   * @f$ \boldsymbol{F} = \boldsymbol{R}\,\boldsymbol{U} @f$, @f$ \boldsymbol{U} = \boldsymbol{C}^{1/2} @f$,
   * @f$ \boldsymbol{C} = \boldsymbol{F}^T\boldsymbol{F} @f$. The wrapped material lives in the frame that rotates with
   * @f$ \boldsymbol{R} @f$ (the frame of the material fibres; Green-Naghdi): it is driven by the increment of the
   * logarithmic (Hencky) strain
   * @f[
   *   \Delta\hat{\boldsymbol{\varepsilon}} = \ln\boldsymbol{U}^{(n+1)} - \ln\boldsymbol{U}^{(n)},
   *   \qquad \ln\boldsymbol{U} = \tfrac12\ln\boldsymbol{C},
   * @f]
   * and its stress is read as the rotated Kirchhoff stress, which is rotated out,
   * @f[
   *   \boldsymbol{\tau} = \boldsymbol{R}\,\hat{\boldsymbol{\tau}}\,\boldsymbol{R}^T .
   * @f]
   *
   * **Why this frame.** A superposed rigid rotation changes @f$ \boldsymbol{R} @f$ and leaves @f$ \boldsymbol{U} @f$ --
   * and hence everything the wrapped material sees -- unchanged, so the Kirchhoff stress rotates exactly, whatever the
   * wrapped material does with its own state. That is what Marmot::Materials::GradientEnhancedHughesWingetWrapper
   * cannot offer a material that recomputes its stress from total strains (GM7 does: its stress is the
   * homogenized microplane stress of its stored strains, and the stress it is handed is overwritten); it rotates the
   * stress it hands over, which such a material discards. Tensor-valued internal variables of the wrapped material
   * (GM7's microplane state, its stored strains) live in the rotating frame too, so they are objective as well.
   *
   * **Why the logarithmic strain.** It sums to a total strain, @f$ \sum\Delta\hat{\boldsymbol{\varepsilon}} =
   * \ln\boldsymbol{U} @f$, independent of the path, so a material that accumulates its strain stores the Hencky
   * strain; and under a pure stretch @f$ \lambda @f$ it hands the material @f$ \ln\lambda = (\lambda-1) + O((\lambda-
   * 1)^2) @f$, i.e. the wrapped model is recovered up to second order in the strain. The rotated Kirchhoff stress is
   * the stress work-conjugate to @f$ \ln\boldsymbol{U} @f$ for an isotropic response, so the elastic and dissipated
   * energies the wrapped material reports are densities per unit reference volume and are passed through unchanged.
   *
   * The nonlocal fields are handed over as they are (see GradientEnhancedHypoElasticWrapper): they are referred to the
   * rotating frame as well -- for GM7, the nonlocal strain is the nonlocal counterpart of
   * @f$ \hat{\boldsymbol{\varepsilon}} @f$ -- and their balance is written in the reference configuration, so they are
   * invariant under a superposed rigid rotation.
   *
   * | this interface reports | from the wrapped material |
   * |---|---|
   * | @c response.tau | @f$ \boldsymbol{R}\,\hat{\boldsymbol{\tau}}\,\boldsymbol{R}^T @f$ |
   * | @c response.L | @c KLocal |
   * | @c response.nonLocalRadius | @f$ \sqrt{c_K} @f$ |
   * | @c tangents.dTau_dN[K] | @f$
   * \boldsymbol{R}\,(\partial\hat{\boldsymbol{\tau}}/\partial\bar{N}_K)\,\boldsymbol{R}^T @f$ | | @c tangents.dL_dF[K]
   * | @f$ \partial K^{\text{local}}_K/\partial\hat{\boldsymbol{\varepsilon}} :
   * \partial\ln\boldsymbol{U}/\partial\boldsymbol{F} @f$ | | @c tangents.dL_dN | @c dKLocalddK |
   *
   * **Tangent.** Consistent and analytic:
   * @f[
   *   \frac{\partial\tau_{ij}}{\partial F_{kL}} =
   *     \frac{\partial R_{ia}}{\partial F_{kL}}\,\hat\tau_{ab}\,R_{jb} + R_{ia}\,\hat\tau_{ab}\,
   *     \frac{\partial R_{jb}}{\partial F_{kL}} + R_{ia}R_{jb}\,\mathbb{C}_{abcd}\,
   *     \frac{\partial(\ln U)_{cd}}{\partial F_{kL}},
   *   \qquad \frac{\partial\boldsymbol{R}}{\partial F_{kL}} = \Bigl(\frac{\partial\boldsymbol{F}}{\partial F_{kL}}
   *     - \boldsymbol{R}\,\frac{\partial\boldsymbol{U}}{\partial F_{kL}}\Bigr)\boldsymbol{U}^{-1},
   * @f]
   * with @f$ \mathbb{C} @f$ the algorithmic tangent of the wrapped material and the derivatives of the isotropic
   * tensor functions @f$ \ln\boldsymbol{U} = \tfrac12\ln\boldsymbol{C} @f$ and @f$ \boldsymbol{U} =
   * \boldsymbol{C}^{1/2} @f$ from the spectral decomposition of @f$ \boldsymbol{C} @f$ (spectralDerivative()). It is
   * exact for any wrapped material whose stress does not depend on the stress it is handed (GM7), and for an
   * incremental one too, since the stress it is handed is the rotated stress of the last call, which does not depend
   * on @f$ \boldsymbol{F} @f$. The explicit path, computeStressExplicit(), skips it.
   *
   * **What it refuses.** @f$ \boldsymbol{C} = \boldsymbol{F}^T\boldsymbol{F} @f$ is the same for
   * @f$ \boldsymbol{F} @f$ and for its mirror image, so the strain does not show that a material point is turned
   * inside out (@f$ \det\boldsymbol{F} < 0 @f$; @f$ \boldsymbol{F} = -\boldsymbol{1} @f$ would be a point without
   * strain or stress, with a reflection for @f$ \boldsymbol{R} @f$). A deformation gradient without a positive
   * determinant is refused (std::invalid_argument).
   *
   * The state of the wrapper: @f$ \ln\boldsymbol{U} @f$ and @f$ \hat{\boldsymbol{\tau}} @f$ of the last call (both in
   * Voigt notation, the strain with engineering shear), the nonlocal fields of the last call, and the state of the
   * wrapped material. A zero state is the undeformed, stress-free one, so no initialization beyond the wrapped
   * material's is needed.
   */
  template < typename BaseMaterialType >
  class GradientEnhancedCorotationalWrapper : public GradientEnhancedHypoElasticWrapper< BaseMaterialType > {

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
    /// Name of the state layout slot holding @f$ \ln\boldsymbol{U} @f$ of the last call (Voigt, engineering shear).
    static constexpr const char* logarithmicStrainSlot = "Corotational_logU_n";
    /// Name of the state layout slot holding the rotated Kirchhoff stress of the last call (Voigt).
    static constexpr const char* rotatedStressSlot = "Corotational_tau_n";
    /// Name of the state layout slot holding the nonlocal fields of the last call.
    static constexpr const char* nonlocalFieldSlot = "Corotational_N_n";

  public:
    /**
     * @brief Construct the wrapper and the wrapped material.
     * @param[in] matProperties_       Material property values; forwarded **unchanged** to the wrapped material.
     * @param[in] nMaterialProperties_ Number of material property values.
     * @param[in] materialNumber_      Unique identifier for this material instance.
     */
    GradientEnhancedCorotationalWrapper( const double* matProperties_, int nMaterialProperties_, int materialNumber_ )
      : Wrapper( matProperties_, nMaterialProperties_, materialNumber_ )
    {
      this->stateLayout.add( logarithmicStrainSlot, 6 );
      this->stateLayout.add( rotatedStressSlot, 6 );
      this->stateLayout.add( nonlocalFieldSlot, nNonlocal );
      this->stateLayout.add( baseMaterialSlot, baseMaterial->getNumberOfRequiredStateVars() );
      this->stateLayout.finalize();
    }

    // The four-argument override below would otherwise hide the five-argument eigen-deformation overload.
    using Wrapper::Interface::computeStress;

    /**
     * @brief Initialise the state: zero logarithmic strain, stress and nonlocal fields of the last call, plus the
     *        wrapped material's own state.
     * @param[in,out] stateVars  Pointer to the state variable array.
     * @param[in]     nStateVars Number of state variables.
     */
    void initializeYourself( double* stateVars, int nStateVars ) override
    {
      std::fill( stateVars, this->baseMaterialState( stateVars ), 0.0 );
      this->initializeBaseMaterial( stateVars );
    }

    /**
     * @brief Corotational stress update of the wrapped material, with the consistent tangent.
     * @param[in,out] response      Constitutive response; carries the state variables and the energies in, the
     *                              Kirchhoff stress, the local driving forces and the nonlocal radii out.
     * @param[out]    tangents      Algorithmic moduli.
     * @param[in]     deformation   Current deformation gradient and nonlocal fields.
     * @param[in]     timeIncrement Current time and time increment.
     */
    void computeStress( ConstitutiveResponse< 3 >& response,
                        AlgorithmicModuli< 3 >&    tangents,
                        const Deformation< 3 >&    deformation,
                        const TimeIncrement&       timeIncrement ) const override
    {
      computeStressCore( response, &tangents, deformation, timeIncrement );
    }

    /**
     * @brief Corotational stress update without the tangent, for explicit time integration.
     * @param[in,out] response      Constitutive response.
     * @param[in]     deformation   Current deformation gradient and nonlocal fields.
     * @param[in]     timeIncrement Current time and time increment.
     */
    void computeStressExplicit( ConstitutiveResponse< 3 >& response,
                                const Deformation< 3 >&    deformation,
                                const TimeIncrement&       timeIncrement ) const override
    {
      computeStressCore( response, nullptr, deformation, timeIncrement );
    }

    /**
     * @brief Derivative of an isotropic function of a symmetric tensor (Daleckii-Krein formula).
     * @param[in] eigenvalues     Eigenvalues @f$ \lambda_a @f$ of the tensor.
     * @param[in] eigenvectors    Its eigenvectors, by column.
     * @param[in] dividedDifference @f$ \theta(\lambda_a, \lambda_b) @f$: @f$ (f(\lambda_a) - f(\lambda_b)) /
     *                            (\lambda_a - \lambda_b) @f$, and @f$ f'(\lambda_a) @f$ for equal arguments.
     * @return @f$ \partial f(\boldsymbol{C})_{ij}/\partial C_{IJ} = \sum_{a,b} \theta(\lambda_a, \lambda_b)\,
     *         Q_{ia} Q_{jb} Q_{Ia} Q_{Jb} @f$, for a symmetric increment of @f$ \boldsymbol{C} @f$.
     *
     * @details The divided difference is continuous where eigenvalues coincide, so the formula holds for repeated
     * eigenvalues as well, whatever eigenvectors the decomposition picks within the eigenspace.
     */
    static Marmot::Math::Tensors::FastorTensors::Tensor3333d spectralDerivative(
      const Eigen::Vector3d&                           eigenvalues,
      const Eigen::Matrix3d&                           eigenvectors,
      const std::function< double( double, double ) >& dividedDifference )
    {
      Eigen::Matrix3d theta;
      for ( int a = 0; a < 3; a++ )
        for ( int b = 0; b < 3; b++ )
          theta( a, b ) = dividedDifference( eigenvalues( a ), eigenvalues( b ) );

      const Eigen::Matrix3d&                            Q = eigenvectors;
      Marmot::Math::Tensors::FastorTensors::Tensor3333d df_dC( 0.0 );
      for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ )
          for ( int I = 0; I < 3; I++ )
            for ( int J = 0; J < 3; J++ )
              for ( int a = 0; a < 3; a++ )
                for ( int b = 0; b < 3; b++ )
                  df_dC( i, j, I, J ) += theta( a, b ) * Q( i, a ) * Q( j, b ) * Q( I, a ) * Q( J, b );
      return df_dC;
    }

  protected:
    /**
     * @brief The corotational update proper.
     * @param[in,out] response      Constitutive response.
     * @param[out]    tangents      Algorithmic moduli, or nullptr for none.
     * @param[in]     deformation   Current deformation gradient and nonlocal fields.
     * @param[in]     timeIncrement Current time and time increment.
     * @throws std::invalid_argument if the determinant of the deformation gradient is not positive: the polar
     *                               decomposition would then give a reflection for the rotation, and
     *                               @f$ \boldsymbol{C} @f$ a strain as if nothing were wrong.
     */
    void computeStressCore( ConstitutiveResponse< 3 >& response,
                            AlgorithmicModuli< 3 >*    tangents,
                            const Deformation< 3 >&    deformation,
                            const TimeIncrement&       timeIncrement ) const
    {
      using namespace Marmot::Math::Tensors::FastorTensors;
      using namespace Marmot::ContinuumMechanics::Voigt;

      Eigen::Map< Marmot::Math::Vector6d >
        logU_n = this->stateLayout.template getAs< Eigen::Map< Marmot::Math::Vector6d > >( response.stateVars,
                                                                                           logarithmicStrainSlot );
      Eigen::Map< Marmot::Math::Vector6d >
        tauHat_n = this->stateLayout.template getAs< Eigen::Map< Marmot::Math::Vector6d > >( response.stateVars,
                                                                                             rotatedStressSlot );
      Eigen::Map< NonlocalVector >
        N_n = this->stateLayout.template getAs< Eigen::Map< NonlocalVector > >( response.stateVars, nonlocalFieldSlot );

      // --- polar decomposition F = R U from the spectral decomposition of C = F^T F ---
      const Eigen::Matrix3d F = toEigen( deformation.F );

      // C = F^T F does not tell an inverted material point from a sound one (F = -1 gives C = 1, no strain and no
      // stress, with a reflection for R), and a comparison with a number that is not one is false: both are refused
      // by asking for a positive determinant.
      if ( !( F.determinant() > 0.0 ) )
        throw std::invalid_argument( Marmot::Utilities::MakeString()
                                     << __PRETTY_FUNCTION__ << ": the determinant of the deformation gradient is "
                                     << F.determinant() << ", not positive: the material point is inverted or "
                                     << "flattened, or its deformation is not a number" );

      const Eigen::SelfAdjointEigenSolver< Eigen::Matrix3d > eigenC( F.transpose() * F );
      const Eigen::Vector3d&                                 lambda = eigenC.eigenvalues();
      const Eigen::Matrix3d&                                 Q      = eigenC.eigenvectors();

      // a positive determinant that round-off turns into a stretch of zero
      if ( !( lambda.minCoeff() > 0.0 ) )
        throw std::invalid_argument( Marmot::Utilities::MakeString()
                                     << __PRETTY_FUNCTION__ << ": the deformation gradient is singular" );

      const Eigen::Vector3d stretch = lambda.cwiseSqrt();
      const Eigen::Matrix3d UInv    = Q * stretch.cwiseInverse().asDiagonal() * Q.transpose();
      const Eigen::Matrix3d logU    = Q * ( 0.5 * lambda.array().log() ).matrix().asDiagonal() * Q.transpose();
      const Eigen::Matrix3d R       = F * UInv;

      // --- the wrapped material, in the rotating frame ---
      const Marmot::Math::Vector6d logUVoigt = voigtFromStrainMatrix< 3 >( logU );

      const typename BaseMaterial::increment inc = this->baseMaterialIncrement( deformation.N,
                                                                                N_n,
                                                                                logUVoigt - logU_n,
                                                                                timeIncrement );

      typename BaseMaterial::response res;
      res.stress = tauHat_n;
      // the small-strain response has no default member initializers
      res.KLocal.setZero();
      res.c.setZero();
      res.elasticEnergyDensity = response.elasticEnergyDensity;
      res.dissipation          = response.dissipation;
      res.stateVars            = this->baseMaterialState( response.stateVars );

      typename BaseMaterial::tangents tan;
      if ( tangents )
        baseMaterial->computeStress( res, tan, inc );
      else
        baseMaterial->computeStressExplicit( res, inc );

      // --- rotate out ---
      const Eigen::Matrix3d tauHat = stressMatrixFromVoigt< 3 >( res.stress );

      response.tau = fromEigen( R * tauHat * R.transpose() );
      response.L   = res.KLocal;
      for ( int K = 0; K < nNonlocal; K++ )
        response.nonLocalRadius( K ) = std::sqrt( std::max( 0.0, res.c( K ) ) );
      response.elasticEnergyDensity = res.elasticEnergyDensity;
      response.dissipation          = res.dissipation;

      if ( tangents )
        computeTangent( *tangents, F, lambda, Q, UInv, R, tauHat, tan );

      // --- commit ---
      logU_n   = logUVoigt;
      tauHat_n = res.stress;
      N_n      = deformation.N;
    }

    /**
     * @brief The consistent tangent of computeStressCore(), see the class description.
     * @param[out] tangents Algorithmic moduli.
     * @param[in]  F        Deformation gradient.
     * @param[in]  lambda   Eigenvalues of @f$ \boldsymbol{C} @f$.
     * @param[in]  Q        Eigenvectors of @f$ \boldsymbol{C} @f$, by column.
     * @param[in]  UInv     Inverse of the right stretch.
     * @param[in]  R        Polar rotation.
     * @param[in]  tauHat   Rotated Kirchhoff stress.
     * @param[in]  tan      Algorithmic tangents of the wrapped material.
     */
    static void computeTangent( AlgorithmicModuli< 3 >&                tangents,
                                const Eigen::Matrix3d&                 F,
                                const Eigen::Vector3d&                 lambda,
                                const Eigen::Matrix3d&                 Q,
                                const Eigen::Matrix3d&                 UInv,
                                const Eigen::Matrix3d&                 R,
                                const Eigen::Matrix3d&                 tauHat,
                                const typename BaseMaterial::tangents& tan )
    {
      using namespace Marmot::Math::Tensors::FastorTensors;
      using namespace Marmot::ContinuumMechanics::Voigt;

      // d(ln U)/dC = 1/2 d(ln C)/dC; theta = (ln a - ln b)/(a - b), evaluated without cancellation by log1p
      const Tensor3333d dLogU_dC = spectralDerivative( lambda, Q, []( double a, double b ) {
        const double x = ( a - b ) / b;
        return 0.5 * ( x == 0.0 ? 1.0 / b : std::log1p( x ) / ( a - b ) );
      } );
      // dU/dC, theta = (sqrt a - sqrt b)/(a - b) = 1/(sqrt a + sqrt b), exact also for a = b
      const Tensor3333d dU_dC = spectralDerivative( lambda, Q, []( double a, double b ) {
        return 1.0 / ( std::sqrt( a ) + std::sqrt( b ) );
      } );

      // dC_MN/dF_kL = delta_ML F_kN + F_kM delta_NL (M, N rather than I, J: I is the Voigt identity in this scope)
      Tensor3333d dC_dF( 0.0 );
      for ( int M = 0; M < 3; M++ )
        for ( int N = 0; N < 3; N++ )
          for ( int k = 0; k < 3; k++ ) {
            dC_dF( M, N, k, M ) += F( k, N );
            dC_dF( M, N, k, N ) += F( k, M );
          }

      Tensor3333d dLogU_dF( 0.0 );
      Tensor3333d dU_dF( 0.0 );
      for ( int c = 0; c < 3; c++ )
        for ( int d = 0; d < 3; d++ )
          for ( int k = 0; k < 3; k++ )
            for ( int L = 0; L < 3; L++ )
              for ( int M = 0; M < 3; M++ )
                for ( int N = 0; N < 3; N++ ) {
                  dLogU_dF( c, d, k, L ) += dLogU_dC( c, d, M, N ) * dC_dF( M, N, k, L );
                  dU_dF( c, d, k, L ) += dU_dC( c, d, M, N ) * dC_dF( M, N, k, L );
                }

      // dR_ia/dF_kL = ( delta_ik delta_pL - R_ip' dU_p'p/dF_kL ) UInv_pa
      Tensor3333d dR_dF( 0.0 );
      for ( int i = 0; i < 3; i++ )
        for ( int a = 0; a < 3; a++ )
          for ( int k = 0; k < 3; k++ )
            for ( int L = 0; L < 3; L++ ) {
              double value = ( i == k ) ? UInv( L, a ) : 0.0;
              for ( int p = 0; p < 3; p++ )
                for ( int q = 0; q < 3; q++ )
                  value -= R( i, p ) * dU_dF( p, q, k, L ) * UInv( q, a );
              dR_dF( i, a, k, L ) = value;
            }

      // dTauHat_ab/dF_kL = C_abcd d(ln U)_cd/dF_kL; voigtToStiffness scatters the engineering-shear tangent so that
      // the sum over both (c,d) and (d,c) of a symmetric strain derivative counts each shear once
      const Tensor3333d C = voigtToStiffness( tan.dStressddStrain );
      Tensor3333d       dTauHat_dF( 0.0 );
      for ( int a = 0; a < 3; a++ )
        for ( int b = 0; b < 3; b++ )
          for ( int k = 0; k < 3; k++ )
            for ( int L = 0; L < 3; L++ )
              for ( int c = 0; c < 3; c++ )
                for ( int d = 0; d < 3; d++ )
                  dTauHat_dF( a, b, k, L ) += C( a, b, c, d ) * dLogU_dF( c, d, k, L );

      // tau_ij = R_ia tauHat_ab R_jb
      const Eigen::Matrix3d tauHatRT = tauHat * R.transpose();
      const Eigen::Matrix3d RTauHat  = R * tauHat;
      for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ )
          for ( int k = 0; k < 3; k++ )
            for ( int L = 0; L < 3; L++ ) {
              double value = 0.0;
              for ( int a = 0; a < 3; a++ ) {
                value += dR_dF( i, a, k, L ) * tauHatRT( a, j ) + RTauHat( i, a ) * dR_dF( j, a, k, L );
                for ( int b = 0; b < 3; b++ )
                  value += R( i, a ) * dTauHat_dF( a, b, k, L ) * R( j, b );
              }
              tangents.dTau_dF( i, j, k, L ) = value;
            }

      for ( int K = 0; K < nNonlocal; K++ ) {
        const Eigen::Matrix3d dTauHat_dN = stressMatrixFromVoigt< 3 >(
          Marmot::Math::Vector6d( tan.dStressddK.col( K ) ) );
        tangents.dTau_dN[K] = fromEigen( R * dTauHat_dN * R.transpose() );

        // the local driving force's derivative with respect to the engineering-shear Voigt strain, scattered like a
        // stress, so that the sum over all nine (c,d) counts each shear once
        const Eigen::Matrix3d dL_dLogU = stressMatrixFromVoigt< 3 >(
          Marmot::Math::Vector6d( tan.dKLocalddStrain.row( K ).transpose() ) );
        Tensor33d dL_dF( 0.0 );
        for ( int k = 0; k < 3; k++ )
          for ( int L = 0; L < 3; L++ )
            for ( int c = 0; c < 3; c++ )
              for ( int d = 0; d < 3; d++ )
                dL_dF( k, L ) += dL_dLogU( c, d ) * dLogU_dF( c, d, k, L );
        tangents.dL_dF[K] = dL_dF;
      }
      tangents.dL_dN = tan.dKLocalddK;
    }
  };

} // namespace Marmot::Materials
