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
#include "Marmot/MarmotDeformationMeasures.h"
#include "Marmot/MarmotEnergyDensityFunctions.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotFiniteStrainPlasticity.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTypedefs.h"
#include <stdexcept>
#include <string>
#include <tuple>

namespace Marmot::Materials {

  using namespace Fastor;
  using namespace FastorStandardTensors;
  using namespace FastorIndices;

  /**
   * @class Marmot::Materials::BergstromBoyce
   * @brief Classical Bergström-Boyce finite-strain viscoelastic-viscoplastic model.
   *
   * Two networks act in parallel:
   *  - Network A (equilibrium, hyperelastic): sees the total deformation directly.
   *  - Network B ("Maxwell-like", viscous): multiplicative split \f$\boldsymbol
   *    F=\boldsymbol F^{\rm e}\boldsymbol F^{\rm v}\f$, with a chain-stretch /
   *    deviatoric-Mandel-stress power-law flow rule (no yield surface -- flow is
   *    always active).
   *
   * Both networks share the same selectable hyperelastic base potential
   * (#HyperelasticBase), network A evaluated on the total \f$\boldsymbol C\f$ and
   * network B's spring on the elastic \f$\boldsymbol C^{\rm e}\f$, each with its own
   * set of coefficients and bulk modulus.
   *
   * @par Hyperelastic base potentials (all four share the same isochoric-volumetric
   * split \f$\Psi=\Psi_{\rm iso}(\bar I_1,\bar I_2)+\frac{\kappa}{8}(\ln\det
   * \boldsymbol C)^2\f$, with \f$\bar I_1=I_1(\det\boldsymbol C)^{-1/3}\f$,
   * \f$\bar I_2=I_2(\det\boldsymbol C)^{-2/3}\f$ the isochoric invariants of
   * \f$\boldsymbol C\f$ -- since these are invariant under \f$\boldsymbol C\to
   * \lambda\boldsymbol C\f$, all four are exactly stress-free at \f$\boldsymbol
   * C=\boldsymbol I\f$ with no linear-shift correction needed)
   * - @b NeoHooke (coefficient 1 = \f$\mu\f$, coefficients 2,3 unused):
   *   \f$\Psi_{\rm iso}=\frac{\mu}{2}(\bar I_1-3)\f$
   * - @b Yeoh (coefficients 1,2,3 = \f$C_{10},C_{20},C_{30}\f$):
   *   \f$\Psi_{\rm iso}=C_{10}(\bar I_1-3)+C_{20}(\bar I_1-3)^2+C_{30}(\bar I_1-3)^3\f$
   * - @b MooneyRivlin (coefficients 1,2 = \f$C_{10},C_{01}\f$, coefficient 3 unused):
   *   \f$\Psi_{\rm iso}=C_{10}(\bar I_1-3)+C_{01}(\bar I_2-3)\f$
   * - @b ArrudaBoyce (coefficients 1,2 = \f$\mu,\lambda_L\f$, coefficient 3 unused):
   *   isochoric 8-chain potential from
   *   Marmot::ContinuumMechanics::EnergyDensityFunctions::FirstOrderDerived::ArrudaBoyce8ChainPotential
   *   (shared with CompressibleFiniteStrainLinearViscoelasticity); reduces exactly
   *   to #NeoHooke's isochoric part as \f$\lambda_L\to\infty\f$ (both being
   *   functions of the same isochoric \f$\bar I_1\f$ now).
   *
   * All four isochoric parts and the shared volumetric term are themselves shared
   * with CompressibleFiniteStrainLinearViscoelasticity via
   * Marmot::ContinuumMechanics::EnergyDensityFunctions::FirstOrderDerived.
   *
   * @par Material parameters
   * - @b #hyperelasticBase -- 0 = NeoHooke, 1 = Yeoh, 2 = MooneyRivlin, 3 = ArrudaBoyce
   * - @b #kappaA -- network A bulk modulus
   * - @b #kappaB -- network B (spring) bulk modulus
   * - @b #A1, #A2, #A3 -- network A hyperelastic coefficients (meaning depends on #hyperelasticBase)
   * - @b #B1, #B2, #B3 -- network B hyperelastic coefficients (meaning depends on #hyperelasticBase)
   * - @b #c1     -- flow-rate prefactor
   * - @b #c2     -- chain-stretch exponent
   * - @b #c3     -- deviatoric-Mandel-stress-magnitude exponent
   * - @b #implementationType -- algorithm selector (see below)
   * - @b #density (optional) -- density
   *
   * @par State variables
   * - @b Fv -- viscous deformation gradient of network B
   *
   * @par Implementation variants (implementationType)
   * - @b 0: CSDA -- Full return mapping; derivatives via complex-step differentiation
   * - @b 1: Full return mapping; all derivatives computed analytically (not yet implemented)
   */
  class BergstromBoyce : public MarmotMaterialFiniteStrain {
  public:
    using MarmotMaterialFiniteStrain::MarmotMaterialFiniteStrain;

    /// @brief Selectable hyperelastic base potential for both networks.
    enum HyperelasticBase { NeoHooke = 0, Yeoh = 1, MooneyRivlin = 2, ArrudaBoyce = 3 };

    /** Hyperelastic base potential selector (read from @c materialProperties[0]) */
    const int hyperelasticBase;
    /** Network A bulk modulus (read from @c materialProperties[1]) */
    const double kappaA;
    /** Network B bulk modulus (read from @c materialProperties[2]) */
    const double kappaB;
    /** Network A hyperelastic coefficient 1 (read from @c materialProperties[3]) */
    const double A1;
    /** Network A hyperelastic coefficient 2 (read from @c materialProperties[4]) */
    const double A2;
    /** Network A hyperelastic coefficient 3 (read from @c materialProperties[5]) */
    const double A3;
    /** Network B hyperelastic coefficient 1 (read from @c materialProperties[6]) */
    const double B1;
    /** Network B hyperelastic coefficient 2 (read from @c materialProperties[7]) */
    const double B2;
    /** Network B hyperelastic coefficient 3 (read from @c materialProperties[8]) */
    const double B3;
    /** Flow-rate prefactor (read from @c materialProperties[9]) */
    const double c1;
    /** Chain-stretch exponent (read from @c materialProperties[10]) */
    const double c2;
    /** Deviatoric-Mandel-stress-magnitude exponent (read from @c materialProperties[11]) */
    const double c3;

    /** Algorithm variant selector (read from @c materialProperties[12]). */
    const int implementationType;
    /** Density (read from @c materialProperties[13]) (if provided). */
    const double density;

    /**
     * @brief Construct the Bergström-Boyce model.
     * @param materialProperties Array with parameters: #hyperelasticBase, #kappaA, #kappaB, #A1, #A2,
     * #A3, #B1, #B2, #B3, #c1, #c2, #c3, #implementationType, #density (optional).
     * @param nMaterialProperties Length of @c materialProperties.
     * @param materialLabel Material label.
     */
    BergstromBoyce( const double* materialProperties, int nMaterialProperties, int materialLabel );

    void computeStress( ConstitutiveResponse< 3 >& response,
                        AlgorithmicModuli< 3 >&    tangents,
                        const Deformation< 3 >&    deformation,
                        const TimeIncrement&       timeIncrement ) const override;

    /** @brief Full return mapping; all derivatives computed analytically.
     *  @note Not yet implemented -- use CSDA (implementationType = 0).
     */
    void computeStressWithFullReturnMapping( ConstitutiveResponse< 3 >& response,
                                             AlgorithmicModuli< 3 >&    tangents,
                                             const Deformation< 3 >&    deformation,
                                             const TimeIncrement&       timeIncrement ) const;

    /** @brief Full return mapping; derivatives via complex-step differentiation approximation (CSDA). */
    void computeStressCSDA( ConstitutiveResponse< 3 >& response,
                            AlgorithmicModuli< 3 >&    tangents,
                            const Deformation< 3 >&    deformation,
                            const TimeIncrement&       timeIncrement ) const;

    /**
     * @brief Get material density.
     * @return Density value.
     */
    double getDensity( const double* stateVars ) const override
    {
      if ( this->nMaterialProperties < 14 ) {
        throw std::runtime_error(
          std::string( MakeString() << __PRETTY_FUNCTION__ << ": No density given! nMaterialProperties < 14." ) );
      }
      return this->density;
    }

    /** @brief Initialize state (sets @f$\boldsymbol F^{\rm v} = \boldsymbol I@f$) */
    void initializeYourself( double* stateVars, int nStateVars ) override;

    /**
     * @brief Compressible neo-Hookean potential and its first derivative w.r.t. C (templated scalar type).
     * @tparam T Scalar type (double, complex).
     * @param C  Right Cauchy-Green-like tensor the potential is evaluated on.
     * @param mu Shear modulus.
     * @param kappa Bulk modulus.
     * @return \f$\{\,\Psi,\;\partial\Psi/\partial\boldsymbol C\,\}\f$.
     *
     * \f[
     *   \Psi(\boldsymbol C) = \frac{\mu}{2}\left(\bar I_1 - 3\right) + \frac{\kappa}{8}(\ln\det\boldsymbol C)^2,
     *   \qquad \bar I_1 = I_1 (\det\boldsymbol C)^{-1/3}
     * \f]
     * (isochoric part shared with CompressibleFiniteStrainLinearViscoelasticity via
     * Marmot::ContinuumMechanics::EnergyDensityFunctions::FirstOrderDerived::NeoHookePotential,
     * volumetric part shared via
     * Marmot::ContinuumMechanics::EnergyDensityFunctions::FirstOrderDerived::VolumetricPenaltyPotential).
     */
    template < typename T >
    std::tuple< T, Tensor33t< T > > neoHookePotential( const Tensor33t< T >& C,
                                                       const double          mu,
                                                       const double          kappa ) const
    {
      auto [psiIso, dPsiIso_dC] = ContinuumMechanics::EnergyDensityFunctions::FirstOrderDerived::NeoHookePotential< T >(
        C, mu );
      auto [psiVol, dPsiVol_dC] =
        ContinuumMechanics::EnergyDensityFunctions::FirstOrderDerived::VolumetricPenaltyPotential< T >( C, kappa );

      return { psiIso + psiVol, evaluate( dPsiIso_dC + dPsiVol_dC ) };
    }

    /**
     * @brief Yeoh hyperelastic potential and its first derivative w.r.t. C (templated scalar type).
     * @tparam T Scalar type (double, complex).
     * @param C   Right Cauchy-Green-like tensor the potential is evaluated on.
     * @param C10 First Yeoh coefficient.
     * @param C20 Second Yeoh coefficient.
     * @param C30 Third Yeoh coefficient.
     * @param kappa Bulk modulus.
     * @return \f$\{\,\Psi,\;\partial\Psi/\partial\boldsymbol C\,\}\f$. Reduces exactly to
     * #neoHookePotential when @c C20=C30=0 and @c C10=mu/2.
     *
     * \f[
     *   \Psi(\boldsymbol C) = C_{10}(\bar I_1-3) + C_{20}(\bar I_1-3)^2 + C_{30}(\bar I_1-3)^3
     *                       + \frac{\kappa}{8}(\ln\det\boldsymbol C)^2, \qquad \bar I_1 = I_1 (\det\boldsymbol C)^{-1/3}
     * \f]
     * (isochoric part shared with CompressibleFiniteStrainLinearViscoelasticity via
     * Marmot::ContinuumMechanics::EnergyDensityFunctions::FirstOrderDerived::YeohPotential,
     * volumetric part shared via
     * Marmot::ContinuumMechanics::EnergyDensityFunctions::FirstOrderDerived::VolumetricPenaltyPotential).
     */
    template < typename T >
    std::tuple< T, Tensor33t< T > > yeohPotential( const Tensor33t< T >& C,
                                                   const double          C10,
                                                   const double          C20,
                                                   const double          C30,
                                                   const double          kappa ) const
    {
      auto [psiIso, dPsiIso_dC] = ContinuumMechanics::EnergyDensityFunctions::FirstOrderDerived::YeohPotential< T >(
        C, C10, C20, C30 );
      auto [psiVol, dPsiVol_dC] =
        ContinuumMechanics::EnergyDensityFunctions::FirstOrderDerived::VolumetricPenaltyPotential< T >( C, kappa );

      return { psiIso + psiVol, evaluate( dPsiIso_dC + dPsiVol_dC ) };
    }

    /**
     * @brief Mooney-Rivlin hyperelastic potential and its first derivative w.r.t. C (templated scalar
     * type).
     * @tparam T Scalar type (double, complex).
     * @param C   Right Cauchy-Green-like tensor the potential is evaluated on.
     * @param C10 First Mooney-Rivlin coefficient.
     * @param C01 Second Mooney-Rivlin coefficient.
     * @param kappa Bulk modulus.
     * @return \f$\{\,\Psi,\;\partial\Psi/\partial\boldsymbol C\,\}\f$. Reduces exactly to
     * #neoHookePotential when @c C01=0 and @c C10=mu/2.
     *
     * \f[
     *   \Psi(\boldsymbol C) = C_{10}(\bar I_1-3) + C_{01}(\bar I_2-3) + \frac{\kappa}{8}(\ln\det\boldsymbol C)^2,
     *   \quad \bar I_1 = I_1 (\det\boldsymbol C)^{-1/3}, \quad \bar I_2 = I_2 (\det\boldsymbol C)^{-2/3},\quad
     *   I_2 = \frac{1}{2}\left(I_1^2 - {\rm tr}(\boldsymbol C^2)\right)
     * \f]
     * (isochoric part shared with CompressibleFiniteStrainLinearViscoelasticity via
     * Marmot::ContinuumMechanics::EnergyDensityFunctions::FirstOrderDerived::MooneyRivlinPotential,
     * volumetric part shared via
     * Marmot::ContinuumMechanics::EnergyDensityFunctions::FirstOrderDerived::VolumetricPenaltyPotential).
     */
    template < typename T >
    std::tuple< T, Tensor33t< T > > mooneyRivlinPotential( const Tensor33t< T >& C,
                                                           const double          C10,
                                                           const double          C01,
                                                           const double          kappa ) const
    {
      auto [psiIso, dPsiIso_dC] =
        ContinuumMechanics::EnergyDensityFunctions::FirstOrderDerived::MooneyRivlinPotential< T >( C, C10, C01 );
      auto [psiVol, dPsiVol_dC] =
        ContinuumMechanics::EnergyDensityFunctions::FirstOrderDerived::VolumetricPenaltyPotential< T >( C, kappa );

      return { psiIso + psiVol, evaluate( dPsiIso_dC + dPsiVol_dC ) };
    }

    /**
     * @brief Arruda-Boyce 8-chain potential (isochoric part, shared with
     * CompressibleFiniteStrainLinearViscoelasticity via
     * Marmot::ContinuumMechanics::EnergyDensityFunctions::FirstOrderDerived::ArrudaBoyce8ChainPotential)
     * plus this material's own \f$\frac{\kappa}{8}(\ln\det\boldsymbol C)^2\f$
     * volumetric term.
     * @tparam T Scalar type (double, complex).
     * @param C Right Cauchy-Green-like tensor the potential is evaluated on.
     * @param mu Shear-modulus-like parameter.
     * @param lambdaL Locking stretch.
     * @param kappa Bulk modulus.
     * @return \f$\{\,\Psi,\;\partial\Psi/\partial\boldsymbol C\,\}\f$. The isochoric
     * part reduces exactly to \f$\frac{\mu}{2}(\bar I_1-3)\f$ as @c lambdaL ->
     * infinity -- NOT to #neoHookePotential's full compressible response, since that
     * is expressed in the raw invariant \f$I_1\f$ rather than the isochoric
     * \f$\bar I_1\f$ used here (the two agree only at \f$\boldsymbol
     * C=\boldsymbol I\f$; see the class-level @par Hyperelastic base potentials
     * block for why the two gradients differ in general). Unlike #yeohPotential /
     * #mooneyRivlinPotential, no linear-shift term is needed to stay stress-free at
     * \f$\boldsymbol C=\boldsymbol I\f$: the shared potential's isochoric invariant
     * \f$\bar I_1\f$ is already exactly stress-free there on its own.
     *
     * \f[
     *   \Psi(\boldsymbol C) = \Psi_{AB}(\bar I_1;\mu,\lambda_L) + \frac{\kappa}{8}(\ln\det\boldsymbol C)^2
     * \f]
     */
    template < typename T >
    std::tuple< T, Tensor33t< T > > arrudaBoycePotential( const Tensor33t< T >& C,
                                                          const double          mu,
                                                          const double          lambdaL,
                                                          const double          kappa ) const
    {
      auto [psiIso, dPsiIso_dC] =
        ContinuumMechanics::EnergyDensityFunctions::FirstOrderDerived::ArrudaBoyce8ChainPotential< T >( C,
                                                                                                        mu,
                                                                                                        lambdaL );
      auto [psiVol, dPsiVol_dC] =
        ContinuumMechanics::EnergyDensityFunctions::FirstOrderDerived::VolumetricPenaltyPotential< T >( C, kappa );

      return { psiIso + psiVol, evaluate( dPsiIso_dC + dPsiVol_dC ) };
    }

    /**
     * @brief Dispatches to the selected #HyperelasticBase potential (templated scalar type).
     * @tparam T Scalar type (double, complex).
     * @param C     Right Cauchy-Green-like tensor the potential is evaluated on.
     * @param base  Selected #HyperelasticBase.
     * @param p1    Coefficient 1 (mu, or C10).
     * @param p2    Coefficient 2 (unused, or C20/C01, or lambdaL for ArrudaBoyce).
     * @param p3    Coefficient 3 (unused, or C30/unused).
     * @param kappa Bulk modulus.
     * @return \f$\{\,\Psi,\;\partial\Psi/\partial\boldsymbol C\,\}\f$.
     */
    template < typename T >
    std::tuple< T, Tensor33t< T > > hyperelasticPotential( const Tensor33t< T >& C,
                                                           const int             base,
                                                           const double          p1,
                                                           const double          p2,
                                                           const double          p3,
                                                           const double          kappa ) const
    {
      switch ( base ) {
      case NeoHooke: return neoHookePotential( C, p1, kappa );
      case Yeoh: return yeohPotential( C, p1, p2, p3, kappa );
      case MooneyRivlin: return mooneyRivlinPotential( C, p1, p2, kappa );
      case ArrudaBoyce: return arrudaBoycePotential( C, p1, p2, kappa );
      default: throw std::invalid_argument( "BergstromBoyce: unknown hyperelasticBase" );
      }
    }

    /**
     * @brief Deviatoric flow direction, flow rate and chain stretch for network B, at a given elastic
     * deformation gradient (templated scalar type).
     * @tparam T Scalar type (double, complex).
     * @param Fe Elastic deformation gradient of network B.
     * @return \f$\{\,\boldsymbol N,\;\rho,\;\dot\gamma\,\}\f$.
     * @details The deviatoric Mandel stress \f$\boldsymbol{\mathcal S}={\rm dev}(2\boldsymbol C^{\rm
     * e}\partial\Psi_B/\partial\boldsymbol C^{\rm e})\f$ is computed from the general hyperelastic
     * potential (#hyperelasticPotential); for the NeoHooke base this reduces to the closed form
     * \f$\boldsymbol{\mathcal S}=\mu_B\,{\rm dev}(\boldsymbol C^{\rm e})\f$, but for Yeoh/Mooney-Rivlin
     * no such shortcut exists since \f$\partial\Psi_B/\partial\boldsymbol C^{\rm e}\f$ is not purely
     * linear in \f$\boldsymbol C^{\rm e}\f$.
     */
    template < typename T >
    std::tuple< Tensor33t< T >, T, T > computeFlowQuantities( const Tensor33t< T >& Fe ) const
    {
      const Tensor33t< T > Ce          = ContinuumMechanics::DeformationMeasures::rightCauchyGreen( Fe );
      const T              I1          = trace( Ce );
      const T              lambdaChain = sqrt( I1 / T( 3.0 ) );

      T              psiB;
      Tensor33t< T > dPsiB_dCe;
      std::tie( psiB, dPsiB_dCe ) = hyperelasticPotential( Ce, hyperelasticBase, B1, B2, B3, kappaB );

      const Tensor33t< T > Mandel = 2. * einsum< IK, KJ, to_IJ >( Ce, dPsiB_dCe );
      const Tensor33t< T > devS   = deviatoric( Mandel );

      // Regularize the RESULT of sqrt, not its argument: at the undeformed reference
      // configuration (devS=0), a complex-step perturbation makes inner(devS,devS)
      // land on a tiny NEGATIVE real number (it is a complex square, not a modulus),
      // and sqrt(negative real) is well-defined and PURELY IMAGINARY -- exactly where
      // the O(h) complex-step derivative information lives. Regularizing the argument
      // (forcing it positive before sqrt) would instead produce a purely REAL sqrt and
      // silently discard that derivative, which previously caused NaNs in the outer
      // Newton solve for stress-controlled steps starting from a fresh, undeformed
      // state (verified: this was traced back to exactly this line during development).
      T normDev = sqrt( Fastor::inner( devS, devS ) );
      if ( Math::makeReal( normDev ) == 0.0 )
        normDev += 1e-15;

      const Tensor33t< T > N   = multiplyFastorTensorWithScalar( devS, T( 1.0 ) / normDev );
      T                    rho = normDev;

      const T stretchTerm = lambdaChain - T( 1.0 );
      // Regularize by SHIFTING the real part up to a small positive floor
      // rather than hard-clamping to a literal T(0.0) (see the normDev
      // regularization above for why a literal discards complex-step
      // derivative information for the CSDA tangent). A hard zero is
      // additionally unsafe here for c2<0: pow(0.0, c2) is +inf for a
      // negative exponent, which happens at lambdaChain=1 -- the reference
      // configuration of every simulation -- propagating to NaN in the
      // outer Newton solve on the very first increment. Adding exactly the
      // shortfall needed to bring the real part up to stretchTermFloor
      // leaves any nonzero imaginary (complex-step) component untouched,
      // and leaves the common case (stretchTerm already above the floor)
      // completely unchanged. The floor value itself was tuned empirically,
      // not derived: 1e-8 avoids the literal inf/NaN but the outer Newton
      // solve still stalls (multi-minute non-convergence, not an exception)
      // for c2 around -0.5 at fast loading rates, since pow(1e-8, c2) is
      // still enormous for c2<0; 1e-2 was the smallest tested value (of
      // 1e-8, 1e-4, 1e-3, 1e-2) that resolved this for every rate group and
      // material tried. Some individual (material, c2, rate) combinations
      // can still fail fast with a clean NaN exception even at this floor
      // (e.g. c2 very close to 0, oddly, more often than c2 near -1) --
      // those are handled by the existing per-residual exception fallback
      // during calibration, not by this floor.
      constexpr double stretchTermFloor = 1e-2;
      const double     stretchTermReal  = Math::makeReal( stretchTerm );
      const T          stretchTermClamped = stretchTermReal > stretchTermFloor
                                               ? stretchTerm
                                               : stretchTerm + T( stretchTermFloor - stretchTermReal );

      const T gammaDot = T( c1 ) * pow( stretchTermClamped, c2 ) * pow( rho, c3 );

      return { N, rho, gammaDot };
    }

  private:
    // helper alias needed inside mooneyRivlinPotential/computeFlowQuantities (declared here since
    // class-scope 'using' of a local OIndex alias must be visible to all member templates)
    using to_IJ = Fastor::OIndex< I_, J_ >;

  public:
    /**
     * @brief Residual vector for the return mapping of network B (templated scalar type).
     * @details Vector X has 10 unknowns: 9 for \f$\boldsymbol F^{\rm e}\f$ (flattened), 1 for \f$\Delta\gamma\f$.
     * @tparam T Scalar type (double, complex).
     * @param X          Current iterate.
     * @param FeTrial    Trial elastic deformation gradient.
     * @param dT         Time increment.
     * @return Residual vector.
     */
    template < typename T >
    VectorXt< T > computeResidualVector( const VectorXt< T >& X, const Tensor33d& FeTrial, const double dT ) const
    {
      using mV9t = Eigen::Map< const Eigen::Matrix< T, 9, 1 > >;
      VectorXt< T > R( 10 );

      const Tensor33t< T > Fe( X.segment( 0, 9 ).data() );
      const T              dGamma = X( 9 );

      Tensor33t< T > N;
      T              rho, gammaDot;
      std::tie( N, rho, gammaDot ) = computeFlowQuantities( Fe );

      const Tensor33t< T > dGp = multiplyFastorTensorWithScalar( N, dGamma );
      const Tensor33t< T > dFv = ContinuumMechanics::FiniteStrain::Plasticity::FlowIntegration::exponentialMap( dGp );

      VectorXt< T > aux = mV9t( Tensor33t< T >( einsum< iJ, JK >( Fe, dFv ) ).data() ) -
                          mV9t( fastorTensorFromDoubleTensor< T >( FeTrial ).data() );

      for ( int i = 0; i < 9; ++i )
        R( i ) = aux( i );

      R( 9 ) = dGamma / T( dT ) - gammaDot;

      return R;
    }
  };

} // namespace Marmot::Materials
