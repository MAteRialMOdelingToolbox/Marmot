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
#include "Marmot/MarmotDeformationMeasures.h"
#include "Marmot/MarmotEnergyDensityFunctions.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotFiniteStrainPlasticity.h"
#include "Marmot/MarmotMaterialGradientEnhancedFiniteStrain.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTypedefs.h"
#include <algorithm>
#include <cmath>
#include <utility>

namespace Marmot::Materials {

  /**
   * @class GradientEnhancedFiniteStrainDruckerPrager
   * @brief Finite-strain Drucker-Prager plasticity with implicit-gradient (nonlocal) damage.
   *
   * The plasticity follows FiniteStrainJ2Plasticity: multiplicative split @f$ \boldsymbol{F} =
   * \boldsymbol{F}^e\boldsymbol{F}^p @f$ with the plastic deformation gradient as state, elasticity of
   * CompressibleNeoHooke (PenceGouPotentialB) in @f$ \boldsymbol{C}^e @f$, a yield function on the Mandel stress
   * @f$ \boldsymbol{M} = \boldsymbol{C}^e\boldsymbol{S} @f$, and the flow integrated with the exponential map,
   * @f$ \boldsymbol{F}^{e,\mathrm{trial}} = \boldsymbol{F}^e \exp( \Delta\lambda\,\partial g/\partial\boldsymbol{M} )
   * @f$. The return mapping solves this equation together with the hardening law and the consistency condition for
   * @f$ \{\boldsymbol{F}^e, \alpha, \Delta\lambda\} @f$ with Newton's method; its Jacobian is computed by the complex
   * step, and the algorithmic tangent follows from the same Jacobian by the implicit function theorem.
   *
   * **Plasticity.** Drucker-Prager cone through the outer edges of the Mohr-Coulomb pyramid, with linear hardening of
   * the cohesion and non-associated flow,
   * @f[
   *   f = \sqrt{J_2} + \eta\,p - \xi\,(c_0 + H\alpha), \qquad g = \sqrt{J_2} + \bar\eta\,p, \qquad
   *   \Delta\alpha = \xi\,\Delta\lambda,
   * @f]
   * @f$ \eta = \frac{6\sin\phi}{\sqrt3\,(3-\sin\phi)} @f$, @f$ \xi = \frac{6\cos\phi}{\sqrt3\,(3-\sin\phi)} @f$,
   * @f$ \bar\eta @f$ as @f$ \eta @f$ with the dilatancy angle. Where the cone return does not exist (a trial state
   * beyond the apex), the state returns to the apex: by isotropy, @f$ \boldsymbol{F}^p @f$ is determined up to a
   * rotation, so @f$ \boldsymbol{F}^e = J_e^{1/3}\boldsymbol{I} @f$ with the unknowns @f$ \{\ln J_e, \alpha\} @f$
   * and @f$ \Delta\alpha = (\xi/\bar\eta)\,\Delta\varepsilon^p_v @f$.
   *
   * **Damage.** The local variable is the accumulated dilatant (volumetric) plastic strain,
   * @f$ \Delta\alpha_\mathrm{local} = \langle\Delta\varepsilon^p_v\rangle @f$ (on the cone
   * @f$ \bar\eta\,\Delta\lambda @f$). It is the source @f$ L @f$ of the nonlocal balance, the damage follows from
   * @f$ \kappa = \max_t( m\bar{N} + (1-m)\alpha_\mathrm{local} ) @f$ with the exponential law
   * @f$ \omega = \min( 1 - e^{-\kappa/\varepsilon_f}, \omega_\max ) @f$, and
   * @f$ \boldsymbol{\tau} = (1-\omega)\,\boldsymbol{\tau}_\mathrm{eff} @f$.
   *
   * Material properties:
   * | idx | symbol                  | meaning                                            |
   * |-----|-------------------------|----------------------------------------------------|
   * | 0   | @f$ K @f$               | bulk modulus                                       |
   * | 1   | @f$ G @f$               | shear modulus                                      |
   * | 2   | @f$ c_0 @f$             | cohesion                                           |
   * | 3   | @f$ \phi @f$            | friction angle [deg]                               |
   * | 4   | @f$ \psi @f$            | dilatancy angle [deg]                              |
   * | 5   | @f$ H @f$               | linear hardening modulus of the cohesion           |
   * | 6   | @f$ \varepsilon_f @f$   | softening modulus of the damage                    |
   * | 7   | @f$ \omega_\max @f$     | maximum damage                                     |
   * | 8   | @f$ l @f$               | nonlocal radius, @f$ c = l^2 @f$                   |
   * | 9   | @f$ m @f$               | weighting of the nonlocal measure in the damage    |
   * | 10  | @f$ \rho @f$            | density in the reference configuration (optional)  |
   *
   * Constraints: @f$ K, G, c_0, \varepsilon_f, l > 0 @f$, @f$ H \geq 0 @f$, @f$ 0 \leq \psi \leq \phi < 90^\circ @f$,
   * @f$ 0 \leq \omega_\max < 1 @f$ and @f$ m \geq 0 @f$ (@f$ m > 1 @f$: over-nonlocal).
   *
   * State variables: @c Fp (9), @c alphaP (hardening variable), @c alphaD (local damage variable),
   * @c kappa (damage history), @c omega (damage).
   *
   * The dissipation is cumulative: the incoming ConstitutiveResponse::dissipation is incremented by
   * @f$ (1-\omega)\,\boldsymbol{M}:\Delta\boldsymbol{\varepsilon}^p + \psi_\mathrm{eff}\,\Delta\omega @f$.
   */
  class GradientEnhancedFiniteStrainDruckerPrager : public MarmotMaterialGradientEnhancedFiniteStrain {

  public:
    template < typename T >
    using Tensor33t    = FastorStandardTensors::Tensor33t< T >;
    using Tensor33d    = FastorStandardTensors::Tensor33d;
    using TensorMap33d = FastorStandardTensors::TensorMap33d;

    /**
     * @brief Construct the material and validate its properties.
     * @param[in] materialProperties Array of the material properties, see the table above.
     * @param[in] nMaterialProperties Length of @p materialProperties (at least 10; 11 with the density).
     * @param[in] materialNumber Material label.
     * @throws std::invalid_argument for a too short property array or a property out of range.
     */
    GradientEnhancedFiniteStrainDruckerPrager( const double* materialProperties,
                                               int           nMaterialProperties,
                                               int           materialNumber );

    using MarmotMaterialGradientEnhancedFiniteStrain::computeStress;

    /**
     * @brief Compute the Kirchhoff stress, the local driving force and the algorithmic tangents.
     * @details Performs the elastic trial; if it is not admissible, the return mapping to the cone or to the apex,
     * then the damage update.
     * @param[in,out] response Kirchhoff stress, local driving force @f$ L @f$, nonlocal radius, energy density,
     * cumulative dissipation and the state variables, updated in place.
     * @param[out] tangents @f$ \partial\boldsymbol{\tau}/\partial\boldsymbol{F} @f$,
     * @f$ \partial\boldsymbol{\tau}/\partial\bar{N} @f$, @f$ \partial L/\partial\boldsymbol{F} @f$ and
     * @f$ \partial L/\partial\bar{N} @f$.
     * @param[in] deformation Deformation gradient and nonlocal field at the end of the increment.
     * @param[in] timeIncrement Time and time increment (not used, the model is rate independent).
     * @throws StressUpdateFailed if the elastic volume ratio is not positive or no admissible return exists.
     */
    void computeStress( ConstitutiveResponse< 3 >& response,
                        AlgorithmicModuli< 3 >&    tangents,
                        const Deformation< 3 >&    deformation,
                        const TimeIncrement&       timeIncrement ) const override;

    /**
     * @brief The density in the reference configuration.
     * @param[in] stateVars State variables (not used).
     * @return The density, property 10.
     * @throws std::runtime_error if the density is not given.
     */
    double getDensity( const double* stateVars ) const override;

    /**
     * @brief Initialize the state variables: zero, except for @f$ \boldsymbol{F}^p = \boldsymbol{I} @f$ (the base
     * default zero-fill would be singular).
     * @param[out] stateVars State variables.
     * @param[in] nStateVars Number of state variables.
     */
    void initializeYourself( double* stateVars, int nStateVars ) override;

    /**
     * @brief Drucker-Prager parameters of the cone through the outer edges of the Mohr-Coulomb pyramid.
     * @param[in] angleInDegrees Friction (or dilatancy) angle in degrees.
     * @return @f$ \{\eta, \xi\} @f$.
     */
    static std::pair< double, double > outerConeParameters( double angleInDegrees );

    /**
     * @brief Second Piola-Kirchhoff stress of the elastic deformation gradient.
     * @tparam T Scalar type (real, or complex for the complex-step Jacobians).
     * @param[in] Fe Elastic deformation gradient.
     * @return @f$ \boldsymbol{S} = 2\,\partial\psi/\partial\boldsymbol{C}^e @f$.
     */
    template < typename T >
    Tensor33t< T > secondPiolaKirchhoff( const Tensor33t< T >& Fe ) const
    {
      using namespace ContinuumMechanics;
      const Tensor33t< T > Ce     = DeformationMeasures::rightCauchyGreen( Fe );
      const auto [psi_, dPsi_dCe] = EnergyDensityFunctions::FirstOrderDerived::PenceGouPotentialB( Ce, K, G );
      return multiplyFastorTensorWithScalar( dPsi_dCe, T( 2.0 ) );
    }

    /**
     * @brief Mandel stress of the elastic deformation gradient.
     * @tparam T Scalar type.
     * @param[in] Fe Elastic deformation gradient.
     * @return @f$ \boldsymbol{M} = \boldsymbol{C}^e\boldsymbol{S} @f$.
     */
    template < typename T >
    Tensor33t< T > mandelStress( const Tensor33t< T >& Fe ) const
    {
      return Tensor33t< T >( ContinuumMechanics::DeformationMeasures::rightCauchyGreen( Fe ) %
                             secondPiolaKirchhoff( Fe ) );
    }

    /**
     * @brief Yield function on the Mandel stress.
     * @tparam T Scalar type.
     * @param[in] M Mandel stress.
     * @param[in] alphaP Hardening variable.
     * @return @f$ f = \sqrt{J_2} + \eta\,p - \xi\,(c_0 + H\alpha) @f$.
     */
    template < typename T >
    T yieldFunction( const Tensor33t< T >& M, const T alphaP ) const
    {
      const Tensor33t< T > dev = deviatoric( M );
      const T              p   = trace( M ) / 3.;
      return sqrt( 0.5 * Fastor::inner( dev, dev ) ) + eta * p - xi * ( c0 + H * alphaP );
    }

    /**
     * @brief Flow direction of the plastic potential.
     * @tparam T Scalar type.
     * @param[in] M Mandel stress, not at the apex, where the cone is singular.
     * @return @f$ \partial g/\partial\boldsymbol{M} @f$.
     */
    template < typename T >
    Tensor33t< T > flowDirection( const Tensor33t< T >& M ) const
    {
      const Tensor33t< T > dev  = deviatoric( M );
      const T              sqJ2 = sqrt( 0.5 * Fastor::inner( dev, dev ) );
      Tensor33t< T >       n    = multiplyFastorTensorWithScalar( dev, T( 0.5 ) / sqJ2 );
      for ( int i = 0; i < 3; i++ )
        n( i, i ) += etaBar / 3.;
      return n;
    }

    /**
     * @brief Residual of the return to the cone.
     * @details The flow rule @f$ \boldsymbol{F}^e\exp(\Delta\lambda\,\partial g/\partial\boldsymbol{M}) -
     * \boldsymbol{F}^{e,\mathrm{trial}} @f$ (with the exponential map by scaling and squaring, which also represents
     * the large plastic increments of intermediate Newton iterates), the hardening law and the yield function scaled
     * by @f$ c_0 @f$.
     * @tparam T Scalar type.
     * @param[in] X Unknowns @f$ \{\boldsymbol{F}^e \text{ (9, row major)}, \alpha, \Delta\lambda\} @f$.
     * @param[in] FeTrial Trial elastic deformation gradient.
     * @param[in] alphaPOld Hardening variable at the beginning of the increment.
     * @return The residual (11).
     */
    template < typename T >
    VectorXt< T > coneResidual( const VectorXt< T >& X, const Tensor33d& FeTrial, const double alphaPOld ) const
    {
      using namespace FastorIndices;
      const Tensor33t< T > Fe( X.segment( 0, 9 ).eval().data() );
      const T              alphaP  = X( 9 );
      const T              dLambda = X( 10 );

      const Tensor33t< T > M   = mandelStress( Fe );
      const Tensor33t< T > dGp = multiplyFastorTensorWithScalar( flowDirection( M ), dLambda );
      const Tensor33t< T >
        dFp = ContinuumMechanics::FiniteStrain::Plasticity::FlowIntegration::exponentialMapScalingAndSquaring( dGp );
      const Tensor33t< T > Fe_dFp = Fastor::einsum< iJ, JK >( Fe, dFp );

      VectorXt< T > R( 11 );
      for ( int i = 0; i < 9; i++ )
        R( i ) = Fe_dFp.data()[i] - T( FeTrial.data()[i] );
      R( 9 )  = alphaP - alphaPOld - xi * dLambda;
      R( 10 ) = yieldFunction( M, alphaP ) / c0;
      return R;
    }

    /**
     * @brief Residual of the return to the apex.
     * @details With @f$ \boldsymbol{F}^e = \exp(\theta/3)\,\boldsymbol{I} @f$: the yield function of the spherical
     * state, scaled by @f$ c_0 @f$, and the hardening law of the volumetric flow.
     * @tparam T Scalar type.
     * @param[in] X Unknowns @f$ \{\theta = \ln J_e, \alpha\} @f$.
     * @param[in] thetaTrial @f$ \ln J_e @f$ of the trial state.
     * @param[in] alphaPOld Hardening variable at the beginning of the increment.
     * @return The residual (2).
     */
    template < typename T >
    VectorXt< T > apexResidual( const VectorXt< T >& X, const double thetaTrial, const double alphaPOld ) const
    {
      const T        theta  = X( 0 );
      const T        alphaP = X( 1 );
      Tensor33t< T > Fe( T( 0.0 ) );
      for ( int i = 0; i < 3; i++ )
        Fe( i, i ) = exp( theta / 3. );
      const T p = trace( mandelStress( Fe ) ) / 3.;

      VectorXt< T > R( 2 );
      R( 0 ) = ( eta * p - xi * ( c0 + H * alphaP ) ) / c0;
      R( 1 ) = alphaP - alphaPOld - xi / etaBar * ( thetaTrial - theta );
      return R;
    }

  protected:
    const double& K;                  ///< bulk modulus
    const double& G;                  ///< shear modulus
    const double& c0;                 ///< cohesion
    const double& frictionAngle;      ///< friction angle in degrees
    const double& dilatancyAngle;     ///< dilatancy angle in degrees
    const double& H;                  ///< linear hardening modulus of the cohesion
    const double& softeningModulus;   ///< softening modulus @f$ \varepsilon_f @f$ of the damage
    const double& maxDamage;          ///< maximum damage
    const double& nonLocalRadius;     ///< nonlocal radius @f$ l @f$
    const double& weightingParameter; ///< weighting @f$ m @f$ of the nonlocal measure in the damage
    const double  eta;                ///< friction parameter @f$ \eta @f$ of the yield function
    const double  xi;                 ///< cohesion parameter @f$ \xi @f$ of the yield function
    const double  etaBar;             ///< dilatancy parameter @f$ \bar\eta @f$ of the plastic potential

    /**
     * @struct ReturnMapping
     * @brief The converged return mapping: the new elastic state, and the sensitivities needed for the tangents.
     */
    struct ReturnMapping {
      bool                               plastic = false; ///< the state has returned to the cone or the apex
      Tensor33d                          Fe;              ///< elastic deformation gradient
      Tensor33d                          FpNew;           ///< plastic deformation gradient
      double                             alphaP;          ///< hardening variable
      double                             dEpVol;          ///< volumetric plastic log strain increment
      double                             plasticWork;     ///< M : dEp
      FastorStandardTensors::Tensor3333d dFe_dF;          ///< d Fe / d F
      Tensor33d                          dDeltaAlphaD_dF; ///< d dalpha_local / d F
    };

    /**
     * @brief The elastic trial and, if it is not admissible, the return to the apex or to the cone.
     * @param[in] F Deformation gradient at the end of the increment.
     * @param[in] FpOld Plastic deformation gradient at the beginning of the increment, a view of the state.
     * @param[in] alphaPOld Hardening variable at the beginning of the increment.
     * @return The converged return mapping.
     * @throws StressUpdateFailed if the elastic volume ratio is not positive or no admissible return exists.
     */
    ReturnMapping returnMapping( const Tensor33d& F, const TensorMap33d& FpOld, double alphaPOld ) const;

    /**
     * @brief The return to the cone.
     * @param[in] FeTrial Trial elastic deformation gradient @f$ \boldsymbol{F}\boldsymbol{F}^{p,-1}_n @f$.
     * @param[in] FpOld Plastic deformation gradient at the beginning of the increment.
     * @param[in] FpOldInv Its inverse.
     * @param[in] alphaPOld Hardening variable at the beginning of the increment.
     * @param[out] converged Whether a solution on the cone was found.
     * @return The converged return mapping (undefined if not @p converged).
     */
    ReturnMapping returnToCone( const Tensor33d&    FeTrial,
                                const TensorMap33d& FpOld,
                                const Tensor33d&    FpOldInv,
                                double              alphaPOld,
                                bool&               converged ) const;

    /**
     * @brief The return to the apex, if it is admissible (otherwise, the state belongs to the cone).
     * @param[in] F Deformation gradient at the end of the increment.
     * @param[in] FeTrial Trial elastic deformation gradient @f$ \boldsymbol{F}\boldsymbol{F}^{p,-1}_n @f$.
     * @param[in] FpOldInv Inverse of the plastic deformation gradient at the beginning of the increment.
     * @param[in] alphaPOld Hardening variable at the beginning of the increment.
     * @param[out] admissible Whether the apex is the solution.
     * @return The converged return mapping (undefined if not @p admissible).
     */
    ReturnMapping returnToApex( const Tensor33d& F,
                                const Tensor33d& FeTrial,
                                const Tensor33d& FpOldInv,
                                double           alphaPOld,
                                bool&            admissible ) const;
  };

} // namespace Marmot::Materials
