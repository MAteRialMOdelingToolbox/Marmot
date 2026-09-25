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
 * Thomas Mader thomas.mader@boku.ac.at
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
#include "Marmot/MarmotTypedefs.h"
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
   * **Damage.** The local variable grows with the volumetric plastic strain over the ductility measure of CDPM2
   * (Grassl et al. 2013), @f$ \Delta\alpha_\mathrm{local} = \Delta\varepsilon^p_v / x_s(R_s) @f$. It is the source @f$
   * L @f$ of the nonlocal balance, the damage follows from @f$ \kappa = \max_t( m\bar{N} + (1-m)\alpha_\mathrm{local} )
   * @f$, and
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
   * | 6   | @f$ A_s @f$             | ductility parameter of the damage                  |
   * | 7   | @f$ \varepsilon_f @f$   | softening modulus of the damage                    |
   * | 8   | @f$ \omega_\max @f$     | maximum damage                                     |
   * | 9   | @f$ l @f$               | nonlocal radius, @f$ c = l^2 @f$                   |
   * | 10  | @f$ m @f$               | weighting of the nonlocal measure in the damage    |
   * | 11  | @f$ \rho @f$            | density in the reference configuration (optional)  |
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
    using Tensor33t = FastorStandardTensors::Tensor33t< T >;
    using Tensor33d = FastorStandardTensors::Tensor33d;

    GradientEnhancedFiniteStrainDruckerPrager( const double* materialProperties,
                                               int           nMaterialProperties,
                                               int           materialNumber );

    using MarmotMaterialGradientEnhancedFiniteStrain::computeStress;

    void computeStress( ConstitutiveResponse< 3 >& response,
                        AlgorithmicModuli< 3 >&    tangents,
                        const Deformation< 3 >&    deformation,
                        const TimeIncrement&       timeIncrement ) const override;

    double getDensity( const double* stateVars ) const override;

    /// Fp must be initialized to the identity (the base default zero-fill would be singular).
    void initializeYourself( double* stateVars, int nStateVars ) override;

    /// Drucker-Prager parameters {eta, xi} for a Mohr-Coulomb angle in degrees (outer cone).
    static std::pair< double, double > outerConeParameters( double angleInDegrees );

    /// second Piola-Kirchhoff stress of the elastic deformation gradient
    template < typename T >
    Tensor33t< T > secondPiolaKirchhoff( const Tensor33t< T >& Fe ) const
    {
      using namespace ContinuumMechanics;
      const Tensor33t< T > Ce     = DeformationMeasures::rightCauchyGreen( Fe );
      const auto [psi_, dPsi_dCe] = EnergyDensityFunctions::FirstOrderDerived::PenceGouPotentialB( Ce, K, G );
      return multiplyFastorTensorWithScalar( dPsi_dCe, T( 2.0 ) );
    }

    /// Mandel stress M = Ce S of the elastic deformation gradient
    template < typename T >
    Tensor33t< T > mandelStress( const Tensor33t< T >& Fe ) const
    {
      return Tensor33t< T >( ContinuumMechanics::DeformationMeasures::rightCauchyGreen( Fe ) %
                             secondPiolaKirchhoff( Fe ) );
    }

    /// effective Kirchhoff stress Fe S Fe^T
    template < typename T >
    Tensor33t< T > kirchhoffStress( const Tensor33t< T >& Fe ) const
    {
      return Tensor33t< T >( Fe % secondPiolaKirchhoff( Fe ) % Fastor::transpose( Fe ) );
    }

    /// yield function on the Mandel stress
    template < typename T >
    T yieldFunction( const Tensor33t< T >& M, const T alphaP ) const
    {
      const Tensor33t< T > dev = deviatoric( M );
      const T              p   = trace( M ) / 3.;
      return sqrt( 0.5 * Fastor::inner( dev, dev ) ) + eta * p - xi * ( c0 + H * alphaP );
    }

    /// flow direction dg/dM of the plastic potential; the cone is singular at its apex
    template < typename T >
    Tensor33t< T > flowDirection( const Tensor33t< T >& M ) const
    {
      const Tensor33t< T > dev  = deviatoric( M );
      const T              sqJ2 = sqrt( 0.5 * Fastor::inner( dev, dev ) );
      Tensor33t< T >       Ivol = multiplyFastorTensorWithScalar( identity< T >(), T( etaBar / 3. ) );
      return Tensor33t< T >( multiplyFastorTensorWithScalar( dev, T( 0.5 ) / sqJ2 ) + Ivol );
    }

    /**
     * Residual of the return to the cone for X = {Fe (9, row major), alphaP, dLambda}:
     * Fe exp( dLambda dg/dM ) - FeTrial, the hardening law and the scaled yield function.
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
      const Tensor33t< T > dFp = ContinuumMechanics::FiniteStrain::Plasticity::FlowIntegration::exponentialMap( dGp );
      const Tensor33t< T > Fe_dFp = Fastor::einsum< iJ, JK >( Fe, dFp );

      VectorXt< T > R( 11 );
      for ( int i = 0; i < 9; i++ )
        R( i ) = Fe_dFp.data()[i] - T( FeTrial.data()[i] );
      R( 9 )  = alphaP - alphaPOld - xi * dLambda;
      R( 10 ) = yieldFunction( M, alphaP ) / c0;
      return R;
    }

    /**
     * Residual of the return to the apex for X = {theta = ln Je, alphaP}, with Fe = exp( theta / 3 ) I: the scaled
     * yield function of the spherical state and the hardening law of the volumetric flow.
     */
    template < typename T >
    VectorXt< T > apexResidual( const VectorXt< T >& X, const double thetaTrial, const double alphaPOld ) const
    {
      const T              theta  = X( 0 );
      const T              alphaP = X( 1 );
      const Tensor33t< T > Fe     = multiplyFastorTensorWithScalar( identity< T >(), exp( theta / 3. ) );
      const T              p      = trace( mandelStress( Fe ) ) / 3.;

      VectorXt< T > R( 2 );
      R( 0 ) = ( eta * p - xi * ( c0 + H * alphaP ) ) / c0;
      R( 1 ) = alphaP - alphaPOld - xi / etaBar * ( thetaTrial - theta );
      return R;
    }

  protected:
    const double& K;
    const double& G;
    const double& c0;
    const double& frictionAngle;
    const double& dilatancyAngle;
    const double& H;
    const double& As;
    const double& softeningModulus;
    const double& maxDamage;
    const double& nonLocalRadius;
    const double& weightingParameter;
    const double  eta, xi, etaBar;

    template < typename T >
    static Tensor33t< T > identity()
    {
      Tensor33t< T > I( T( 0.0 ) );
      for ( int i = 0; i < 3; i++ )
        I( i, i ) = T( 1.0 );
      return I;
    }

    /// the converged return mapping: the new elastic state, and the sensitivities needed for the tangents
    struct ReturnMapping {
      bool                        plastic = false;
      Tensor33d                   Fe;              ///< elastic deformation gradient
      Tensor33d                   FpNew;           ///< plastic deformation gradient
      double                      alphaP;          ///< hardening variable
      Fastor::Tensor< double, 3 > dEpPrincipal;    ///< principal plastic log strain increment (for the damage)
      double                      plasticWork;     ///< M : dEp
      Eigen::MatrixXd             dFe_dF;          ///< 9 x 9, row-major flattening of both
      Eigen::MatrixXd             dDeltaAlphaD_dF; ///< 1 x 9, sensitivity of the local damage increment
    };

    ReturnMapping returnMapping( const Tensor33d& F, const Tensor33d& FpOld, double alphaPOld ) const;
    ReturnMapping returnToCone( const Tensor33d& F, const Tensor33d& FpOld, double alphaPOld, bool& converged ) const;
    ReturnMapping returnToApex( const Tensor33d& F, const Tensor33d& FpOld, double alphaPOld ) const;

    /// the increment of the local damage variable for principal plastic log strain increments
    double deltaAlphaLocal( const Fastor::Tensor< double, 3 >& dEpPrincipal ) const;

    /// the ductility measure of CDPM2 (Grassl et al. 2013)
    double ductility( double Rs ) const;
  };

} // namespace Marmot::Materials
