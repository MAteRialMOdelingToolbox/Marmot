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
   * Both networks use the same compressible neo-Hookean potential
   * \f$\Psi(\boldsymbol C)=\frac{\mu}{2}(I_1-3-\ln\det\boldsymbol C)+\frac{\kappa}{8}
   * (\ln\det\boldsymbol C)^2\f$, evaluated on the total \f$\boldsymbol C\f$ for network
   * A and on the elastic \f$\boldsymbol C^{\rm e}\f$ for network B's spring.
   *
   * @par Material parameters
   * - @b #muA    -- network A shear modulus
   * - @b #kappaA -- network A bulk modulus
   * - @b #muB    -- network B (spring) shear modulus
   * - @b #kappaB -- network B (spring) bulk modulus
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

    /** Network A shear modulus (read from @c materialProperties[0]) */
    const double muA;
    /** Network A bulk modulus (read from @c materialProperties[1]) */
    const double kappaA;
    /** Network B shear modulus (read from @c materialProperties[2]) */
    const double muB;
    /** Network B bulk modulus (read from @c materialProperties[3]) */
    const double kappaB;
    /** Flow-rate prefactor (read from @c materialProperties[4]) */
    const double c1;
    /** Chain-stretch exponent (read from @c materialProperties[5]) */
    const double c2;
    /** Deviatoric-Mandel-stress-magnitude exponent (read from @c materialProperties[6]) */
    const double c3;

    /** Algorithm variant selector (read from @c materialProperties[7]). */
    const int implementationType;
    /** Density (read from @c materialProperties[8]) (if provided). */
    const double density;

    /**
     * @brief Construct the Bergström-Boyce model.
     * @param materialProperties Array with parameters: #muA, #kappaA, #muB, #kappaB, #c1, #c2, #c3,
     * #implementationType, #density (optional).
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
      if ( this->nMaterialProperties < 9 ) {
        throw std::runtime_error(
          std::string( MakeString() << __PRETTY_FUNCTION__ << ": No density given! nMaterialProperties < 9." ) );
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
     *   \Psi(\boldsymbol C) = \frac{\mu}{2}\left(I_1 - 3 - \ln\det\boldsymbol C\right)
     *                       + \frac{\kappa}{8}\left(\ln\det\boldsymbol C\right)^2
     * \f]
     */
    template < typename T >
    std::tuple< T, Tensor33t< T > > neoHookePotential( const Tensor33t< T >& C,
                                                       const double          mu,
                                                       const double          kappa ) const
    {
      const T I1      = trace( C );
      const T detC     = determinant( C );
      const T lnDetC   = log( detC );
      const T psi      = mu / 2. * ( I1 - 3. - lnDetC ) + kappa / 8. * lnDetC * lnDetC;

      const Tensor33t< T > I     = fastorTensorFromDoubleTensor< T >( Spatial3D::I );
      const Tensor33t< T > CInv  = inverse( C );
      const T              coeff = kappa / 4. * lnDetC - mu / 2.;

      const Tensor33t< T > dPsi_dC = multiplyFastorTensorWithScalar( I, T( mu / 2. ) ) +
                                     multiplyFastorTensorWithScalar( CInv, coeff );

      return { psi, dPsi_dC };
    }

    /**
     * @brief Deviatoric flow direction, flow rate and chain stretch for network B, at a given elastic
     * deformation gradient (templated scalar type).
     * @tparam T Scalar type (double, complex).
     * @param Fe Elastic deformation gradient of network B.
     * @return \f$\{\,\boldsymbol N,\;\rho,\;\dot\gamma\,\}\f$.
     */
    template < typename T >
    std::tuple< Tensor33t< T >, T, T > computeFlowQuantities( const Tensor33t< T >& Fe ) const
    {
      const Tensor33t< T > Ce         = ContinuumMechanics::DeformationMeasures::rightCauchyGreen( Fe );
      const T              I1         = trace( Ce );
      const T              lambdaChain = sqrt( I1 / T( 3.0 ) );

      Tensor33t< T > devCe   = deviatoric( Ce );
      T              normDev = sqrt( Fastor::inner( devCe, devCe ) );
      if ( Math::makeReal( normDev ) == 0.0 )
        normDev += 1e-15;

      const Tensor33t< T > N   = multiplyFastorTensorWithScalar( devCe, T( 1.0 ) / normDev );
      T                    rho = T( muB ) * normDev;
      if ( Math::makeReal( rho ) == 0.0 )
        rho += 1e-15;

      const T stretchTerm        = lambdaChain - T( 1.0 );
      const T stretchTermClamped = Math::makeReal( stretchTerm ) > 0.0 ? stretchTerm : T( 0.0 );

      const T gammaDot = T( c1 ) * pow( stretchTermClamped, c2 ) * pow( rho, c3 );

      return { N, rho, gammaDot };
    }

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
