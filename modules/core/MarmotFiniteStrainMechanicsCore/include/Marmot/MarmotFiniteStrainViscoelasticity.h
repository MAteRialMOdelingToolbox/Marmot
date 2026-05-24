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
 * Alexander Dummer alexander.dummer@uibk.ac.at
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
#include "Marmot/MarmotMath.h"

namespace Marmot {
  namespace ContinuumMechanics::FiniteStrain::Viscoelasticity {

    /**
     * @brief Struct to hold properties of a generalized maxwell model
     */
    struct MaxwellProperties {
      int                   nMaxwell; ///< Number of maxwell elements
      std::vector< double > gamma;    ///< Weighting factors of maxwell elements
      double                sumGamma; ///< Sum of weighting factors
      std::vector< double > tau;      ///< Relaxation times of maxwell elements
    };

    MaxwellProperties createMaxwellProperties( int nMaxwell, const double* gammaTauPairVector );

    /**
     * @brief Evaluate generalized maxwell model contribution to stress and tangent
     * @param[in,out] stress             Stress tensor to be updated (input is the instantaneous hyperelastic stress)
     * @param[in,out] tangent            Tangent tensor to be updated (input is the instantaneous hyperelasticelastic
     * tangent)
     * @param[in]  dTangent_dDeformation Derivative of tangent w.r.t. deformation
     * @param[in]  initialCompliance     Initial compliance tensor
     * @param[in]  dStress               Increment of the initial stress tensor
     * @param[in]  dT                    Time increment
     * @param[in]  maxwellProperties     Properties of the generalized maxwell model
     * @param[in,out] stateVars          State variables array (size: nMaxwell * 9)
     *
     * @details This update function uses the approach in:
     * Liu et al. 2021, A continuum and computational framework for viscoelastodynamics: I. Finite deformation linear
     * models, Comput. Meth. Appl. Mech. Engrg. 385, 114059
     */
    void evaluateGeneralizedMaxwellModel( FastorStandardTensors::Tensor33d&           stress,
                                          FastorStandardTensors::Tensor3333d&         tangent,
                                          const FastorStandardTensors::Tensor333333d& dTangent_dDeformation,
                                          const FastorStandardTensors::Tensor3333d&   initialCompliance,
                                          const FastorStandardTensors::Tensor33d&     dStress,
                                          const double                                dT,
                                          const MaxwellProperties&                    maxwellProperties,
                                          double*                                     stateVars );

    /**
     * @brief Evaluate generalized maxwell model contribution to stress and tangent
     * @param[in,out] stress             Stress tensor to be updated (input is the instantaneous hyperelastic stress)
     * tangent)
     * @param[in,out] tangent            Tangent tensor to be updated (input is the instantaneous hyperelasticelastic
     * tangent)
     * @param[in]  dStress               Increment of the initial stress tensor
     * @param[in]  dT                    Time increment
     * @param[in]  maxwellProperties     Properties of the generalized maxwell model
     * @param[in,out] stateVars          State variables array (size: nMaxwell * 9)
     *
     * @details This update correspond to the approach by Simo (1987)
     * and is used in the implementation of the generalized maxwell model in the finite element code
     */
    void evaluateGeneralizedMaxwellModel( FastorStandardTensors::Tensor33d&       stress,
                                          FastorStandardTensors::Tensor3333d&     tangent,
                                          const FastorStandardTensors::Tensor33d& dStress,
                                          const double                            dT,
                                          const MaxwellProperties&                maxwellProperties,
                                          double*                                 stateVars );

    /**
     * @brief Evaluate generalized maxwell model contribution to stress without tangent update
     * @param[in,out] stress             Stress tensor to be updated (input is the instantaneous hyperelastic stress)
     * tangent)
     * @param[in]  dStress               Increment of the initial stress tensor
     * @param[in]  dT                    Time increment
     * @param[in]  maxwellProperties     Properties of the generalized maxwell model
     * @param[in,out] stateVars          State variables array (size: nMaxwell * 9)
     *
     * @details This update correspond to the approach by Simo (1987)
     * and is used in the implementation of the generalized maxwell model in the finite element code
     */
    template < typename T = double >
    void evaluateGeneralizedMaxwellModel( FastorStandardTensors::Tensor33t< T >&       stress,
                                          const FastorStandardTensors::Tensor33t< T >& dStress,
                                          const double                                 dT,
                                          const MaxwellProperties&                     maxwellProperties,
                                          double*                                      stateVars )
    {

      if ( maxwellProperties.nMaxwell == 0 )
        return;

      using namespace Fastor;
      using namespace Marmot::FastorStandardTensors;
      using namespace Marmot::FastorIndices;

      // scale equilibrium stress contribution
      stress = multiplyFastorTensorWithScalar( stress, T( 1.0 - maxwellProperties.sumGamma ) );

      for ( int i = 0; i < maxwellProperties.nMaxwell; ++i ) {
        // get old  maxewell element stress from state variables
        const Tensor33d& Q_n = Tensor33d( stateVars + i * 9 );

        // get parameters of maxwell element
        const double& tau   = maxwellProperties.tau[i];
        const double& gamma = maxwellProperties.gamma[i];

        const double dT_tau    = std::max( dT / tau, 1e-15 );
        const double expFactor = Math::exp( -dT_tau );

        double alpha = expFactor;
        double beta  = gamma / dT_tau * ( 1.0 - expFactor );

        if ( dT_tau < 1e-6 ) {
          // use taylor expansion for small dt/tau
          alpha = 1.0 - dT_tau + 0.5 * dT_tau * dT_tau;
          beta  = gamma * ( 1.0 - 0.5 * dT_tau + 1.0 / 6.0 * dT_tau * dT_tau );
        }

        // compute new stress in maxwell element
        const Tensor33t< T > Q_np = makeOtherScalarType< T >( evaluate( alpha * Q_n ) ) +
                                    multiplyFastorTensorWithScalar( dStress, T( beta ) );

        // add contribution to stress
        stress += Q_np;

        const Tensor33d Q_np_real = makeReal( Q_np );

        // update state variables
        memcpy( stateVars + i * 9, Q_np_real.data(), 9 * sizeof( double ) );
      }
    }

    /**
     * @brief Evaluate generalized maxwell model contribution to stress without tangent update
     * @param[in,out] stress             Stress tensor to be updated (input is the instantaneous hyperelastic stress)
     * tangent)
     * @param[in]  tangent               Tangent tensor to be used for stress update
     * @param[in]  initialCompliance     Initial compliance tensor
     * @param[in]  dStress               Increment of the initial stress tensor
     * @param[in]  dT                    Time increment
     * @param[in]  maxwellProperties     Properties of the generalized maxwell model
     * @param[in,out] stateVars          State variables array (size: nMaxwell * 9)
     *
     * @details This update function uses the approach in:
     * Liu et al. 2021, A continuum and computational framework for viscoelastodynamics: I. Finite deformation linear
     * models, Comput. Meth. Appl. Mech. Engrg. 385, 114059
     */

    template < typename T = double >
    void evaluateGeneralizedMaxwellModel( FastorStandardTensors::Tensor33t< T >&         stress,
                                          const FastorStandardTensors::Tensor3333t< T >& tangent,
                                          const FastorStandardTensors::Tensor3333t< T >& initialCompliance,
                                          const FastorStandardTensors::Tensor33t< T >&   dStress,
                                          const double                                   dT,
                                          const MaxwellProperties&                       maxwellProperties,
                                          double*                                        stateVars )
    {

      if ( maxwellProperties.nMaxwell == 0 )
        return;

      using namespace Fastor;
      using namespace Marmot::FastorStandardTensors;
      using namespace Marmot::FastorIndices;

      // scale equilibrium stress contribution
      stress = stress * ( 1.0 - maxwellProperties.sumGamma );

      for ( int i = 0; i < maxwellProperties.nMaxwell; ++i ) {
        // get old  maxewell element stress from state variables
        const Tensor33d& Q_n = Tensor33d( stateVars + i * 9 );

        // get parameters of maxwell element
        const double& tau   = maxwellProperties.tau[i];
        const double& gamma = maxwellProperties.gamma[i];

        const double dT_tau    = std::max( dT / tau, 1e-15 );
        const double expFactor = Math::exp( -dT_tau );

        double alpha = expFactor;
        double beta  = gamma / dT_tau * ( 1.0 - expFactor );

        if ( dT_tau < 1e-6 ) {
          // use taylor expansion for small dt/tau
          alpha = 1.0 - dT_tau + 0.5 * dT_tau * dT_tau;
          beta  = gamma * ( 1.0 - 0.5 * dT_tau + 1.0 / 6.0 * dT_tau * dT_tau );
        }

        // compute new stress in maxwell element
        const Tensor33t< T > Q_np = alpha * makeOtherScalarType< T >( Q_n ) + beta * dStress;

        // add contribution to stress
        const Tensor33t< T > H_np = einsum< ijkl, kl >( initialCompliance, Q_np );
        stress += einsum< ij, ijkl >( H_np, tangent );

        // update state variables
        Tensor33d Q_np_real = makeReal( Q_np );
        memcpy( stateVars + i * 9, Q_np_real.data(), 9 * sizeof( double ) );
      }
    }

  } // namespace ContinuumMechanics::FiniteStrain::Viscoelasticity
} // namespace Marmot
