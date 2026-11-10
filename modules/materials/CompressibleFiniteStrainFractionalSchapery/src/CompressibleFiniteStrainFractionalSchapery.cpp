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

#include "Marmot/CompressibleFiniteStrainFractionalSchapery.h"
#include "Marmot/MarmotAutomaticDifferentiationForFastor.h"
#include "Marmot/MarmotDeformationMeasures.h"
#include "Marmot/MarmotEnergyDensityFunctions.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotStressMeasures.h"
#include <Fastor/expressions/linalg_ops/unary_trans_op.h>
#include <Fastor/tensor/Tensor.h>
#include <Fastor/tensor_algebra/einsum_explicit.h>
#include <Fastor/tensor_algebra/indicial.h>
#include <algorithm>
#include <autodiff/forward/dual/dual.hpp>
#include <cmath>
#include <cstring>
#include <functional>

namespace Marmot::Materials {

  using namespace Fastor;
  using namespace Marmot::Math::Tensors::FastorIndices;
  using namespace Marmot::Math::Tensors::FastorTensors;

  CompressibleFiniteStrainFractionalSchapery::CompressibleFiniteStrainFractionalSchapery(
    const double* materialProperties,
    int           nMaterialProperties,
    int           materialLabel )
    : MaterialFiniteStrain( materialProperties, nMaterialProperties, materialLabel ),
      hyperelasticBase( static_cast< HyperelasticBase >( materialProperties[0] ) ),
      onlyShearCreep( materialProperties[1] ),
      elasticProperties( &materialProperties[2], nElasticPropertiesMap.at( hyperelasticBase ) ),
      alpha( materialProperties[2 + nElasticPropertiesMap.at( hyperelasticBase )] ),
      tau0( materialProperties[3 + nElasticPropertiesMap.at( hyperelasticBase )] ),
      g_inf( materialProperties[4 + nElasticPropertiesMap.at( hyperelasticBase )] ),
      a1( materialProperties[5 + nElasticPropertiesMap.at( hyperelasticBase )] ),
      a2( materialProperties[6 + nElasticPropertiesMap.at( hyperelasticBase )] ),
      b_S( materialProperties[7 + nElasticPropertiesMap.at( hyperelasticBase )] )
  {
    initializeFractionalSpectrum();
    initializeStateLayout();
  }

  void CompressibleFiniteStrainFractionalSchapery::initializeFractionalSpectrum()
  {
    // Discretize the continuous fractional power-law spectrum over 5 logarithmic decades around tau0
    const double alpha_safe = std::clamp( alpha, 0.01, 0.99 );
    const double tau0_safe  = std::max( tau0, 1e-4 );
    const double g_inf_safe = std::clamp( g_inf, 0.01, 0.99 );

    // Spanning 5 decades: 1e-2 * tau0, 1e-1 * tau0, tau0, 1e1 * tau0, 1e2 * tau0
    const double factors[nFractionalModes] = { 0.01, 0.1, 1.0, 10.0, 100.0 };
    double       sumWeight                 = 0.0;
    double       weights[nFractionalModes];

    for ( int m = 0; m < nFractionalModes; ++m ) {
      tau[m] = std::max( tau0_safe * factors[m], 1e-6 );
      // Power-law spectral distribution: weight_m ~ (tau_m / tau0)^(-alpha)
      weights[m] = std::pow( factors[m], -alpha_safe );
      sumWeight += weights[m];
    }

    const double sumGamma = 1.0 - g_inf_safe;
    for ( int m = 0; m < nFractionalModes; ++m ) {
      gamma[m] = sumGamma * ( weights[m] / sumWeight );
    }
  }

  double CompressibleFiniteStrainFractionalSchapery::getDensity( const double* stateVars ) const
  {
    (void)stateVars;
    const int nRequired = 2 + nElasticPropertiesMap.at( hyperelasticBase ) + 6;
    if ( nMaterialProperties > nRequired ) {
      return materialProperties[nMaterialProperties - 1];
    }
    return 0.0;
  }

  std::tuple< double, Tensor33d, Tensor3333d, Tensor333333d > CompressibleFiniteStrainFractionalSchapery::
    computeEnergyDensityAndDerivatives( const Tensor33d& C ) const
  {
    using namespace Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions;

    std::function< autodiff::dual3rd( const Tensor33t< autodiff::dual3rd >& ) > energyDensityFunction;

    switch ( hyperelasticBase ) {
    case NeoHooke:
      energyDensityFunction = [this]( const Tensor33t< autodiff::dual3rd >& C_ad ) {
        return NeoHookePotential< autodiff::dual3rd >( C_ad, elasticProperties[0] ) +
               VolumetricPenaltyPotential< autodiff::dual3rd >( C_ad, elasticProperties[1] );
      };
      break;
    case Yeoh:
      energyDensityFunction = [this]( const Tensor33t< autodiff::dual3rd >& C_ad ) {
        return YeohPotential< autodiff::dual3rd >( C_ad,
                                                   elasticProperties[0],
                                                   elasticProperties[1],
                                                   elasticProperties[2] ) +
               VolumetricPenaltyPotential< autodiff::dual3rd >( C_ad, elasticProperties[3] );
      };
      break;
    case PenceGouNeoHooke:
      energyDensityFunction = [this]( const Tensor33t< autodiff::dual3rd >& C_ad ) {
        return PenceGouPotentialB< autodiff::dual3rd >( C_ad, elasticProperties[0], elasticProperties[1] );
      };
      break;
    case MooneyRivlin:
      energyDensityFunction = [this]( const Tensor33t< autodiff::dual3rd >& C_ad ) {
        return MooneyRivlinPotential< autodiff::dual3rd >( C_ad, elasticProperties[0], elasticProperties[1] ) +
               VolumetricPenaltyPotential< autodiff::dual3rd >( C_ad, elasticProperties[2] );
      };
      break;
    case ArrudaBoyce:
      energyDensityFunction = [this]( const Tensor33t< autodiff::dual3rd >& C_ad ) {
        return ArrudaBoyce8ChainPotential< autodiff::dual3rd >( C_ad, elasticProperties[0], elasticProperties[1] ) +
               VolumetricPenaltyPotential< autodiff::dual3rd >( C_ad, elasticProperties[2] );
      };
      break;
    default:
      throw std::runtime_error(
        "CompressibleFiniteStrainFractionalSchapery::computeEnergyDensityAndDerivatives: Unknown hyperelastic base." );
    }

    return Marmot::Math::Differentiation::Automatic::ThirdOrder::d3f_dT3< 3 >( energyDensityFunction, C, true );
  }

  void CompressibleFiniteStrainFractionalSchapery::computeStress( ConstitutiveResponse< 3 >& response,
                                                                  AlgorithmicModuli< 3 >&    tangents,
                                                                  const Deformation< 3 >&    deformation,
                                                                  const TimeIncrement&       timeIncrement ) const
  {
    const auto& F = deformation.F;

    using namespace Marmot::ContinuumMechanics;
    const auto [C, dC_dF] = Kinematics::DeformationMeasures::FirstOrder::rightCauchyGreen( F );

    const auto [psi_, dPsi_dC, d2Psi_dCdC, d3Psi_dCdCdC] = computeEnergyDensityAndDerivatives( C );

    // Instantaneous PK2 stress
    const Tensor33d PK2zero = 2. * dPsi_dC;

    // Split into volumetric and isochoric/deviatoric parts
    const Tensor33d PK2vol = ( onlyShearCreep / 3.0 ) * trace( PK2zero ) * Spatial3D::I;
    const Tensor33d PK2iso = PK2zero - PK2vol;

    // First invariant of unimodular C
    const double detC        = determinant( C );
    const double J23         = std::cbrt( std::max( detC, 1e-12 ) );
    const double Ibar1       = trace( C ) / J23;
    const double trC_minus_3 = std::max( Ibar1 - 3.0, 0.0 );

    // Schapery nonlinear strain functions
    const double g1 = std::max( 1.0 + a1 * trC_minus_3, 0.01 );
    const double g2 = std::max( 1.0 + a2 * trC_minus_3, 0.01 );
    const double aS = std::clamp( std::exp( -b_S * trC_minus_3 ), 1e-6, 1e6 );

    // Reduced time increment
    const double d_psi = ( timeIncrement.dT > 0.0 ) ? ( timeIncrement.dT / aS ) : 0.0;

    // Effective driving stress
    const Tensor33d S_eff = g2 * PK2iso;

    // Old effective stress from state variables
    TensorMap33d    S_eff_old = stateLayout.getAs< TensorMap33d >( response.stateVars, "S0_old" );
    const Tensor33d dS_eff    = S_eff - S_eff_old;
    std::memcpy( S_eff_old.data(), S_eff.data(), 9 * sizeof( double ) );

    double* creepStateVars = stateLayout.getPtr( response.stateVars, "creepStateVars" );

    Tensor33d S_neq;
    S_neq.zeros();
    double sumBeta = 0.0;

    for ( int m = 0; m < nFractionalModes; ++m ) {
      const Tensor33d Q_m_old( creepStateVars + m * 9 );
      const double&   tau_m   = tau[m];
      const double&   gamma_m = gamma[m];

      const double dT_tau    = std::max( d_psi / tau_m, 1e-15 );
      const double expFactor = std::exp( -dT_tau );

      double alpha_fact = expFactor;
      double beta_fact  = gamma_m / dT_tau * ( 1.0 - expFactor );

      if ( dT_tau < 1e-6 ) {
        alpha_fact = 1.0 - dT_tau + 0.5 * dT_tau * dT_tau;
        beta_fact  = gamma_m * ( 1.0 - 0.5 * dT_tau + ( 1.0 / 6.0 ) * dT_tau * dT_tau );
      }

      const Tensor33d Q_m_new = alpha_fact * Q_m_old + beta_fact * dS_eff;
      S_neq += Q_m_new;
      sumBeta += beta_fact;

      std::memcpy( creepStateVars + m * 9, Q_m_new.data(), 9 * sizeof( double ) );
    }

    // Assemble total PK2 stress: volumetric + equilibrium isochoric + nonequilibrium isochoric
    const double    g_inf_safe   = std::clamp( g_inf, 0.01, 0.99 );
    const Tensor33d PK2iso_total = g_inf_safe * PK2iso + g1 * S_neq;
    const Tensor33d PK2          = evaluate( PK2vol + PK2iso_total );

    // Push forward to Kirchhoff stress
    const auto [tau_k, dTau_dPK2, dTau_dF] = StressMeasures::FirstOrder::KirchhoffStressFromPK2( PK2, F );
    response.tau                           = tau_k;
    response.elasticEnergyDensity          = psi_;

    // Algorithmic tangent operator dPK2 / dE
    const Tensor3333d dPK2_dE    = 4. * d2Psi_dCdC;
    const Tensor3333d dPK2vol_dE = onlyShearCreep / 3.0 *
                                   outer( Spatial3D::I, einsum< ijkl, kl >( dPK2_dE, Spatial3D::I ) );
    const Tensor3333d dPK2iso_dE = dPK2_dE - dPK2vol_dE;

    // Effective modulus multiplier: g_inf + g1 * g2 * sumBeta
    const double      beta_eff      = g_inf_safe + g1 * g2 * sumBeta;
    const Tensor3333d dPK2_total_dE = evaluate( dPK2vol_dE + beta_eff * dPK2iso_dE );

    // Algorithmic tangent operator dtau / dF
    tangents.dTau_dF = einsum< ijKL, KLMN >( einsum< ijKL, KLMN >( dTau_dPK2, dPK2_total_dE ), 0.5 * dC_dF ) + dTau_dF;
  }

} // namespace Marmot::Materials
