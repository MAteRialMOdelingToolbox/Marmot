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
 * Matthias Neuner matthias.neuner@uibk.ac.at
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
#include "Marmot/MarmotAutomaticDifferentiationForFastor.h"
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include <Eigen/Core>
#include <functional>

/**
 * @brief Base class for finite strain materials using automatic differentiation.
 *
 * This class extends MarmotMaterialFiniteStrain by providing an interface for
 * computing stresses and algorithmic tangent moduli using automatic
 * differentiation. Derived material models implement computeStressAD() to
 * supply the AD-based constitutive response.
 */
class MarmotMaterialFiniteStrainAD : public MarmotMaterialFiniteStrain {

public:
  MarmotMaterialFiniteStrainAD( const double* matProperties_, int nMaterialProperties_, int materialNumber_ )
    : MarmotMaterialFiniteStrain( matProperties_, nMaterialProperties_, materialNumber_ )
  {
  }
  /**
   * @struct ConstitutiveResponseAD
   * @brief Constitutive response of a material at given state.
   * @tparam nDim Number of spatial dimensions (2 or 3).
   *
   *  Contains stress, density and elastic energy density.
   */
  template < int nDim >
  struct ConstitutiveResponseAD {
    Fastor::Tensor< autodiff::dual, nDim, nDim > tau;                  ///< Kirchhoff stress
    double                                       elasticEnergyDensity; ///< elastic energy per unit volume
    double*                                      stateVars;            ///< pointer to state variables
  };

  /**
   * @struct DeformationAD
   * @brief Represents the deformation state of a material.
   *
   * This struct holds the deformation gradient \f$\boldsymbol{F}\f$ which describes the local deformation of a
   * material.
   */
  template < int nDim >
  struct DeformationAD {
    Fastor::Tensor< autodiff::dual, nDim, nDim > F; ///< deformation gradient
  };

  /**
   * @brief Compute Kirchhoff stress and algorithmic tangent using automatic differentiation.
   *
   * This method overrides MarmotMaterialFiniteStrain::computeStress() and provides a generic
   * AD-based implementation for finite strain material models.
   *
   * The child class must implement computeStressAD(), which evaluates the constitutive response
   * (Kirchhoff stress and any state updates) for a dual-valued deformation gradient. This base
   * implementation then uses Marmot::AutomaticDifferentiation::dF_dT() to:
   *  - compute the Kirchhoff stress \f$\boldsymbol{\tau}\f$ for the given deformation, and
   *  - compute the consistent algorithmic tangent \f$\partial \boldsymbol{\tau} / \partial \boldsymbol{F}\f$.
   *
   * @note Automatic differentiation requires multiple seeded evaluations of computeStressAD().
   *       Therefore, the state variables are reset to their "old" values before each seeded call.
   *       Implementations of computeStressAD() must update state based only on primal values of
   *       the dual inputs and must not branch or accumulate history based on derivative components.
   *
   * @param[out] response     Constitutive response container (stress, density, energy, state variables).
   * @param[out] tangents     Algorithmic moduli container; filled with \f$\partial \tau / \partial F\f$.
   * @param[in]  deformation  Deformation state (deformation gradient \f$\boldsymbol{F}\f$).
   * @param[in]  timeIncrement Time increment information for rate-/history-dependent models.
   */
  void computeStress( ConstitutiveResponse< 3 >& response,
                      AlgorithmicModuli< 3 >&    tangents,
                      const Deformation< 3 >&    deformation,
                      const TimeIncrement&       timeIncrement ) const override
  {
    using namespace Marmot;
    using namespace FastorStandardTensors;
    using scalar = autodiff::dual;

    Eigen::Map< Eigen::VectorXd > stateVars( response.stateVars, this->getNumberOfRequiredStateVars() );

    // Remember old state to reset during potential multiple AD evaluations
    const Eigen::VectorXd stateVarsOld = stateVars;
    const Tensor33d       tauOld       = response.tau;

    std::function< Tensor33t< scalar >( const Tensor33t< scalar >& ) > computeTauAD =
      [&]( const Tensor33t< scalar >& F_ ) {
        // Reset stateVars to old state
        stateVars = stateVarsOld;

        // Explicitly convert old (double-valued) stress to a dual-valued tensor for AD
        Tensor33t< scalar > tauAD = makeDual( tauOld );

        // Construct AD state
        DeformationAD< 3 > deformationAD;
        deformationAD.F = F_;

        ConstitutiveResponseAD< 3 > responseAD;
        responseAD.tau                  = tauAD;
        responseAD.elasticEnergyDensity = response.elasticEnergyDensity;
        responseAD.stateVars            = stateVars.data();

        // Compute stress utilizing the child class's AD implementation
        this->computeStressAD( responseAD, deformationAD, timeIncrement );

        // Extract and return the dual stress tensor for differentiation
        Tensor33t< scalar > res = responseAD.tau;

        return res;
      };

    // Compute Kirchhoff stress (tau) and algorithmic tangent (dTau_dF) with autodiff
    std::tie( response.tau, tangents.dTau_dF ) = Marmot::AutomaticDifferentiation::dF_dT( computeTauAD, deformation.F );
  }

  /**
   * @brief Compute the constitutive response in automatic differentiation (AD) form.
   *
   * This pure virtual function must be implemented by derived material models to compute
   * the Kirchhoff stress tensor in \p response.tau for a given deformation in dual number
   * form. The base class uses this AD formulation together with
   * Marmot::AutomaticDifferentiation utilities to obtain algorithmic tangents by
   * differentiating the stress with respect to the deformation gradient.
   *
   * Implementations may update the provided state variables via \p response.stateVars
   * according to the given time increment.
   *
   * @param[out] response    Constitutive response in AD form (3D), with \p response.tau
   *                         to be filled with the computed Kirchhoff stress.
   * @param[in]  deformation Deformation state in AD form (3D), providing the deformation
   *                         gradient \f$\boldsymbol{F}\f$ as dual numbers.
   * @param[in]  timeIncrement Current time increment information for the constitutive update.
   */
  virtual void computeStressAD( ConstitutiveResponseAD< 3 >& response,
                                const DeformationAD< 3 >&    deformation,
                                const TimeIncrement&         timeIncrement ) const = 0;
};
