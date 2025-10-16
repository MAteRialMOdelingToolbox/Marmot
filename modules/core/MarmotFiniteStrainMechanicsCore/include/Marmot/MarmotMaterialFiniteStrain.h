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
#include "Fastor/Fastor.h"
#include "Marmot/MarmotMaterial.h"
#include <Fastor/tensor/Tensor.h>

class MarmotMaterialFiniteStrain : public MarmotMaterial {

  /**
   * @class MarmotMaterialFiniteStrain
   * @brief Abstract basic class for mechanical materials in the finite strain regime
   */
public:
  /**
   * @struct ConstitutiveResponse
   * @brief Constitutive response of a material at given state.
   * @tparam nDim Number of spatial dimensions (2 or 3).
   *
   *  Contains stress, density and elastic energy density.
   */
  template < int nDim >
  struct ConstitutiveResponse {
    Fastor::Tensor< double, nDim, nDim > tau;                  ///< Kirchhoff stress
    double                               rho;                  ///< mass density
    double                               elasticEnergyDensity; ///< elastic energy per unit volume
  };

  /**
   * @struct AlgorithmicModuli
   * @brief Algorithmic tangent moduli of a material.
   * @tparam nDim Number of spatial dimensions (2 or 3).
   *
   * Contains the algorithmic tangent moduli \f$\frac{\partial \boldsymbol{\tau}}{\partial \boldsymbol{F}}\f$
   * with respect to the deformation gradient \f$\boldsymbol{F}\f$.
   * */
  template < int nDim >
  struct AlgorithmicModuli {
    Fastor::Tensor< double, nDim, nDim, nDim, nDim > dTau_dF; ///< tangent operator w.r.t. deformation gradient
  };

  /**
   * @struct Deformation
   * @brief Represents the deformation state of a material.
   *
   * This struct holds the deformation gradient \f$\boldsymbol{F}\f$ which describes the local deformation of a
   * material.
   */
  template < int nDim >
  struct Deformation {
    Fastor::Tensor< double, nDim, nDim > F; ///< deformation gradient
  };

  /**
   * @struct TimeIncrement
   * @brief Represents a time increment in a simulation.
   *
   * This struct holds information about the current time
   * and the time step size (dT).
   */
  struct TimeIncrement {
    const double time; ///< time at the beginning of the increment
    const double dT;   ///< size of the time increment
  };

  using MarmotMaterial::MarmotMaterial;

  /**
   * @brief Updates the material state.
   * @param[inout] response ConstitutiveResponse instance
   * @param[out] tangents AlgorithmicModuli instance
   * @param[in] deformation Deformation instance
   * @param[in] timeIncrement TimeIncrement instance
   *
   * Computes the Kirchhoff \f$\boldsymbol{\tau}\f$ stress
   * from the deformation gradient \f$\boldsymbol{F}\f$ and the
   * time increment \f$\Delta t\f$.
   * It further updates the mass density \f$\rho\f$ and the elastic energy density.
   * Additionally, computes the algorithmic tangent moduli
   * \f$\frac{\partial \boldsymbol{\tau}}{\partial \boldsymbol{F}}\f$.
   *
   * */
  virtual void computeStress( ConstitutiveResponse< 3 >& response,
                              AlgorithmicModuli< 3 >&    tangents,
                              const Deformation< 3 >&,
                              const TimeIncrement& ) = 0;
  /**
   * @brief Computes the Kirchhoff stress given the deformation, time increment, and eigen deformation.
   * @param[inout] response ConstitutiveResponse instance
   * @param[out] tangents AlgorithmicModuli instance
   * @param[in] deformation Deformation instance
   * @param[in] timeIncrement TimeIncrement instance
   * @param[in] eigenDeformation Tuple representing eigen deformation in each spatial direction.
   *
   */
  virtual void computeStress( ConstitutiveResponse< 3 >&                  response,
                              AlgorithmicModuli< 3 >&                     tangents,
                              const Deformation< 3 >&                     deformation,
                              const TimeIncrement&                        timeIncrement,
                              const std::tuple< double, double, double >& eigenDeformation );

  /**
   * @brief Compute stress under plane strain conditions.
   * @param[inout] response ConstitutiveResponse instance
   * @param[out] algorithmicModuli AlgorithmicModuli instance
   * @param[in] deformation Deformation instance
   * @param[in] timeIncrement TimeIncrement instance
   *
   * It uses the general 3D computeStress function for a plane strain Deformation.
   * The algorithmic tangent is modified according to plane strain conditions.
   */
  virtual void computePlaneStrain( ConstitutiveResponse< 3 >& response,
                                   AlgorithmicModuli< 3 >&    algorithmicModuli,
                                   const Deformation< 3 >&    deformation,
                                   const TimeIncrement&       timeIncrement );
  /**
   * @brief Compute stress under plane strain conditions with eigen deformation.
   * @param[inout] response ConstitutiveResponse instance
   * @param[out] algorithmicModuli AlgorithmicModuli instance
   * @param[in] deformation Deformation instance
   * @param[in] timeIncrement TimeIncrement instance
   * @param[in] eigenDeformation Tuple representing eigen deformation in each spatial direction.
   *
   * It uses the general 3D computeStress function for a plane strain Deformation.
   * The algorithmic tangent is modified according to plane strain conditions.
   */
  virtual void computePlaneStrain( ConstitutiveResponse< 3 >&                  response,
                                   AlgorithmicModuli< 3 >&                     algorithmicModuli,
                                   const Deformation< 3 >&                     deformation,
                                   const TimeIncrement&                        timeIncrement,
                                   const std::tuple< double, double, double >& eigenDeformation );
  /**
   * @brief Compute stress under plane stress conditions.
   * @param[inout] response ConstitutiveResponse instance
   * @param[out] algorithmicModuli AlgorithmicModuli instance
   * @param[in] deformation Deformation instance
   * @param[in] timeIncrement TimeIncrement instance
   *
   * It uses the general 3D computeStress function and iteratively finds the out-of-plane deformation.
   * The algorithmic tangent is modified according to plane stress conditions.
   */
  virtual void computePlaneStress( ConstitutiveResponse< 2 >& response,
                                   AlgorithmicModuli< 2 >&    algorithmicModuli,
                                   const Deformation< 2 >&    deformation,
                                   const TimeIncrement&       timeIncrement );

  /**
   * @brief Find the eigen deformation that corresponds to a given eigen stress.
   * @param initialGuess Initial guess for the eigen deformation.
   * @param eigenStress Target eigen stress.
   * @return Eigen deformation that corresponds to the given eigen stress.
   *
   * This function iteratively finds the eigen deformation that corresponds to a given eigen stress.
   * This is used e.g. for geostatic stress initialization.
   */
  std::tuple< double, double, double > findEigenDeformationForEigenStress(
    const std::tuple< double, double, double >& initialGuess,
    const std::tuple< double, double, double >& eigenStress );
};
