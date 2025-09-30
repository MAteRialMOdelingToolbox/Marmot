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

  /*
     Abstract basic class for mechanical materials in the finite strain regime
  */
public:
  /** @struct ConstitutiveResponse
   * @brief Constitutive response of a material at given state.
   *
   *  Contains stress, density and elastic energy density.
   * */
  template < int nDim >
  struct ConstitutiveResponse {
    Fastor::Tensor< double, nDim, nDim > tau;                  // kirchhoff stress
    double                               rho;                  // density
    double                               elasticEnergyDensity; // elastic energy per unit volume
  };

  template < int nDim >
  struct AlgorithmicModuli {
    Fastor::Tensor< double, nDim, nDim, nDim, nDim > dTau_dF; // tangent operator w.r.t. deformation gradient
  };

  /** @struct Deformation
   * @brief Represents the deformation state of a material.
   *
   * This struct holds the deformation gradient tensor (F) which describes the local deformation of a material.
   */
  template < int nDim >
  struct Deformation {
    Fastor::Tensor< double, nDim, nDim > F;
  };

  /** @struct TimeIncrement
   * @brief Represents a time increment in a simulation.
   *
   * This struct holds information about the current time and the time step size (dT).
   */
  struct TimeIncrement {
    const double time;
    const double dT;
  };

  using MarmotMaterial::MarmotMaterial;

  /** @brief Computes the Kirchhoff stress given the deformation and time increment.
   * @param[inout] response Constitutive response of the material
   * @param[out] tangents Algorithmic moduli
   * @param[in] deformation Deformation state
   * @param[in] timeIncrement Time increment
   *
   *
   *
   *  Pure virtual function, must be implemented by derived classes.
   * */
  virtual void computeStress( ConstitutiveResponse< 3 >& response,
                              AlgorithmicModuli< 3 >&    tangents,
                              const Deformation< 3 >&,
                              const TimeIncrement& ) = 0;

  /**
   * Compute Stress, but account for eigen deformations (e.g, geostatic stress states). Modifies algorithmic tangent
   * accordingly.*/
  virtual void computeStress( ConstitutiveResponse< 3 >&                  response,
                              AlgorithmicModuli< 3 >&                     tangents,
                              const Deformation< 3 >&                     deformation,
                              const TimeIncrement&                        timeIncrement,
                              const std::tuple< double, double, double >& eigenDeformation );

  virtual void computePlaneStrain( ConstitutiveResponse< 3 >& response,
                                   AlgorithmicModuli< 3 >&    algorithmicModuli,
                                   const Deformation< 3 >&    deformation,
                                   const TimeIncrement&       timeIncrement );
  /***/
  /* * Compute stress under plane strain conditions, but account for eigen deformations (e.g, geostatic stress
   * states).*/
  /* * Modifies algorithmic tangent accordingly.*/
  virtual void computePlaneStrain( ConstitutiveResponse< 3 >&                  response,
                                   AlgorithmicModuli< 3 >&                     algorithmicModuli,
                                   const Deformation< 3 >&                     deformation,
                                   const TimeIncrement&                        timeIncrement,
                                   const std::tuple< double, double, double >& eigenDeformation );

  virtual void computePlaneStress( ConstitutiveResponse< 2 >& response,
                                   AlgorithmicModuli< 2 >&    algorithmicModuli,
                                   const Deformation< 2 >&    deformation,
                                   const TimeIncrement&       timeIncrement );
  /** Find for given eigen stress the appropriate eigen deformation.
   * */
  std::tuple< double, double, double > findEigenDeformationForEigenStress(
    const std::tuple< double, double, double >& initialGuess,
    const std::tuple< double, double, double >& eigenStress );
};
