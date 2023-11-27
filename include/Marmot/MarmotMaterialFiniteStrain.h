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
  template < int nDim >
  struct ConstitutiveResponse {
    Fastor::Tensor< double, nDim, nDim > S;         // kirchhoff stress
    double                               rho = 1.0; // density
  };

  template < int nDim >
  struct AlgorithmicModuli {
    Fastor::Tensor< double, nDim, nDim, nDim, nDim > dS_dF; // kirchhoff stress
  };

  template < int nDim >
  struct DeformationIncrement {
    Fastor::Tensor< double, nDim, nDim > F_n;
    Fastor::Tensor< double, nDim, nDim > F_np;
  };

  struct TimeIncrement {
    const double* timeOld;
    const double  dT;
  };

  using MarmotMaterial::MarmotMaterial;

  virtual void computeStress( ConstitutiveResponse< 3 >& response,
                              AlgorithmicModuli< 3 >&    tangents,
                              const DeformationIncrement< 3 >&,
                              const TimeIncrement&,
                              double& pNewDT ) = 0;

  /**
   * Compute Stress, but account for eigen deformations (e.g, geostatic stress states). Modifies algorithmic tangent
   * accordingly.*/
  virtual void computeStress( ConstitutiveResponse< 3 >& response,
                              AlgorithmicModuli< 3 >&    tangents,
                              const DeformationIncrement< 3 >&,
                              const TimeIncrement&,
                              double&                                     pNewDT,
                              const std::tuple< double, double, double >& eigenDeformation );

  /*virtual void computePlaneStrain( Fastor::Tensor< double, 3, 3 >&       kirchhoffStress,*/
  /*                                 Fastor::Tensor< double, 3, 3, 3, 3 >& dS_dF,*/
  /*                                 const DeformationIncrement< 3 >&,*/
  /*                                 const TimeIncrement&,*/
  /*                                 double& pNewDT );*/

  /***/
  /* * Compute stress under plane strain conditions, but account for eigen deformations (e.g, geostatic stress
   * states).*/
  /* * Modifies algorithmic tangent accordingly.*/
  /*virtual void computePlaneStrain( Fastor::Tensor< double, 3, 3 >&       kirchhoffStress,*/
  /*                                 Fastor::Tensor< double, 3, 3, 3, 3 >& dS_dF,*/
  /*                                 const DeformationIncrement< 3 >&,*/
  /*                                 const TimeIncrement&,*/
  /*                                 double&                                     pNewDT,*/
  /*                                 const std::tuple< double, double, double >& eigenDeformation );*/

  /*virtual void computePlaneStress( Fastor::Tensor< double, 3, 3 >&       kirchhoffStress,*/
  /*                                 Fastor::Tensor< double, 3, 3, 3, 3 >& dS_dF,*/
  /*                                 const DeformationIncrement< 2 >&,*/
  /*                                 const TimeIncrement&,*/
  /*                                 double& pNewDT );*/

  /** Find for given eigen stress the appropriate eigen deformation.
   * */
  std::tuple< double, double, double > findEigenDeformationForEigenStress(
    const std::tuple< double, double, double >& initialGuess,
    const std::tuple< double, double, double >& eigenStress );
};
