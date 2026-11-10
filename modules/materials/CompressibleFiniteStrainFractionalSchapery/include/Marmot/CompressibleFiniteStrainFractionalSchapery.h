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
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include <Eigen/Dense>
#include <Fastor/Fastor.h>
#include <array>
#include <map>
#include <string>
#include <tuple>

namespace Marmot::Materials {

  /**
   * @class Marmot::Materials::CompressibleFiniteStrainFractionalSchapery
   * @brief Finite-strain compressible Fractional Schapery viscoelastic material model.
   *
   * @details Combines a selectable hyperelastic base potential (Neo-Hooke, Yeoh, Mooney-Rivlin,
   * Pence-Gou, Arruda-Boyce) with a fractional-order power-law relaxation kernel and
   * Schapery non-linear strain-dependent modifier functions (g1, g2, a_S).
   *
   * @par Material parameters
   * - @b hyperelasticBase - selector: 0=NeoHooke, 1=Yeoh, 2=MooneyRivlin, 3=PenceGouNeoHooke, 4=ArrudaBoyce
   * - @b onlyShearCreep - flag: 1.0 restricts viscoelastic creep to isochoric part, 0.0 applies to full stress
   * - @b elasticProperties - shape parameters and bulk modulus for the chosen base model
   * - @b alpha - fractional relaxation exponent (0 < alpha < 1)
   * - @b tau0 - reference relaxation time scale in seconds
   * - @b g_inf - long-term equilibrium stiffness ratio (0 < g_inf < 1)
   * - @b a1 - Schapery g1 parameter: g1 = 1 + a1 * (Ibar1 - 3)
   * - @b a2 - Schapery g2 parameter: g2 = 1 + a2 * (Ibar1 - 3)
   * - @b b_S - Schapery a_S shift parameter: a_S = exp(-b_S * (Ibar1 - 3))
   * - @b rho - optional material density (last parameter if present)
   *
   * @par State variables
   * - @c S0_old - previous effective driving PK2 stress (9 entries)
   * - @c creepStateVars - internal hereditary history variables for 5 discrete Maxwell modes (45 entries)
   */
  class CompressibleFiniteStrainFractionalSchapery : public Marmot::Materials::MaterialFiniteStrain {
  public:
    using Marmot::Materials::MaterialFiniteStrain::MaterialFiniteStrain;

    CompressibleFiniteStrainFractionalSchapery( const double* materialProperties,
                                                int           nMaterialProperties,
                                                int           materialLabel );

    void computeStress( ConstitutiveResponse< 3 >& response,
                        AlgorithmicModuli< 3 >&    tangents,
                        const Deformation< 3 >&    deformation,
                        const TimeIncrement&       timeIncrement ) const override;

    double getDensity( const double* stateVars ) const override;

  protected:
    enum HyperelasticBase { NeoHooke = 0, Yeoh = 1, MooneyRivlin = 2, PenceGouNeoHooke = 3, ArrudaBoyce = 4 };

    const std::map< HyperelasticBase, int > nElasticPropertiesMap = {
      { NeoHooke, 2 },
      { Yeoh, 4 },
      { MooneyRivlin, 3 },
      { PenceGouNeoHooke, 2 },
      { ArrudaBoyce, 3 },
    };

    const HyperelasticBase                    hyperelasticBase;
    const double                              onlyShearCreep;
    const Eigen::Map< const Eigen::VectorXd > elasticProperties;

    // Fractional Schapery parameters
    const double alpha; ///< Fractional power-law exponent (0 < alpha < 1)
    const double tau0;  ///< Reference relaxation time scale
    const double g_inf; ///< Long-term equilibrium fraction (0 < g_inf < 1)
    const double a1;    ///< Schapery g1 parameter: g1 = 1 + a1*(Ibar1 - 3)
    const double a2;    ///< Schapery g2 parameter: g2 = 1 + a2*(Ibar1 - 3)
    const double b_S;   ///< Schapery a_S shift parameter: a_S = exp(-b_S*(Ibar1 - 3))

    static constexpr int                   nFractionalModes = 5;
    std::array< double, nFractionalModes > tau;
    std::array< double, nFractionalModes > gamma;

    void initializeFractionalSpectrum();

    void initializeStateLayout()
    {
      stateLayout.add( "S0_old", 9 );
      stateLayout.add( "creepStateVars", nFractionalModes * 9 );
      stateLayout.finalize();
    }

    std::tuple< double,
                Marmot::Math::Tensors::FastorTensors::Tensor33d,
                Marmot::Math::Tensors::FastorTensors::Tensor3333d,
                Marmot::Math::Tensors::FastorTensors::Tensor333333d >
    computeEnergyDensityAndDerivatives( const Marmot::Math::Tensors::FastorTensors::Tensor33d& C ) const;
  };

} // namespace Marmot::Materials
