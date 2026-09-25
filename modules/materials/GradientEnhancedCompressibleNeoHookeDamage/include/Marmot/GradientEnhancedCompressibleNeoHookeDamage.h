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
#include "Marmot/MarmotMaterialGradientEnhancedFiniteStrain.h"

namespace Marmot::Materials {

  /**
   * @class GradientEnhancedCompressibleNeoHookeDamage
   * @brief Compressible neo-Hookean solid with isotropic, nonlocal-field driven damage.
   *
   * A minimal, fully consistent reference model for the MarmotMaterialGradientEnhancedFiniteStrain interface,
   * in the spirit of the implicit-gradient damage model of Peerlings et al. (1996), carried over to finite
   * strains:
   * @f[
   *   \boldsymbol{\tau} = (1 - D(\kappa))\,\boldsymbol{\tau}_0(\boldsymbol{F}), \qquad
   *   \kappa = \max_{t}\bigl(\kappa_0,\,\bar{N}\bigr), \qquad
   *   L = \sqrt{2\,\psi_0(\boldsymbol{F}) / E},
   * @f]
   * where @f$ \psi_0 @f$ and @f$ \boldsymbol{\tau}_0 @f$ are the energy density and Kirchhoff stress of the
   * compressible neo-Hookean potential (PenceGouPotentialB, as in CompressibleNeoHooke), @f$ E = 9KG/(3K+G) @f$,
   * so that @f$ L @f$ is an energy-equivalent strain which reduces to the axial strain in small-strain uniaxial
   * stress. The nonlocal field @f$ \bar{N} @f$ solves @f$ \bar{N} - l^2\nabla^2\bar{N} = L @f$, and
   * @f[
   *   D(\kappa) = 1 - \frac{\kappa_0}{\kappa}\exp\!\left(-\frac{\kappa - \kappa_0}{\kappa_f - \kappa_0}\right)
   *   \quad \text{for} \quad \kappa > \kappa_0, \qquad D = 0 \quad \text{otherwise}.
   * @f]
   *
   * Material properties:
   * | idx | symbol           | meaning                                   |
   * |-----|------------------|-------------------------------------------|
   * | 0   | @f$ K @f$        | bulk modulus                              |
   * | 1   | @f$ G @f$        | shear modulus                             |
   * | 2   | @f$ \kappa_0 @f$ | damage threshold (equivalent strain)      |
   * | 3   | @f$ \kappa_f @f$ | softening parameter, @f$ \kappa_f > \kappa_0 @f$ |
   * | 4   | @f$ l @f$        | nonlocal radius, @f$ c = l^2 @f$          |
   * | 5   | @f$ \rho @f$     | density in the reference configuration    |
   *
   * State variables: @c kappa, the history maximum of the nonlocal field.
   *
   * The dissipation is cumulative: the incoming ConstitutiveResponse::dissipation is incremented by the energy
   * released by the damage increment, @f$ \psi_0\,\Delta D @f$.
   */
  class GradientEnhancedCompressibleNeoHookeDamage : public MarmotMaterialGradientEnhancedFiniteStrain {
  public:
    GradientEnhancedCompressibleNeoHookeDamage( const double* materialProperties,
                                                int           nMaterialProperties,
                                                int           materialNumber );

    using MarmotMaterialGradientEnhancedFiniteStrain::computeStress;

    void computeStress( ConstitutiveResponse< 3 >& response,
                        AlgorithmicModuli< 3 >&    tangents,
                        const Deformation< 3 >&    deformation,
                        const TimeIncrement&       timeIncrement ) const override;

    double getDensity( const double* stateVars ) const override;

    /**
     * @brief Damage and its derivative for a given history variable.
     * @param[in] kappa History maximum of the nonlocal field.
     * @return @f$ \{D,\ \partial D/\partial\kappa\} @f$.
     */
    std::pair< double, double > damage( double kappa ) const;

  private:
    const double& K;
    const double& G;
    const double& kappa0;
    const double& kappaF;
    const double& nonLocalRadius;
  };

} // namespace Marmot::Materials
