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
 *
 * Alexandros Stathas alexandros.stathas@boku.ac.at
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
#include "Marmot/MarmotInterfaceMaterialHypoElastic.h"

namespace Marmot::Materials {
  /**
   * \brief BMGu variant of the linear elastic interface material for 3D stress states.
   *
   * This model is intentionally different from \c LINEARELASTICINTERFACE. It requires the material of the outer blocks
   * toparticipate in the mechanical behavior of the interface. Both models use the same kinematic helper tensors (see
   * \ref linearelasticinterface), but they differ in parameterization and in the
   * constitutive coefficients assigned to the interface operators.
   *
   * Differences to \c LINEARELASTICINTERFACE:
   * - \b Parameterization:
   *   - \c LINEARELASTICINTERFACE uses 3 parameters: $E_0, \nu_0, h$.
   *   - \c LINEARELASTICINTERFACEBMGU uses 7 parameters with the following layout:
   *     - \c materialProperties[0] = $E_M$
   *     - \c materialProperties[1] = $\nu_M$
   *     - \c materialProperties[2] = $E_I$
   *     - \c materialProperties[3] = $\nu_I$
   *     - \c materialProperties[4] = $E_0$
   *     - \c materialProperties[5] = $\nu_0$
   *     - \c materialProperties[6] = $h$
   *     - \c materialProperties[1] and \c materialProperties[3] are currently unused. We assume the materials have the
   * same Poisson's ratio, this simplifies the material description for the sake of the argument, but we keep these
   * parameters for potential future use.
   * - \b Constitutive scaling:
   *   - Define
   *     $H_{\mathrm{bar}} = \frac{2}{E_0} - \frac{1}{E_M} - \frac{1}{E_I}$ and
   *     $Z_{\mathrm{bar}} = E_M + E_I - 2 E_0$.
   *   - \c LINEARELASTICINTERFACEBMGU sets
   *     $Z_{ijkl} = -\frac{h}{2} Z_{\mathrm{bar}}\, \hat{Z}_{ijkl}$,
   *     Due to the Poisson ration being the same between the three materials the following simplifications apply to the
   * other two operators: $Y_{n}H^{-1}F_{n} = 0$, $H^{-1}_{ij} = \frac{2}{h\,H_{\mathrm{bar}}}\,\hat{H}^{-1}_{ij}$, and
   * $H^{-1}nF = 0$.
   *   - \c LINEARELASTICINTERFACE instead uses the direct $E_0$-scaled terms for
   *     all four operators, and these operators do not vanish anymore.
   *
   * For further information see \ref linearelasticinterface.
   */
  class LinearElasticInterfaceBMGu : public MarmotInterfaceMaterialHypoElastic {
  public:
    using MarmotInterfaceMaterialHypoElastic::MarmotInterfaceMaterialHypoElastic;

    LinearElasticInterfaceBMGu( const double* materialProperties, int nMaterialProperties, int materialNumber );

    void computeStress( State&               state,
                        Tangents&            tangents,
                        const Deformation&   deformation,
                        const TimeIncrement& timeIncrement ) override;

    int getNumberOfRequiredStateVars() const override { return 0; }

  private:
    const double& E_M;
    const double& nu_M;
    const double& E_I;
    const double& nu_I;
    const double& E_0;
    const double& nu_0;
    const double& h;
    const double  Hbar = ( 2. / E_0 ) - ( 1. / E_M ) - ( 1. / E_I );
    const double  Zbar = ( E_M ) + ( E_I ) - ( 2. * E_0 );
  };
} // namespace Marmot::Materials
