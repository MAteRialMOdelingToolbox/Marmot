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
#include "Marmot/MarmotTransportMaterial.h"
#include "Marmot/MarmotTypedefs.h"

namespace Marmot::Materials {

  /**
   * @class Marmot::Materials::HeatConduction
   * @brief Fourier heat conduction with an anisotropic, temperature dependent conductivity and a temperature
   * dependent specific heat capacity.
   *
   * One field, the temperature \f$T\f$. With \f$\theta = T - T_\mathrm{ref}\f$:
   * \f[
   *   \mathbf{q} = -\mathbf{k}(T)\,\nabla T,\qquad \mathbf{k}(T) = \mathbf{k}_0\,(1 + \beta_k\,\theta),\qquad
   *   c(T) = c_0\,(1 + \beta_c\,\theta),
   * \f]
   * and the content is the enthalpy per unit volume
   * \f[
   *   H(T) = \rho \int_{T_\mathrm{ref}}^{T} c(\tilde T)\,\mathrm{d}\tilde T
   *        = \rho\,c_0 \left(\theta + \tfrac{1}{2}\beta_c\,\theta^2\right),
   * \f]
   * whose increment over the time increment is returned (see Marmot::Materials::TransportMaterial).
   *
   * Material properties, in this order (3, 6 or 11 of them):
   *
   * | # | symbol | meaning | unit (SI) |
   * |---|--------|---------|-----------|
   * | 0 | \f$\rho\f$ | density | kg/m³ |
   * | 1 | \f$c_0\f$ | specific heat capacity at \f$T_\mathrm{ref}\f$ | J/(kg K) |
   * | 2 | \f$k_{11}\f$ (or \f$k\f$) | conductivity at \f$T_\mathrm{ref}\f$ | W/(m K) |
   * | 3 | \f$k_{22}\f$ | | W/(m K) |
   * | 4 | \f$k_{33}\f$ | | W/(m K) |
   * | 5 | \f$k_{12}\f$ | | W/(m K) |
   * | 6 | \f$k_{13}\f$ | | W/(m K) |
   * | 7 | \f$k_{23}\f$ | | W/(m K) |
   * | 8 | \f$\beta_k\f$ | linear temperature coefficient of the conductivity | 1/K |
   * | 9 | \f$\beta_c\f$ | linear temperature coefficient of the heat capacity | 1/K |
   * | 10 | \f$T_\mathrm{ref}\f$ | reference temperature | K or °C |
   *
   * The short forms are: 3 properties \f$[\rho, c_0, k]\f$ (isotropic, constant), and 6 properties
   * \f$[\rho, c_0, k, \beta_k, \beta_c, T_\mathrm{ref}]\f$ (isotropic, temperature dependent).
   * A zero heat capacity \f$\rho c_0 = 0\f$ gives steady state heat conduction.
   *
   * @throws std::invalid_argument for another number of properties, a negative \f$\rho c_0\f$ or a conductivity
   * which is not positive definite.
   */
  class HeatConduction : public TransportMaterial {
  public:
    /**
     * @brief Construct the material.
     * @param materialProperties Pointer to the material properties, see the class description.
     * @param nMaterialProperties Number of material properties: 3, 6 or 11.
     * @param materialNumber Material identifier.
     */
    HeatConduction( const double* materialProperties, int nMaterialProperties, int materialNumber );

    /// @return 1, the temperature
    int getNumberOfFields() const override { return 1; }

    /// @return False for a zero heat capacity \f$\rho c_0\f$: steady state heat conduction
    bool hasStorage( [[maybe_unused]] int field ) const override { return rhoC0 > 0.0; }

    /**
     * @brief Enthalpy increment, heat flux and their algorithmic tangent.
     * @param[inout] response The response.
     * @param[out] tangent The algorithmic tangent.
     * @param[in] fieldState Temperatures and temperature gradient.
     * @param[in] timeIncrement Time information (not used).
     * @throws std::runtime_error if the conductivity factor \f$1 + \beta_k\theta\f$ is not positive.
     */
    void computeTransport( Response&            response,
                           AlgorithmicTangent&  tangent,
                           const FieldState&    fieldState,
                           const TimeIncrement& timeIncrement ) const override;

    /**
     * @brief The enthalpy per unit volume \f$H(T)\f$ relative to \f$T_\mathrm{ref}\f$.
     * @param T The temperature.
     * @return \f$H(T)\f$.
     */
    double enthalpy( double T ) const;

  private:
    double                 rhoC0;      ///< \f$\rho c_0\f$
    Marmot::Math::Matrix3d k0;         ///< conductivity at the reference temperature
    double                 betaK;      ///< \f$\beta_k\f$
    double                 betaC;      ///< \f$\beta_c\f$
    double                 TReference; ///< \f$T_\mathrm{ref}\f$
  };

} // namespace Marmot::Materials
