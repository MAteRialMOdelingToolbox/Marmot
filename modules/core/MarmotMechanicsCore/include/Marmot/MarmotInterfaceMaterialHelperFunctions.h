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
#include "Marmot/MarmotFastorTensorBasics.h"

#include <Eigen/Core>
#include <tuple>

namespace Marmot::Materials {

  namespace InterfaceMaterialHelperFunctions {

    using Marmot::FastorStandardTensors::Tensor3333d;
    using Marmot::FastorStandardTensors::Tensor333d;
    using Marmot::FastorStandardTensors::Tensor33d;
    using Marmot::FastorStandardTensors::Tensor3d;

    /**
     * @brief Computes the geometry-system coupling tensors for an interface.
     *
     * Forms the contracted matrix \f$Q = L:N\f$ and its inverse \f$G\f$,
     * then computes the coupling tensor \f$A = G \otimes N\f$ and the
     * condensed constitutive tensor \f$B = L - L:A:L\f$.
     *
     * @param[in] N Normal projector \f$\boldsymbol{n} \otimes \boldsymbol{n}\f$.
     * @param[in] L Fourth-order constitutive tensor.
     *
     * @return Tuple containing, in order:
     *         - \f$B\f$: condensed fourth-order constitutive tensor,
     *         - \f$L\f$: unchanged input constitutive tensor,
     *         - \f$A\f$: fourth-order geometry-system coupling tensor,
     *         - \f$G = Q^{-1}\f$: inverse contracted constitutive matrix.
     */
    std::tuple< Tensor3333d, const Tensor3333d, Tensor3333d, Tensor33d > interfaceGeometrySystemCouplings(
      const Tensor33d&   N,
      const Tensor3333d& L );

    /**
     * @brief Computes the four constitutive operators used by an interface material.
     *
     * Condenses the bulk constitutive tensor with respect to the interface
     * normal direction and derives the operators coupling displacement jumps,
     * surface strains, interface forces, and surface stresses.
     *
     * @param[in] normal Unit normal vector of the interface.
     * @param[in] N Normal projector \f$\boldsymbol{n} \otimes \boldsymbol{n}\f$.
     * @param[in] C_nu_aibj Fourth-order constitutive tensor.
     *
     * @return Tuple containing, in order:
     *         - \f$Z\f$: condensed fourth-order surface-strain operator,
     *         - \f$H^{-1}\f$: second-order displacement-jump operator,
     *         - \f$H^{-1}nF\f$: third-order jump/surface-strain coupling operator,
     *         - \f$nYH^{-1}Fn\f$: fourth-order coupling correction.
     */
    std::tuple< Tensor3333d, Tensor33d, Tensor333d, Tensor3333d > calculateMaterialMatrices(
      const Tensor3d&    normal,
      const Tensor33d&   N,
      const Tensor3333d& C_nu_aibj );

    /**
     * @brief Calculates the FY tensor components for interface material formulation.
     *
     * This function computes six fourth and second-order tensor quantities that represent
     * the FY components used in the interface element material formulation. These tensors
     * are derived from the transformation matrices and elastic stiffness tensor.
     *
     * @param[in] N           3x3 direction cosine matrix for the interface normal direction
     * @param[in] C_nu_aibj   4th-order elastic stiffness tensor (in Voigt-like notation)
     *
     * @return A tuple containing six tensors:
     *         - Tensor3333d: 1st FY component (4th-order tensor)
     *         - Tensor3333d: 2nd FY component (4th-order tensor)
     *         - Tensor3333d: 3rd FY component (4th-order tensor)
     *         - Tensor3333d: 4th FY component (4th-order tensor)
     *         - Tensor33d:   5th FY component (2nd-order tensor)
     *         - Tensor3333d: 6th FY component (4th-order tensor)
     */
    std::tuple< Tensor3333d, Tensor3333d, Tensor3333d, Tensor3333d, Tensor33d, Tensor3333d > calculateFY(
      const Tensor33d&   N,
      const Tensor3333d& C_nu_aibj );

    /**
     * @brief Calculates interface material parameters from a normal vector and Poisson's ratio.
     *
     * This overload computes interface material matrices assuming an isotropic elastic material
     * defined by a normal direction and Poisson's ratio. The elastic properties are derived
     * from these parameters.
     *
     * @param[in] normal   3-component unit normal vector to the interface
     * @param[in] nu_0     Poisson's ratio of the interface material
     *
     * @return A tuple containing four tensors:
     *         - Tensor3333d: 1st material matrix parameter (4th-order tensor)
     *         - Tensor33d:   2nd material matrix parameter (2nd-order tensor)
     *         - Tensor333d:  3rd material matrix parameter (3rd-order tensor)
     *         - Tensor3333d: 4th material matrix parameter (4th-order tensor)
     */
    std::tuple< Tensor3333d, Tensor33d, Tensor333d, Tensor3333d > calculateInterfaceMaterialParameters(
      const Tensor3d& normal,
      const double&   nu_0 );

    /**
     * @brief Calculates interface material parameters from a normal vector and elastic stiffness matrix.
     *
     * This overload computes interface material matrices given a normal direction and a 6x6 Voigt
     * representation of the elastic stiffness tensor. This allows specification of arbitrary
     * (including anisotropic) elastic properties for the interface material.
     *
     * @param[in] normal      3-component unit normal vector to the interface
     * @param[in] C_ep_voigt  6x6 elastic stiffness matrix in Voigt notation
     *
     * @return A tuple containing four tensors:
     *         - Tensor3333d: 1st material matrix parameter (4th-order tensor)
     *         - Tensor33d:   2nd material matrix parameter (2nd-order tensor)
     *         - Tensor333d:  3rd material matrix parameter (3rd-order tensor)
     *         - Tensor3333d: 4th material matrix parameter (4th-order tensor)
     */
    std::tuple< Tensor3333d, Tensor33d, Tensor333d, Tensor3333d > calculateInterfaceMaterialParameters(
      const Tensor3d&                      normal,
      const Eigen::Matrix< double, 6, 6 >& C_ep_voigt );
  } // namespace InterfaceMaterialHelperFunctions
} // namespace Marmot::Materials
