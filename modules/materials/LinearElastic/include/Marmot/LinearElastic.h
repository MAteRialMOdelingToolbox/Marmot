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
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotTypedefs.h"
#include <iostream>
#include <string>
#include <vector>

namespace Marmot::Materials {
  /**
   * \brief Implementation of a linear elastic material
   * for 3D stress states.
   *
   * For further information see \ref linearelastic.
   */
  class LinearElastic : public MarmotMaterialHypoElastic {
  public:
    using MarmotMaterialHypoElastic::MarmotMaterialHypoElastic;

    LinearElastic( const double* materialProperties, int nMaterialProperties, int materialNumber );

  protected:
    /// \brief Type of isotropic and anisotropic behavior.
    /** #Type is an enum class which involves the following case:*/

    enum class Type {
      Isotropic           = 2,  /**< Number of materialProperties equals 2.*/
      TransverseIsotropic = 11, /**< Number of materialProperties equals 11. 
                                 *   5 parameters describe the material behavior and 
                                 *   6 parameters correspond to two normal vectors that span the isotropic plane.
                                 *   The two vectors lie in the isotropic plane (directions \f$x_2\f$ and \f$x_3\f$).
                                 *   The principal axis \f$x_1\f$ is the unit normal to that plane (e.g. obtained as the normalized cross product of the two vectors).
                                 */
      Orthotropic = 15          /**< Number of materialProperties equals 15. 
                                 * 	 9 parameters describe the material behaviour and
                                 *   6 parameters are the components of two vectors that define the principal material directions.
                                 *   The two vectors correspond to the directions \f$x_1\f$ and \f$x_2\f$.
                                 *   The third principal axis \f$x_3\f$ is the unit normal to the plane spanned by those vectors (e.g. the normalized cross product).
                                 */
    } anisotropicType;          /**< #anisotropicType represents the directional dependence of the material behavior.
                                 *  	It is a representative of the enum class #Type. */

    /// \brief Young's modulus in \f$x_1\f$-direction
    /** \f$E_1\f$ represents the Young's modulus effective in the direction of the \f$x_1\f$-axis of the user-defined local coordinate system.
     * It is a reference variable to #materialProperties[0]. 
     */
    const double& E1;

    /// \brief Young's modulus in \f$x_2\f$-direction
    /** \f$E_2\f$ represents the Young's modulus effective in the direction of the \f$x_2\f$-axis of the user-defined local coordinate system.
     * In case of isotropic behavior it is set to \f$E_1\f$ and does not need to be specified.
     * Otherwise it is a reference variable to #materialProperties[1]. 
     */
    const double& E2;

    /// \brief Young's modulus in \f$x_3\f$-direction
    /** \f$E_3\f$ represents the Young's modulus effective in the direction of the \f$x_3\f$-axis of the user-defined local coordinate system.
     * - For orthotropic behavior it must be specified and is a reference variable
     *   to #materialProperties[2].
     * - For transverse-isotropic or isotropic behavior it is automatically set
     *   to \f$E_2\f$ and does not need to be specified. 
     */
    const double& E3;

    /// \brief Poisson's ratio between axes \f$x_1\f$ and \f$x_2\f$
    /** \details \f$\nu_{12}\f$ represents the Poisson's ratio describing the strain in direction \f$x_1\f$ caused by a load in direction \f$x_2\f$ of the user-defined local coordinate system.
     * Depending on the anisotropic type it is a reference variable to:
     *   - Isotropic:           #materialProperties[1]
     *   - TransverseIsotropic: #materialProperties[2]
     *   - Orthotropic:         #materialProperties[3] 
     */
    const double& nu12;

    /// \brief Poisson's ratio between axes \f$x_2\f$ and \f$x_3\f$
    /** \details \f$\nu_{23}\f$ represents the Poisson's ratio describing the strain in direction \f$x_2\f$ caused by a load in direction \f$x_3\f$ of the user-defined local coordinate system.
     * Depending on the anisotropic type it is a reference variable to:
     *   - TransverseIsotropic: #materialProperties[3]
     *   - Orthotropic:         #materialProperties[4]
     * For isotropic behaviour it is set to #nu12 and does not need to be specified.
     */
     const double& nu23;

    /// \brief Poisson's ratio between axes \f$x_1\f$ and \f$x_3\f$
    /** \details \f$\nu_{13}\f$ represents the Poisson's ratio describing the strain in direction \f$x_1\f$ caused by a load in direction \f$x_3\f$ of the user-defined local coordinate system.
     * - For orthotropic behaviour it must be specified (#materialProperties[5]).
     * - For transverse-isotropic and isotropic behaviour it is set to #nu23 and does not need to be specified.
     */
     const double& nu13;

    /// \brief Shear modulus between axes \f$x_1\f$ and \f$x_2\f$
    /** \details \f$G_{12}\f$ represents the shear modulus effective between the \f$x_1\f$- and \f$x_2\f$-axes of the user-defined local coordinate system.
     * In case of isotropic behavior it does not need to be specified and is set to \f$\tfrac{E_1}{2 \cdot (1 + \nu_{12})}\f$.
     * Otherwise, depending on the anisotropic type it is a reference variable to:
     *   - TransverseIsotropic: #materialProperties[4]  
     *   - Orthotropic:         #materialProperties[6]
     */
    const double G12;

    /// \brief Shear modulus between axes \f$x_2\f$ and \f$x_3\f$
    /** \details \f$G_{23}\f$ represents the shear modulus effective between the \f$x_2\f$- and \f$x_3\f$-axes of the user-defined local coordinate system.
     * For isotropic behavior it is set to #G12.
     * For transverse-isotropic behavior it is set to \f$\tfrac{E_2}{2 \cdot (1 + \nu_{23})}\f$.
     * For orthotropic behavior it is a reference variable to #materialProperties[7].
     */
    const double G23;

    /// \brief Shear modulus between axes \f$x_1\f$ and \f$x_3\f$
    /** \details \f$G_{13}\f$ represents the shear modulus effective between the \f$x_1\f$- and \f$x_3\f$-axes of the user-defined local coordinate system.
     * - For isotropic and transverse-isotropic behavior it is set to #G12.
     * - For orthotropic behavior it is a reference variable to #materialProperties[8].
     */
     const double& G13;

    /// \brief Material stiffness tensor.
    /** #globalStiffnessTensor represents the materials stiffness tensor in voigt notation
     * in the global coordinate system.
     * It is calculated by the functions implemented in *MarmotElasticity.h*.
     */

    Matrix6d globalStiffnessTensor;

    void computeStress( double* stress,
                        double* dStressDDStrain,

                        const double* dStrain,
                        const double* timeOld,
                        const double  dT,
                        double&       pNewDT );

    StateView getStateView( const std::string& result ) { return { nullptr, 0 }; };

    int getNumberOfRequiredStateVars() { return 0; }

    double getDensity();
  };
} // namespace Marmot::Materials
