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
      TransverseIsotropic = 11, /**< Number of materialProperties equals 8. 5 to describe
                                 *   the material behavior and 3 components of a normal
                                 *   vector. The latter corresponds to the x1 - axis of the principal material
                                 *   directions and defines the materials plane of isotropy (directions x2 and x3).
                                 */
      Orthotropic = 15          /**< Number of materialProperties equals 12. 9 to describe material behavior and
                                 * 	3 components of a orthogonal vectors. The latter corresponds to the
                                 *   x1 and x2 of the principal material directions.
                                 */
    } anisotropicType;          /**< #anisotropicType represents the directional dependence of the material behavior.
                                 *  	It is a representative of the enum class #Type. */

    /// \brief Young's modulus in x1 - direction
    /** #E1 represents the Young's modulus effective in the direction of the x1 - axis
     * of the user defined local coordinate system.
     * It is a reference variable to #materialProperties[0]. */

    const double& E1;

    /// \brief Young's modulus in x2 - direction
    /** #E2 represents the Young's modulus effective in the direction of the x2 - axis
     * of a user defined local coordinate system. In case of isotropic behavior it is
     * set to #E1 and does not need to be specified.
     * Otherwise it is a reference variable to #materialProperties[1]. */

    const double& E2;

    /// \brief Young's modulus in x3 - direction
    /** #E3 represents the Young's modulus effective in the direction of the x3 - axis
     * of a user defined local coordinate system. In case of orthotropic behavior it needs
     * to be specified and is a reference variable to #materialProperties[2].
     *
     * In case of transverse isotropic and isotropic behavior it is automatically set
     * to #E2 and does not need to be specified. */

    const double& E3;

    /// \brief Poisson's ratio between axis x1 and x2.
    /** #nu12 represents the Poisson's ratio effective between the x1 and x2 - axis of
     * the user defined local coordinate system. Depending on the anisotropic type it is a
     * reference variable to:
     * 	- Isotropic:		#materialProperties[1]
     *  	- TransverseIsotropic:	#materialProperties[2]
     *  	- Orthotropic:		#materialProperties[3]. */

    const double& nu12;

    /// \brief Poisson's ratio between axis x2 and x3.
    /** #nu23 represents the Poisson's ratio effective between the x2 and x3 - axis of
     * the user defined local coordinate system. Depending on the anisotropic type it is a
     * reference variable to:
     *  	- TransverseIsotropic:	#materialProperties[3]
     *  	- Orthotropic:		#materialProperties[4].
     *
     * In case of isotropic behavior it is set to #nu12 and does not need to be specified. */

    const double& nu23;

    /// \brief Poisson's ratio between axis x1 and x3.
    /** #nu13 represents the Poisson's ratio effective between the x1 and x3 - axis of
     * the user defined local coordinate system. In case of orthotropic behavior
     * it needs to be specified and is a reference variable to #materialProperties[5].
     *
     * In case of transverse isotropic and isotropic behavior it is set to #nu23 and does not need to be
     * specified.*/

    const double& nu13;

    /// \brief Shear modulus between axis x1 and x2.
    /** #G12 represents the Shear modulus effective between the x1 and x2 - axis of the user defined
     * local coordinate system. In case of isotropic behavior, it does not need to be specified and is
     * set to \f$ \frac{E_1}{2 \cdot ( 1 + \nu_{12})}\f$.
     * Otherwise, depending on the anisotropic type it is a reference variable to:
     *  	- TransverseIsotropic:	#materialProperties[4].
     *  	- Orthotropic:		#materialProperties[6].*/

    const double G12;

    /// \brief Shear modulus between axis x2 and x3.
    /** #G23 represents the Shear modulus effective between the x2 and x3 - axis of the user defined
     * local coordinate system. Neither isotropic nor transversal isotropic behavior needs to be specified.
     * In case of isotropic behavior it is set to #G12 and in case of transverse isotropic behavior
     * it is set to \f$ \frac{E_2}{2 \cdot ( 1 + \nu_{23})}\f$
     * In case of orthotropic behavior it is a reference variable to #materialProperties[7].*/

    const double G23;

    /// \brief Shear modulus between axis x1 and x3.
    /** #G13 represents the Shear modulus effective between the x1 and x3 - axis of the user defined
     * local coordinate system. Neither isotropic nor transversal isotropic behavior needs to be specified and it
     * is set in both cases to #G12.
     * In case of orthotropic behavior it is a reference variable to #materialProperties[8].*/

    const double& G13;

    /// \brief Material stiffness tensor.
    /** #globalStiffnessTensor represents the materials stiffness tensor in voigt notation
     * in the global coordinate system.
     * It is calculated by the functions implemented in \ref MarmotElasticity.h .*/

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
