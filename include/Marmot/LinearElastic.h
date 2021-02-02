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
#include <iostream>
#include <string>
#include <vector>
#include "Marmot/MarmotTypedefs.h"

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
        /// \brief Young's modulus
        /**< #E represents Young's modulus for isotropic linear elasticity.
         * It is a reference variable to #materialProperties[0]. */
        /// \brief Poisson's ratio
        /**< #nu represents Poisson's ratio for isotropic linear elasticity.
         * It is a reference variable to #materialProperties[1]. */
	//std::vector<double> materialProps;
 	
	enum class Type {Isotropic = 2 , TransverseIsotropic = 8 , Orthotropic = 12} anisotropicType; 

     	const double& E1;
     	const double& E2;
     	const double& E3;
     	const double& nu12;
     	const double& nu23;
     	const double& nu13;
     	const double G12;
     	const double G23;
     	const double& G13;

	Matrix6d C;

        void computeStress( double* stress,
                            double* dStressDDStrain,

                            const double* dStrain,
                            const double* timeOld,
                            const double  dT,
                            double&       pNewDT );

        PermanentResultLocation getPermanentResultPointer( const std::string& result ) { return { nullptr, 0 }; };

        int getNumberOfRequiredStateVars() { return 0; }
    };
} // namespace Marmot::Materials
