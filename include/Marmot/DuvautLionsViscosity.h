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
#include "Marmot/MarmotTypedefs.h"

namespace Marmot {
  namespace ContinuumMechanics::CommonConstitutiveModels {
    /**
     * \brief Implementation of Duvaut-Lions viscosity for a material with \ref nMatTangentSize internal degrees of
     * freedom
     * @todo: Update member names to more descriptive ones
     */
    template < int nMatTangentSize >
    class DuvautLionsViscosity {
    private:
      /*
       * Viscosity parameter for Duvault-Lions viscosity*/
      const double viscosity;

    public:
      typedef Eigen::Matrix< double, nMatTangentSize, nMatTangentSize > TangentSizedMatrix;

      DuvautLionsViscosity( double viscosity );

      /**
       * Apply viscosity on a scalar internal variable depending on the extremal solutions for t=0 (trialState) and
       * t=\f$\infty\f$, and timestep dt */
      double applyViscosityOnStateVar( double stateVarTrial, double StateVarInf, double dT );
      /**
       * Apply viscosity on voigt sized rank two tensor depending on the extremal solutions for t=0 (trialState) and
       * t=\f$\infty\f$, and timestep dt */
      Marmot::Vector6d applyViscosityOnStress( const Marmot::Vector6d& trialStress,
                                               const Marmot::Vector6d& stressInf,
                                               double                  dT );
      /**
       * Apply viscosity on the inverse (algorithmic) material tangent depending on the extremal solutions for t=0
       * (trialState) and t=\f$\infty\f$, and timestep dt
       * @todo: Check if application to inverse can be replaced by application to non-inverse in general*/
      TangentSizedMatrix applyViscosityOnMatTangent( const TangentSizedMatrix& matTangentInv, double dT );
    };
  } // namespace ContinuumMechanics::CommonConstitutiveModels
} // namespace Marmot

namespace Marmot {
  namespace ContinuumMechanics::CommonConstitutiveModels {
    template < int s >
    DuvautLionsViscosity< s >::DuvautLionsViscosity( double viscosity ) : viscosity( viscosity )
    {
    }

    template < int s >
    double DuvautLionsViscosity< s >::applyViscosityOnStateVar( double stateVarTrial, double StateVarInf, double dT )
    {
      return ( stateVarTrial + ( dT / viscosity ) * StateVarInf ) / ( dT / viscosity + 1 );
    }

    template < int s >
    Marmot::Vector6d DuvautLionsViscosity< s >::applyViscosityOnStress( const Marmot::Vector6d& trialStress,
                                                                        const Marmot::Vector6d& stressInf,
                                                                        double                  dT )
    {
      return ( trialStress + ( dT / viscosity ) * stressInf ) / ( dT / viscosity + 1 );
    }

    template < int s >
    typename DuvautLionsViscosity< s >::TangentSizedMatrix DuvautLionsViscosity< s >::applyViscosityOnMatTangent(
      const TangentSizedMatrix& matTangentInv,
      double                    dT )
    {
      return ( 1 / ( 1 + dT / viscosity ) ) * ( TangentSizedMatrix::Identity() + dT / viscosity * matTangentInv );
    }
  } // namespace ContinuumMechanics::CommonConstitutiveModels
} // namespace Marmot
