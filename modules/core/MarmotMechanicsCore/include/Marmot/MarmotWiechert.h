/* ---------------------------------------------------------------------
 *                                       _
 *  _ __ ___   __ _ _ __ ___ ___   ___ | |_
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

#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotViscoelasticity.h"

#include <cmath>
#include <functional>

namespace Marmot::Materials {

  namespace Wiechert {

    using Properties        = Marmot::ContinuumMechanics::Viscoelasticity::DiscreteSpectrum::Properties;
    using mapProperties     = Marmot::ContinuumMechanics::Viscoelasticity::DiscreteSpectrum::mapProperties;
    using StateVarMatrix    = Marmot::ContinuumMechanics::Viscoelasticity::DiscreteSpectrum::StateVarMatrix;
    using mapStateVarMatrix = Marmot::ContinuumMechanics::Viscoelasticity::DiscreteSpectrum::mapStateVarMatrix;

    template < int k >
    double evaluatePostWidderFormula( std::function< autodiff::Real< k, double >( autodiff::Real< k, double > ) > psi,
                                      double                                                                      tau )
    {
      using Marmot::ContinuumMechanics::Viscoelasticity::DiscreteSpectrum::PostWidderCoefficientSign;
      return Marmot::ContinuumMechanics::Viscoelasticity::DiscreteSpectrum::evaluatePostWidderFormula<
        k >( psi, tau, PostWidderCoefficientSign::Positive );
    }

    template < int k >
    Properties computeElasticModuli( std::function< autodiff::Real< k, double >( autodiff::Real< k, double > ) > psi,
                                     const Properties& relaxationTimes,
                                     double            spacing,
                                     bool              gaussQuadrature = false )
    {
      Properties elasticModuli( relaxationTimes.size() );
      for ( int i = 0; i < relaxationTimes.size(); ++i ) {
        const double tau = relaxationTimes( i );
        if ( !gaussQuadrature ) {
          elasticModuli( i ) = std::log( spacing ) * evaluatePostWidderFormula< k >( psi, tau );
        }
        else {
          elasticModuli(
            i ) = std::log( spacing ) / 2. *
                  ( evaluatePostWidderFormula< k >( psi, tau * std::pow( spacing, -std::sqrt( 3. ) / 6. ) ) +
                    evaluatePostWidderFormula< k >( psi, tau * std::pow( spacing, std::sqrt( 3. ) / 6. ) ) );
        }
      }
      return elasticModuli;
    }

    Properties generateRelaxationTimes( int n, double min, double spacing );

    void updateStateVarMatrix( const double                 dT,
                               const Properties&            elasticModuli,
                               const Properties&            relaxationTimes,
                               Eigen::Ref< StateVarMatrix > stateVars,
                               const Marmot::Vector6d&      dStrain,
                               const Marmot::Matrix6d&      unitD_ijkl );

    void evaluateWiechert( const double                 dT,
                           const Properties&            elasticModuli,
                           const Properties&            relaxationTimes,
                           Eigen::Ref< StateVarMatrix > stateVars,
                           double&                      uniaxialStiffness,
                           Marmot::Vector6d&            dStress,
                           const double                 factor );

  } // namespace Wiechert
} // namespace Marmot::Materials
