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
 * Alexander Dummer alexander.dummer@uibk.ac.at
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
#include "Marmot/MarmotNumericalIntegration.h"
#include "Marmot/MarmotTypedefs.h"
#include "autodiff/forward/real.hpp"
#include <functional>

namespace Marmot::Materials {

  namespace KelvinChain {

    typedef Eigen::VectorXd          Properties;
    typedef Eigen::Map< Properties > mapProperties;

    typedef Eigen::Matrix< double, 6, Eigen::Dynamic > StateVarMatrix;
    typedef Eigen::Map< StateVarMatrix >               mapStateVarMatrix;

    template < int N >
    struct Factorial {
      enum { value = N * Factorial< N - 1 >::value };
    };

    template <>
    struct Factorial< 0 > {
      enum { value = 1 };
    };

    template < int k >
    double evaluatePostWidderFormula( std::function< autodiff::Real< k, double >( autodiff::Real< k, double > ) > phi,
                                      double                                                                      tau )
    {

      autodiff::Real< k, double > tau_( tau * k );

      double val = -pow( -tau * k, k ) / double( Factorial< k - 1 >::value );
      val *= autodiff::derivatives( phi, autodiff::along( 1. ), autodiff::at( tau_ ) )[k];
      return val;
    }

    template < int k >
    double approximateZerothCompliance( std::function< autodiff::Real< k, double >( autodiff::Real< k, double > ) > phi,
                                        double tauMin,
                                        double spacing = 10. )
    {
      NumericalAlgorithms::Integration::scalar_to_scalar_function_type f = [&]( double tau ) {
        double                      val_ = -pow( -k, k ) * pow( tau, k - 1 ) / double( Factorial< k - 1 >::value );
        autodiff::Real< k, double > tau_( tau * k );
        val_ *= autodiff::derivatives( phi, autodiff::along( 1. ), autodiff::at( tau_ ) )[k];
        return val_;
      };

      double
        val = NumericalAlgorithms::Integration::integrateScalarFunction( f,
                                                                         { 1e-14, tauMin / sqrt( spacing ) },
                                                                         100,
                                                                         NumericalAlgorithms::Integration::simpson );
      return val;
    }

    template < int k >
    Properties computeElasticModuli( std::function< autodiff::Real< k, double >( autodiff::Real< k, double > ) > phi,
                                     Properties retardationTimes,
                                     bool       gaussQuadrature = false )
    {
      Properties elasticModuli( retardationTimes.size() );
      double     spacing = retardationTimes( 1 ) / retardationTimes( 0 );

      for ( int i = 0; i < retardationTimes.size(); i++ ) {
        double tau = retardationTimes( i );
        if ( !gaussQuadrature ) {
          elasticModuli( i ) = 1. / ( log( spacing ) * evaluatePostWidderFormula< k >( phi, tau ) );
        }
        else {
          elasticModuli( i ) = 1. /
                               ( log( spacing ) / 2. *
                                 ( evaluatePostWidderFormula< k >( phi, tau * pow( spacing, -sqrt( 3. ) / 6. ) ) +
                                   evaluatePostWidderFormula< k >( phi, tau * pow( spacing, sqrt( 3. ) / 6. ) ) ) );
        }
      }

      return elasticModuli;
    }

    Properties generateRetardationTimes( int n, double min, double spacing );

    void updateStateVarMatrix( const double                 dT,
                               Properties                   elasticModuli,
                               Properties                   retardationTimes,
                               Eigen::Ref< StateVarMatrix > stateVars,
                               const Marmot::Vector6d&      dStress,
                               const Marmot::Matrix6d&      unitComplianceMatrix );

    void evaluateKelvinChain( const double      dT,
                              Properties        elasticModuli,
                              Properties        retardationTimes,
                              StateVarMatrix    stateVars,
                              double&           uniaxialCompliance,
                              Marmot::Vector6d& dStrain,
                              const double      factor );

    void computeLambdaAndBeta( double dT, double tau, double& lambda, double& beta );

  } // namespace KelvinChain
} // namespace Marmot::Materials
