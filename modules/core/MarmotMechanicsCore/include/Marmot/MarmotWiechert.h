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
#include "Marmot/MarmotKelvinChain.h"
#include "Marmot/MarmotNumericalIntegration.h"
#include "Marmot/MarmotTypedefs.h"
#include "autodiff/forward/real.hpp"
#include <fstream>  // **needed for std::ofstream**
#include <functional>
#include <iostream> // for std::cout / std::cerr
//
namespace Marmot::Materials {

  namespace Wiechert {

    typedef Eigen::VectorXd          Properties;
    typedef Eigen::Map< Properties > mapProperties;

    typedef Eigen::Matrix< double, 6, Eigen::Dynamic > StateVarMatrix;

    typedef Eigen::Map< StateVarMatrix > mapStateVarMatrix;

    // template < int k >
    // Properties computeElasticModuli_Ru( std::function< autodiff::Real< k, double >( autodiff::Real< k, double > ) >
    // phi,
    //                                  Properties relaxationTimes_Ru,
    //                                  bool       gaussQuadrature = false )
    //{
    //   Properties elasticModuli_Ru( relaxationTimes_Ru.size() );
    //   double     spacing = relaxationTimes_Ru( 1 ) / relaxationTimes_Ru( 0 );
    //   for ( int i = 0; i < relaxationTimes_Ru.size(); i++ ) {
    //     double tau = relaxationTimes_Ru( i );
    //     if ( !gaussQuadrature ) {
    //
    //       elasticModuli_Ru( i ) = 1. / ( log( spacing ) * KelvinChain::evaluatePostWidderFormula< k >( phi, tau ) );
    //
    //    else {
    //      elasticModuli_Ru( i ) = 1. /
    //                           ( log( spacing ) / 2. *
    //                             ( KelvinChain::evaluatePostWidderFormula< k >( phi, tau * pow( spacing, -sqrt( 3. )
    //                             / 6. ) ) +
    //                               KelvinChain::evaluatePostWidderFormula< k >( phi, tau * pow( spacing, sqrt( 3. )
    //                               / 6. ) ) ) );
    //    }
    //
    //      return elasticModuli_Ru;
    //    }

    // template < int k >
    // Properties computeElasticModuli_Rs( std::function< autodiff::Real< k, double >( autodiff::Real< k, double > ) >
    // phi,
    //                                  Properties relaxationTimes_Rs,
    //                                  bool       gaussQuadrature = false )
    //{
    //   Properties elasticModuli_Rs( relaxationTimes_Rs.size() );
    //   double     spacing = relaxationTimes_Rs( 1 ) / relaxationTimes_Rs( 0 );
    //
    //      for ( int i = 0; i < relaxationTimes_Rs.size(); i++ ) {
    //        double tau = relaxationTimes_Rs( i );
    //        if ( !gaussQuadrature ) {
    //          elasticModuli_Rs( i ) = 1. / ( log( spacing ) * KelvinChain::evaluatePostWidderFormula< k >( phi, tau )
    //          );
    //       }
    //       else {
    //         elasticModuli_Rs( i ) = 1. /
    //                              ( log( spacing ) / 2. *
    //                               ( KelvinChain::evaluatePostWidderFormula< k >( phi, tau * pow( spacing, -sqrt( 3. )
    //                               / 6. ) ) +
    //                                 KelvinChain::evaluatePostWidderFormula< k >( phi, tau * pow( spacing, sqrt( 3. )
    //                                 / 6. ) ) ) );
    //       }
    //     }
    //
    //     return elasticModuli_Rs;
    //  }

    // Properties generateRelaxationTimes( int n, double min, double spacing );

    Properties initializeElasticModuli( int nMaxwell, double n );

    Properties initializeRelaxationTimes( int nMaxwell, double m );

    void updateStateVarMatrix( const double                 dT,
                               Properties                   elasticModuli,
                               Properties                   relaxationTimes,
                               Eigen::Ref< StateVarMatrix > stateVars,
                               const Marmot::Vector6d&      dStress,
                               const Marmot::Matrix6d&      unitD_ijkl );

    void evaluateWiechert( const double      dT,
                           Properties        elasticModuli,
                           Properties        relaxationTimes,
                           StateVarMatrix    stateVars,
                           double&           uniaxialStiffness,
                           Marmot::Vector6d& dStress,
                           const double      factor );

    void computeLambdaAndBeta( double dT, double tau, double& lambda, double& beta );

  } // namespace Wiechert
} // namespace Marmot::Materials
