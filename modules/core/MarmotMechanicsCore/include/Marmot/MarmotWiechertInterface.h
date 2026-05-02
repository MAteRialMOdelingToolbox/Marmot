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

  namespace WiechertInterface {

    typedef Eigen::VectorXd          Properties;
    typedef Eigen::Map< Properties > mapProperties;

    typedef Eigen::Matrix< double, 3, Eigen::Dynamic > StateVarMatrix_force_uu;
    typedef Eigen::Matrix< double, 3, Eigen::Dynamic > StateVarMatrix_force_us;
    typedef Eigen::Matrix< double, 9, Eigen::Dynamic > StateVarMatrix_surface_stress_Z;
    typedef Eigen::Matrix< double, 9, Eigen::Dynamic > StateVarMatrix_surface_stress_Y;
    typedef Eigen::Matrix< double, 9, Eigen::Dynamic > StateVarMatrix_surface_stress_us;

    typedef Eigen::Map< StateVarMatrix_force_uu >          mapStateVarMatrix_force_uu;
    typedef Eigen::Map< StateVarMatrix_force_us >          mapStateVarMatrix_force_us;
    typedef Eigen::Map< StateVarMatrix_surface_stress_Z >  mapStateVarMatrix_surface_stress_Z;
    typedef Eigen::Map< StateVarMatrix_surface_stress_Y >  mapStateVarMatrix_surface_stress_Y;
    typedef Eigen::Map< StateVarMatrix_surface_stress_us > mapStateVarMatrix_surface_stress_us;

    // template < int k >
    // Properties computeElasticModuli_Ru( std::function< autodiff::Real< k, double >( autodiff::Real< k, double > ) >
    // phi,
    //                                  Properties retardationTimes_Ru,
    //                                  bool       gaussQuadrature = false )
    //{
    //   Properties elasticModuli_Ru( retardationTimes_Ru.size() );
    //   double     spacing = retardationTimes_Ru( 1 ) / retardationTimes_Ru( 0 );
    //   for ( int i = 0; i < retardationTimes_Ru.size(); i++ ) {
    //     double tau = retardationTimes_Ru( i );
    //     if ( !gaussQuadrature ) {
    //
    //       elasticModuli_Ru( i ) = 1. / ( log( spacing ) * KelvinChain::evaluatePostWidderFormula< k >( phi, tau ) );
    //
    //    }
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
    //                                  Properties retardationTimes_Rs,
    //                                  bool       gaussQuadrature = false )
    //{
    //   Properties elasticModuli_Rs( retardationTimes_Rs.size() );
    //   double     spacing = retardationTimes_Rs( 1 ) / retardationTimes_Rs( 0 );
    //
    //      for ( int i = 0; i < retardationTimes_Rs.size(); i++ ) {
    //        double tau = retardationTimes_Rs( i );
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

    void updateStateVarMatrix_force_uu( const double                          dT,
                                        Properties                            elasticModuli,
                                        Properties                            relaxationTimes,
                                        Eigen::Ref< StateVarMatrix_force_uu > stateVars_force_uu,
                                        const Marmot::Vector3d&               dforce_uu,
                                        const Marmot::Matrix3d&               unitH_inv_ij );

    void updateStateVarMatrix_force_us( const double                          dT,
                                        Properties                            elasticModuli,
                                        Properties                            relaxationTimes,
                                        Eigen::Ref< StateVarMatrix_force_us > stateVars_force_us,
                                        const Marmot::Vector9d&               dforce_us,
                                        const Eigen::Matrix< double, 3, 9 >&  unitH_inv_nF_ijk );

    void updateStateVarMatrix_surface_stress_Z(
      const double                                  dT,
      Properties                                    elasticModuli,
      Properties                                    relaxationTimes,
      Eigen::Ref< StateVarMatrix_surface_stress_Z > stateVars_surface_stress_Z,
      const Marmot::Vector9d&                       dsurfaceStress_Z,
      const Marmot::Matrix9d&                       unitZ_ijkl );

    void updateStateVarMatrix_surface_stress_Y(
      const double                                  dT,
      Properties                                    elasticModuli,
      Properties                                    relaxationTimes,
      Eigen::Ref< StateVarMatrix_surface_stress_Y > stateVars_surface_stress_Y,
      const Marmot::Vector9d&                       dsurfaceStress_Y,
      const Marmot::Matrix9d&                       unitYn_H_inv_Fn_ijkl );

    void updateStateVarMatrix_surface_stress_us(
      const double                                   dT,
      Properties                                     elasticModuli,
      Properties                                     relaxationTimes,
      Eigen::Ref< StateVarMatrix_surface_stress_us > stateVars_surface_stress_us,
      const Marmot::Vector3d&                        djumpU,
      const Eigen::Matrix< double, 3, 9 >&           unitH_inv_nF_ijk );

    void evaluateWiechert( const double                     dT,
                           Properties                       elasticModuli,
                           Properties                       relaxationTimes,
                           StateVarMatrix_force_uu          stateVars_force_uu,
                           StateVarMatrix_force_us          stateVars_force_us,
                           StateVarMatrix_surface_stress_Z  stateVars_surface_stress_Z,
                           StateVarMatrix_surface_stress_Y  stateVars_surface_stress_Y,
                           StateVarMatrix_surface_stress_us stateVars_surface_stress_us,
                           double&                          uniaxialStiffness,
                           Marmot::Vector3d&                dforce_uu,
                           Marmot::Vector3d&                dforce_us,
                           Marmot::Vector9d&                dsurfaceStress_Z,
                           Marmot::Vector9d&                dsurfaceStress_Y,
                           Marmot::Vector9d&                dsurfaceStress_us,
                           const double                     factor );

    void computeLambdaAndBeta( double dT, double tau, double& lambda, double& beta );

  } // namespace WiechertInterface
} // namespace Marmot::Materials
