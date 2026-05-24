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

#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include "Marmot/MarmotNumericalDifferentiation.h"
#include <cmath>
#include <cstring>
#include <memory>
#include <stdexcept>

namespace Marmot::Materials {

  using namespace Fastor;

  /**
   * @class MarmotMaterialFiniteStrainSubstepped
   * @brief A decorator that applies time substepping with analytical tangent accumulation.
   * @tparam BaseMaterialType The concrete material class to wrap.
   * Properties: [nSubsteps, Prop1, Prop2, ...]
   */
  template < typename BaseMaterialType >
  class MarmotMaterialFiniteStrainSubstepped : public MarmotMaterialFiniteStrain {
  protected:
    std::unique_ptr< BaseMaterialType > baseMaterial;
    int                                 nSubsteps;

  public:
    /**
     * @struct StateSensitivities
     * @brief Sensitivity matrices required for analytical substepping.
     * All matrices represent derivatives of flattened arrays and are stored in Eigen's column-major format.
     */
    struct StateSensitivities {
      /** @brief Jacobian of State_new w.r.t Deformation Gradient F_new.
       * Size: nStateVars x 9
       */
      Eigen::MatrixXd dState_dF;

      /** @brief Jacobian of State_new w.r.t State_old.
       * Size: nStateVars x nStateVars
       */
      Eigen::MatrixXd dState_dStateOld;

      /** @brief Jacobian of Stress w.r.t State_old.
       * Size: 9 x nStateVars
       */
      Eigen::MatrixXd dStress_dStateOld;
    };

    MarmotMaterialFiniteStrainSubstepped( const double* matProperties_, int nMaterialProperties_, int materialNumber_ )
      : MarmotMaterialFiniteStrain( matProperties_, nMaterialProperties_, materialNumber_ )
    {
      // 1. Parse Substepping parameters
      if ( nMaterialProperties_ < 1 ) {
        throw std::invalid_argument( "MarmotMaterialFiniteStrainSubstepped requires at least 1 property (nSubsteps)" );
      }
      this->nSubsteps = static_cast< int >( matProperties_[0] );
      if ( this->nSubsteps < 1 )
        this->nSubsteps = 1;

      baseMaterial = std::make_unique< BaseMaterialType >( matProperties_ + 1,
                                                           nMaterialProperties_ - 1,
                                                           materialNumber_ );
      initializeStateLayout();
    }

    virtual ~MarmotMaterialFiniteStrainSubstepped() = default;

    double getDensity( const double* stateVars ) const override { return baseMaterial->getDensity( stateVars ); }

    void initializeStateLayout()
    {

      // Add our state variable: Deformation Gradient at start of GLOBAL step (F_n)
      this->stateLayout.add( "Substepping_F_n", 9 );

      // Add Base Material state variables
      this->stateLayout.add( "materialstate", baseMaterial->getNumberOfRequiredStateVars() );

      this->stateLayout.finalize();
    }

    void initializeYourself( double* stateVars, int nStateVars ) override
    {
      // 1. Initialize Base Material part
      int baseVarsCount = baseMaterial->getNumberOfRequiredStateVars();

      baseMaterial->initializeYourself( stateLayout.getPtr( stateVars, "materialstate" ), baseVarsCount );

      FastorStandardTensors::Tensor33d&
        Fn = this->stateLayout.getAs< FastorStandardTensors::Tensor33d& >( stateVars, "Substepping_F_n" );
      memcpy( Fn.data(), FastorStandardTensors::Spatial3D::I.data(), 9 * sizeof( double ) );
    }

    /**
     * @brief Extended stress update computing sensitivities for substepping.
     * * The default implementation uses Forward Finite Differences.
     * @param response Constitutive response to be computed.
     * @param tangents Algorithmic moduli to be computed.
     * @param sensitivities Sensitivity matrices to be computed.
     * @param deformation Current deformation state.
     * @param timeIncrement Time increment information.
     *
     * @note This function can be overridden in derived classes to provide analytical sensitivities.
     * @note The base material's computeStress function is called within this function.
     */
    virtual void computeStressWithSensitivities( ConstitutiveResponse< 3 >& response,
                                                 AlgorithmicModuli< 3 >&    tangents,
                                                 StateSensitivities&        sensitivities,
                                                 const Deformation< 3 >&    deformation,
                                                 const TimeIncrement&       timeIncrement ) const
    {
      using namespace Eigen;
      using namespace Marmot::NumericalAlgorithms::Differentiation;

      int nState = baseMaterial->getNumberOfRequiredStateVars();

      std::vector< double > stateOld( nState );
      std::memcpy( stateOld.data(), response.stateVars, nState * sizeof( double ) );

      baseMaterial->computeStress( response, tangents, deformation, timeIncrement );

      if ( nState == 0 )
        return;

      sensitivities.dState_dF.resize( nState, 9 );
      sensitivities.dState_dStateOld.resize( nState, nState );
      sensitivities.dStress_dStateOld.resize( 9, nState );

      auto func_dStressAndState_dF = [&]( const VectorXd& F_vec ) -> VectorXd {
        Deformation< 3 > defPerturbed;
        std::memcpy( defPerturbed.F.data(), F_vec.data(), 9 * sizeof( double ) );

        ConstitutiveResponse< 3 > respPert;
        std::vector< double >     stateTemp = stateOld;
        respPert.stateVars                  = stateTemp.data();

        AlgorithmicModuli< 3 > tanPert;

        baseMaterial->computeStress( respPert, tanPert, defPerturbed, timeIncrement );

        // return [stress; state_new] as a single vector
        VectorXd res( 9 + nState );
        // Fill Stress part
        std::memcpy( res.data(), respPert.tau.data(), 9 * sizeof( double ) );
        // Fill State_new part
        std::memcpy( res.data() + 9, respPert.stateVars, nState * sizeof( double ) );
        return res;
      };

      auto func_dStressAndState_dStateOld = [&]( const VectorXd& StateOld_vec ) -> VectorXd {
        ConstitutiveResponse< 3 > respPert;
        std::vector< double >     stateTemp( nState );
        std::memcpy( stateTemp.data(), StateOld_vec.data(), nState * sizeof( double ) );
        respPert.stateVars = stateTemp.data();

        AlgorithmicModuli< 3 > tanPert; // dummy

        // Run update with perturbed old state
        baseMaterial->computeStress( respPert, tanPert, deformation, timeIncrement );

        // return [stress; state_new] as a single vector
        VectorXd res( 9 + nState );
        // Fill Stress part
        std::memcpy( res.data(), respPert.tau.data(), 9 * sizeof( double ) );
        // Fill State_new part
        std::memcpy( res.data() + 9, respPert.stateVars, nState * sizeof( double ) );

        return res;
      };

      // Compute dState/dF_new
      const MatrixXd dStressAndState_dF = forwardDifference( func_dStressAndState_dF,
                                                             Map< const VectorXd >( deformation.F.data(), 9 ) );

      // Extract dStress/dF_new
      MatrixXd dStress_dF_flat = dStressAndState_dF.block( 0, 0, 9, 9 );

      // possibiblity to take fd tangent directly
      tangents.dTau_dF = reshape< 3, 3, 3, 3 >(
        Fastor::transpose( FastorStandardTensors::Tensor99d( dStress_dF_flat.data() ) ) );

      // Extract dState/dF_new
      sensitivities.dState_dF = dStressAndState_dF.block( 9, 0, nState, 9 );

      // Compute [dStress/dStateOld; dState/dStateOld]
      const MatrixXd dStressAndState_dStateOld = forwardDifference( func_dStressAndState_dStateOld,
                                                                    Map< const VectorXd >( stateOld.data(), nState ) );

      // Extract dState/dStateOld
      sensitivities.dState_dStateOld = dStressAndState_dStateOld.block( 9, 0, nState, nState );

      // Extract dStress/dStateOld
      sensitivities.dStress_dStateOld = dStressAndState_dStateOld.block( 0, 0, 9, nState );
    }

    /**
     * @brief Compute stress with time substepping and analytical tangent accumulation.
     * @param response Constitutive response to be computed.
     * @param tangents Algorithmic moduli to be computed.
     * @param deformation Current deformation state.
     * @param timeIncrement Time increment information.
     *
     * @details The deformation increment is divided into nSubsteps sub-increments.
     * For each sub-increment, the base material's stress update is called,
     * and the sensitivities are used to accumulate the overall tangent.
     */
    void computeStress( ConstitutiveResponse< 3 >& response,
                        AlgorithmicModuli< 3 >&    tangents,
                        const Deformation< 3 >&    deformation,
                        const TimeIncrement&       timeIncrement ) const override
    {
      using namespace Eigen;
      using namespace FastorStandardTensors;

      Tensor33d&      Fn_ref = this->stateLayout.getAs< Tensor33d& >( response.stateVars, "Substepping_F_n" );
      const Tensor33d Fn     = Fn_ref;
      const Tensor33d Fn1    = deformation.F;

      int nBaseState = baseMaterial->getNumberOfRequiredStateVars();

      // J_accum = d(State_i)/d(F_global). Size: nState x 9
      MatrixXd J_accum = MatrixXd::Zero( nBaseState, 9 );

      double dt_sub = timeIncrement.dT / static_cast< double >( nSubsteps );
      double t_curr = timeIncrement.time;

      ConstitutiveResponse< 3 > subResponse = { Tensor33d( 0.0 ), 0.0, nullptr };
      subResponse.stateVars                 = this->stateLayout.getPtr( response.stateVars, "materialstate" );

      AlgorithmicModuli< 3 > subTangents      = { Tensor3333d( 0.0 ) };
      StateSensitivities     subSensitivities = { MatrixXd(), MatrixXd(), MatrixXd() };
      Deformation< 3 >       subDef           = { Tensor33d( 0.0 ) };

      for ( int i = 1; i <= nSubsteps; ++i ) {
        double xi = static_cast< double >( i ) / static_cast< double >( nSubsteps );

        // Linear Interpolation: F(xi) = (1-xi)*F_n + xi*F_n1
        subDef.F                    = ( 1.0 - xi ) * Fn + xi * Fn1;
        const double dFsub_dFglobal = xi;

        TimeIncrement subTime = { t_curr, dt_sub };

        computeStressWithSensitivities( subResponse, subTangents, subSensitivities, subDef, subTime );

        // Update accumulated Jacobian for history
        if ( i < nSubsteps ) {
          if ( nBaseState > 0 ) {
            J_accum = subSensitivities.dState_dF * dFsub_dFglobal + subSensitivities.dState_dStateOld * J_accum;
          }
        }
        else {
          // Final Step: Compute overall Tangent dTau/dF_global
          MatrixXd C_hist = MatrixXd::Zero( 9, 9 );
          if ( nBaseState > 0 ) {
            C_hist = subSensitivities.dStress_dStateOld * J_accum;
          }

          // build global tangent
          // transpose because second term is from Eigen's world with column-major layout
          Fastor::Tensor< double, 9, 9 > C_global = subTangents.dTau_dF * dFsub_dFglobal +
                                                    transpose( Fastor::Tensor< double, 9, 9 >( C_hist.data() ) );

          std::memcpy( tangents.dTau_dF.data(), C_global.data(), 81 * sizeof( double ) );
        }
        t_curr += dt_sub;
      }

      // set final response
      response.tau                  = subResponse.tau;
      response.elasticEnergyDensity = subResponse.elasticEnergyDensity;

      // Update stored F_n to F_n1
      Fn_ref = deformation.F;
    }
  };

} // namespace Marmot::Materials
