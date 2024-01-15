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
#include "Marmot/MarmotAutomaticDifferentiationForFastor.h"
#include "Marmot/MarmotDeformationMeasures.h"
#include "Marmot/MarmotEnergyDensityFunctions.h"
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include "Marmot/MarmotMicromorphicTensorBasics.h"
#include "Marmot/MarmotMicropolarGeneralizedInvariants.h"
#include "Marmot/MarmotMicropolarPlasticity.h"
#include "Marmot/MarmotStateVarVectorManager.h"
#include "Marmot/MarmotTypedefs.h"
#include <string>

namespace Marmot::Materials {

  using namespace Fastor;
  using namespace FastorStandardTensors;
  using namespace FastorIndices;

  class FiniteStrainJ2Plasticity : public MarmotMaterialFiniteStrain {
  public:
    using MarmotMaterialFiniteStrain::MarmotMaterialFiniteStrain;

    // elastic constants
    const double K, G;

    // plasticity parameters
    const double fy, fyInf, eta, H;

    // implementation
    const int implementationType;

    FiniteStrainJ2Plasticity( const double* materialProperties, int nMaterialProperties, int materialLabel );

    void computeStress( ConstitutiveResponse< 3 >& response,
                        AlgorithmicModuli< 3 >&    tangents,
                        const Deformation< 3 >&    deformation,
                        const TimeIncrement&       timeIncrement );
    void computeStressWithScalarReturnMapping( ConstitutiveResponse< 3 >& response,
                                               AlgorithmicModuli< 3 >&    tangents,
                                               const Deformation< 3 >&    deformation,
                                               const TimeIncrement&       timeIncrement );
    void computeStressWithFullReturnMapping( ConstitutiveResponse< 3 >& response,
                                             AlgorithmicModuli< 3 >&    tangents,
                                             const Deformation< 3 >&    deformation,
                                             const TimeIncrement&       timeIncrement );
    int  getNumberOfRequiredStateVars() { return FiniteStrainJ2PlasticityStateVarManager::layout.nRequiredStateVars; }

    class FiniteStrainJ2PlasticityStateVarManager : public MarmotStateVarVectorManager {

    public:
      inline const static auto layout = makeLayout( {
        { .name = "Fp", .length = 9 },
        { .name = "alphaP", .length = 1 },
      } );

      Fastor::TensorMap< double, 3, 3 > Fp;
      double&                           alphaP;

      FiniteStrainJ2PlasticityStateVarManager( double* theStateVarVector )
        : MarmotStateVarVectorManager( theStateVarVector, layout ), Fp( &find( "Fp" ) ), alphaP( find( "alphaP" ) ){};
    };
    std::unique_ptr< FiniteStrainJ2PlasticityStateVarManager > stateVars;

    void assignStateVars( double* stateVars, int nStateVars );

    StateView getStateView( const std::string& result );

    void initializeYourself();

    std::tuple< double, Tensor33d, double > yieldFunction( const Tensor33d& Fe, const double betaP )
    {

      Tensor33d   mandelStress;
      Tensor3333d dMandel_dFe;
      std::tie( mandelStress, dMandel_dFe ) = computeMandelStress( Fe );

      double      f;
      Tensor33d   df_dMandel;
      Tensor3333d d2f_dMandel_dMandel;
      Tensor33d   df_dFe;
      double      df_dBetaP;

      std::tie( f, df_dMandel, d2f_dMandel_dMandel, df_dBetaP ) = yieldFunctionFromStress( mandelStress, betaP );

      df_dFe = einsum< IJ, IJKL >( df_dMandel, dMandel_dFe );

      return { f, df_dFe, df_dBetaP };
    }
    template < typename T >
    T yieldFunctionFromStress( const Tensor33t< T >& mandelStress, const T betaP )
    {
      Tensor33t< T > dev = ContinuumMechanics::Micropolar::GeneralizedInvariants::deviatoric( mandelStress );
      T              rho = sqrt( Fastor::inner( dev, dev ) );
      if ( double( rho ) == 0.0 )
        rho += 1e-15;
      T f = 1. / fy * ( rho - Marmot::Constants::sqrt2_3 * betaP );
      return f;
    }

    std::tuple< double, Tensor33d, Tensor3333d, double > yieldFunctionFromStress( const Tensor33d& mandelStress,
                                                                                  const double     betaP )
    {
      Tensor33d    dev = ContinuumMechanics::Micropolar::GeneralizedInvariants::deviatoric( mandelStress );
      const double rho = std::max( sqrt( Fastor::inner( dev, dev ) ), 1e-15 );
      double       f   = 1. / fy * ( rho - Marmot::Constants::sqrt2_3 * betaP );

      Tensor33d dRho_dMandel = 1. / rho * dev;

      /* Tensor3333d d2Rho_dMandel_dMandel = -1.0 / std::pow( rho, 3 ) * outer( dev, dev ) + */
      /*                                     1. / rho * Spatial3D::Deviatoric; */

      /* Tensor3333d d2f_dMandel_dMandel_analytical = 1. / fy * d2Rho_dMandel_dMandel; */
      Tensor33d df_dMandel_analytical = 1. / fy * dRho_dMandel;
      Tensor33d df_dMandel;

      Tensor3333d d2f_dMandel_dMandel;

      // Workaround using autodiff for 2nd derivative of f wrt mandel stress
      std::tie( f, df_dMandel, d2f_dMandel_dMandel ) = AutomaticDifferentiation::SecondOrder::d2f_dTensor_dTensor<
        3 >( [&]( const FastorStandardTensors::Tensor33t< autodiff::dual2nd >&
                    mandel_ ) { return yieldFunctionFromStress( mandel_, autodiff::dual2nd( betaP ) ); },
             mandelStress );
      /* std::cout << "diff 1st: " << df_dMandel_analytical - df_dMandel << std::endl; */
      /* std::cout << "diff 2nd: " << d2f_dMandel_dMandel_analytical - d2f_dMandel_dMandel << std::endl; */

      return { f, df_dMandel_analytical, d2f_dMandel_dMandel, -Constants::sqrt2_3 / fy };
    }

    bool isYielding( const Tensor33d& Fe, const double betaP )
    {
      double    f, df_dBetaP;
      Tensor33d df_dFe;
      std::tie( f, df_dFe, df_dBetaP ) = yieldFunction( Fe, betaP );
      if ( f > 0.0 )
        return true;
      else
        return false;
    }

    std::tuple< Tensor33d, Tensor3333d > computeMandelStress( const Tensor33d& Fe )
    {
      using namespace Marmot::ContinuumMechanics;
      Tensor33d   Ce;
      Tensor3333d dCe_dFe;
      std::tie( Ce, dCe_dFe ) = DeformationMeasures::FirstOrderDerived::CauchyGreen( Fe );

      double      psi_;
      Tensor33d   dPsi_dCe;
      Tensor3333d d2Psi_dCedCe, dMandel_dCe;

      std::tie( psi_, dPsi_dCe, d2Psi_dCedCe ) = EnergyDensityFunctions::SecondOrderDerived::PenceGouPotentialB( Ce,
                                                                                                                 K,
                                                                                                                 G );
      Tensor33d PK2                            = 2.0 * dPsi_dCe;
      /* const Tensor33d mandel = Ce % PK2; */
      const Tensor33d mandel = einsum< Ii, iJ >( Ce, PK2 );
      dMandel_dCe            = einsum< Ii, iJKL, to_IJKL >( Ce, 2. * d2Psi_dCedCe ) +
                    einsum< IK, iL, iJ, to_IJKL >( Spatial3D::I, Spatial3D::I, PK2 );
      Tensor3333d dMandel_dFe = einsum< IJKL, KLMN >( dMandel_dCe, dCe_dFe );
      return { mandel, dMandel_dFe };
    }

    std::tuple< double, double > computeBetaP( const double alphaP )
    {
      const double beta           = fyInf + ( fy - fyInf ) * exp( -alphaP * eta ) + alphaP * H;
      const double dBetaP_dAlphaP = -eta * ( fy - fyInf ) * exp( -alphaP * eta ) + H;
      return { beta, dBetaP_dAlphaP };
    }

    std::tuple< double, double > g( const Tensor33d& Fe_trial, const Tensor33d& df_dS, double dLambda, double betaP )
    {
      Tensor33d dFp, ddFp_ddLambda;
      std::tie( dFp, ddFp_ddLambda ) = computePlasticIncrement( df_dS, dLambda );

      double    f, df_dBetaP;
      Tensor33d df_dFe;

      Tensor33d   dFpInv   = Fastor::inverse( dFp );
      Tensor3333d dFe_ddFp = einsum< iI, iJKL >( Fe_trial, einsum< LI, JK, to_IJKL >( -dFpInv, dFpInv ) );

      std::tie( f, df_dFe, df_dBetaP ) = yieldFunction( Tensor33d( Fe_trial % dFpInv ), betaP );

      double df_ddLambda = einsum< IJ, IJKL, KL >( df_dFe, dFe_ddFp, ddFp_ddLambda ).toscalar();

      return { f, df_ddLambda };
    }

    template < typename T >
    Tensor33t< T > computeReturnMappingDirection( Tensor33t< T > Fe, T betaP )
    {

      Tensor33d   mandelStressTrial;
      Tensor3333d dMandel_dFe;
      std::tie( mandelStressTrial, dMandel_dFe ) = computeMandelStress( Fe );
      double    f, df_dBetaP;
      Tensor33d df_dS;
      std::tie( f, df_dS, df_dBetaP ) = yieldFunctionFromStress( mandelStressTrial, betaP );

      return df_dS;
    }

    template < typename T >
    std::tuple< Tensor33d, Tensor33d > computePlasticIncrement( Tensor33t< T > df_dS, T dLambda )
    {

      Tensor33d   dGp = multiplyFastorTensorWithScalar( df_dS, dLambda );
      Tensor33d   dFp;
      Tensor3333d ddFp_ddGp;
      std::tie( dFp, ddFp_ddGp ) = ContinuumMechanics::Micropolar::Plasticity::FlowIntegration::FirstOrderDerived::
        exponentialMap( dGp );

      return { dFp, einsum< IJKL, KL >( ddFp_ddGp, df_dS ) };
    }

    std::tuple< Eigen::VectorXd, Eigen::MatrixXd > computeResidualVector( const Eigen::VectorXd& X,
                                                                          const Tensor33d&       FeTrial,
                                                                          const double           alphaPTrial )
    {

      const int idxA = 9;
      const int idxF = 10;
      using namespace Eigen;
      using mV9d = Eigen::Map< Eigen::Matrix< double, 9, 1 > >;
      using mM9d = Eigen::Map< Eigen::Matrix< double, 9, 9 > >;
      VectorXd R( 11 );
      MatrixXd dR_dX( 11, 11 );
      dR_dX.setZero();
      // initialize residual
      R.segment< 9 >( 0 ) = -mV9d( FeTrial.data() );
      R( 9 )              = -alphaPTrial;

      Tensor33d Fe( X.segment( 0, 9 ).data() );

      const double dLambda = X( 10 );
      const double alphaP  = X( 9 );

      /* std::cout << "dLambda: " << dLambda << std::endl; */
      double betaP, dBetaP_dAlphaP;
      std::tie( betaP, dBetaP_dAlphaP ) = computeBetaP( alphaP );

      // compute mandel stress
      Tensor33d   mandelStress;
      Tensor3333d dMandel_dFe, d2f_dMandel_dMandel;
      std::tie( mandelStress, dMandel_dFe ) = computeMandelStress( Fe );

      double    f, df_dBetaP;
      Tensor33d df_dMandel;
      std::tie( f, df_dMandel, d2f_dMandel_dMandel, df_dBetaP ) = yieldFunctionFromStress( mandelStress, betaP );

      Tensor33d   dGp = dLambda * df_dMandel;
      Tensor33d   dFp;
      Tensor3333d ddFp_ddGp;
      std::tie( dFp, ddFp_ddGp ) = ContinuumMechanics::Micropolar::Plasticity::FlowIntegration::FirstOrderDerived::
        exponentialMap( dGp );

      Tensor3333d ddGp_dFe      = dLambda * einsum< ijmn, mnkL >( d2f_dMandel_dMandel, dMandel_dFe );
      Tensor33d   ddFp_ddLambda = einsum< IJKL, KL >( ddFp_ddGp, df_dMandel );

      Tensor3333d ddFp_dFe     = einsum< iImn, mnkL >( ddFp_ddGp, ddGp_dFe );
      Tensor3333d dFeTrial_dFe = einsum< iIkL, iJ, to_IJkL >( ddFp_dFe, Fe ) +
                                 einsum< iI, ik, JL, to_IJkL >( dFp, Spatial3D::I, Spatial3D::I );

      /* std::cout << "ddFeTrial_dFe: " << dFeTrial_dFe << std::endl; */

      Tensor33d dFe_ddLambda = einsum< iI, iJ >( ddFp_ddLambda, Fe );

      df_dBetaP = -Constants::sqrt2_3 / fy;

      Tensor33d df_dFe                 = einsum< IJ, IJKL >( df_dMandel, dMandel_dFe );
      std::tie( f, df_dFe, df_dBetaP ) = yieldFunction( Fe, betaP );

      R.segment< 9 >( 0 ) += mV9d( Tensor33d( einsum< iI, iJ >( dFp, Fe ) ).data() );
      R( idxA ) += ( alphaP + dLambda * df_dBetaP );
      R( idxF ) = f;

      dR_dX.block< 9, 9 >( 0, 0 ) = mM9d( dFeTrial_dFe.data() ).transpose();

      dR_dX.block< 1, 9 >( idxF, 0 ) = mV9d( df_dFe.data() );
      dR_dX.block< 9, 1 >( 0, idxF ) = mV9d( dFe_ddLambda.data() );
      dR_dX( idxF, idxA )            = df_dBetaP * dBetaP_dAlphaP;
      dR_dX( idxA, idxF )            = df_dBetaP;
      dR_dX( idxA, idxA )            = 1.0;

      return { R, dR_dX };
    }
  };

} // namespace Marmot::Materials
