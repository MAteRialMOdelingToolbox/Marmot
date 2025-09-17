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
#include "Marmot/MarmotDeformationMeasures.h"
#include "Marmot/MarmotEnergyDensityFunctions.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotFiniteStrainPlasticity.h"
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotStateVarVectorManager.h"
#include "Marmot/MarmotTypedefs.h"
#include <string>

namespace Marmot::Materials {

  using namespace Fastor;
  using namespace FastorStandardTensors;
  using namespace FastorIndices;

  /**
   * \class Marmot::Materials::FiniteStrainJ2Plasticity
   * \brief Finite-strain J2 plasticity with isotropic hardening, multiplicative split
   * \f$\mathbf{F}=\mathbf{F}_e\mathbf{F}_p\f$.
   *
   * The elastic response uses the compressible Neo-Hookean Pence–Gou potential (variant B) with \c K and \c G.
   * Plasticity is formulated in Mandel stress with an associated flow rule. The yield function is
   * \f[
   *   f = \frac{\rho}{f_y} - \sqrt{\frac{2}{3}}\,\frac{\beta_p}{f_y},
   * \f]
   * where \f$\rho=\|\mathrm{dev}(\mathbf{M})\|\f$ and \f$\beta_p\f$ is the current flow stress (isotropic hardening).
   *
   * \par Material parameters (indices in #materialProperties)
   * - \c K   (#materialProperties[0]) — bulk modulus [Pa]
   * - \c G   (#materialProperties[1]) — shear modulus [Pa]
   * - \c fy  (#materialProperties[2]) — initial yield stress (reference) [Pa]
   * - \c fyInf (#materialProperties[3]) — saturated (asymptotic) flow stress [Pa]
   * - \c eta (#materialProperties[4]) — saturation rate parameter [-]
   * - \c H   (#materialProperties[5]) — linear hardening modulus [Pa]
   * - \c implementationType (#materialProperties[6]) — algorithm selector (see below)
   * - \c density (optional, #materialProperties[7]) — mass density [kg/m^3]
   *
   * \par State variables
   * - \c Fp (9 components) — plastic deformation gradient
   * - \c alphaP (1) — accumulated plastic strain (hardening variable)
   *
   * \par Implementation variants (implementationType)
   * - 0: scalar return mapping (not implemented)
   * - 1: full return mapping with analytic linearization for \c dTau_dF
   * - 2: FDAF — finite-difference (forward) linearization of the residual
   * - 3: FDAC — finite-difference (central) linearization of the residual
   * - 4: CSDA — complex-step differentiation (forward) of the residual
   *
   * \ingroup materials_plasticity
   */

  class FiniteStrainJ2Plasticity : public MarmotMaterialFiniteStrain {
  public:
    using MarmotMaterialFiniteStrain::MarmotMaterialFiniteStrain;

    // elastic constants
    const double K, G; /**< Elastic moduli: bulk \c K and shear \c G (from #materialProperties[0..1]). */

    // plasticity parameters
    const double fy, fyInf, eta, H; /**< Isotropic hardening parameters (from #materialProperties[2..5]). */

    // implementation
    const int implementationType; /**< Algorithm variant selector (from #materialProperties[6]). */

    // mass properties;
    const double density; /**< Mass density (optional; #materialProperties[7] if provided). */

    /**
     * \brief Construct the finite-strain J2 plasticity model.
     * \param materialProperties Array with parameters: K, G, fy, fyInf, eta, H, implementationType, [density].
     * \param nMaterialProperties Length of \p materialProperties.
     * \param materialLabel User-defined material label (passed to base).
     */
    FiniteStrainJ2Plasticity( const double* materialProperties, int nMaterialProperties, int materialLabel );

    /**
     * \brief Compute stress and tangent for the current step (dispatches by \c implementationType).
     * \param[out] response  Populated with \c tau, \c elasticEnergyDensity, \c rho.
     * \param[out] tangents  Populated with \c dTau_dF (algorithmic tangent).
     * \param[in]  deformation Uses \c deformation.F.
     * \param[in]  timeIncrement Current time step information.
     * \throws std::invalid_argument If \c implementationType is not supported.
     */
    void computeStress( ConstitutiveResponse< 3 >& response,
                        AlgorithmicModuli< 3 >&    tangents,
                        const Deformation< 3 >&    deformation,
                        const TimeIncrement&       timeIncrement );

    /** \brief Scalar-return mapping variant (not yet implemented). */

    void computeStressWithScalarReturnMapping( ConstitutiveResponse< 3 >& response,
                                               AlgorithmicModuli< 3 >&    tangents,
                                               const Deformation< 3 >&    deformation,
                                               const TimeIncrement&       timeIncrement );

    /** \brief Full return mapping; derivatives computed analytically. */

    void computeStressWithFullReturnMapping( ConstitutiveResponse< 3 >& response,
                                             AlgorithmicModuli< 3 >&    tangents,
                                             const Deformation< 3 >&    deformation,
                                             const TimeIncrement&       timeIncrement );

    /** \brief Full return mapping; derivatives approximated via forward finite differences. */

    void computeStressFDAF( ConstitutiveResponse< 3 >& response,
                            AlgorithmicModuli< 3 >&    tangents,
                            const Deformation< 3 >&    deformation,
                            const TimeIncrement&       timeIncrement );

    /** \brief Full return mapping; derivatives approximated via central finite differences. */

    void computeStressFDAC( ConstitutiveResponse< 3 >& response,
                            AlgorithmicModuli< 3 >&    tangents,
                            const Deformation< 3 >&    deformation,
                            const TimeIncrement&       timeIncrement );

    /** \brief Full return mapping; derivatives via complex-step differentiation. */

    void computeStressCSDA( ConstitutiveResponse< 3 >& response,
                            AlgorithmicModuli< 3 >&    tangents,
                            const Deformation< 3 >&    deformation,
                            const TimeIncrement&       timeIncrement );

    /** \brief Number of required state variables (\c Fp:9 + \c alphaP:1 = 10). */

    int getNumberOfRequiredStateVars() { return FiniteStrainJ2PlasticityStateVarManager::layout.nRequiredStateVars; }

    /** \brief Return the material density if provided in the material parameters */

    double getDensity() { return density; }

    /** \brief State variable manager for \c Fp and \c alphaP. */

    class FiniteStrainJ2PlasticityStateVarManager : public MarmotStateVarVectorManager {

    public:
      inline const static auto layout = makeLayout( {
        { .name = "Fp", .length = 9 },
        { .name = "alphaP", .length = 1 },
      } );

      Fastor::TensorMap< double, 3, 3 > Fp;
      double&                           alphaP;

      FiniteStrainJ2PlasticityStateVarManager( double* theStateVarVector )
        : MarmotStateVarVectorManager( theStateVarVector, layout ), Fp( &find( "Fp" ) ), alphaP( find( "alphaP" ) ) {};
    };
    std::unique_ptr< FiniteStrainJ2PlasticityStateVarManager > stateVars;

    /** \brief Bind external storage for internal variables (\c Fp, \c alphaP).
     *  \throws std::invalid_argument if storage is insufficient.
     */

    void assignStateVars( double* stateVars, int nStateVars );

    /** \brief Expose a named view into the state vector.
     *  \param result One of: \c "Fp" (length 9), \c "alphaP" (length 1).
     *  \throws std::out_of_range if \p result is unknown.
     */

    StateView getStateView( const std::string& result );

    /** \brief Initialize (\c Fp = I). */

    void initializeYourself();

    /** \brief Evaluate yield function and its first derivatives w.r.t. \c Fe and \c betaP.
     *  \param Fe     Elastic deformation gradient.
     *  \param betaP  Hardening parameter.
     */

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

    /**
     * \brief Yield function \c f(\c mandelStress,\c betaP) w.r.t. Mandel stress (templated).
     * \tparam T Scalar type (double, autodiff dual, complex).
     * \param mandelStress Mandel stress \c M.
     * \param betaP        Hardening measure.
     */

    template < typename T >
    T yieldFunctionFromStress( const Tensor33t< T >& mandelStress, const T betaP )
    {
      const Tensor33t< T > dev = deviatoric( mandelStress );
      T                    rho = sqrt( Fastor::inner( dev, dev ) );
      if ( double( rho ) == 0.0 )
        rho += 1e-15;
      const T f = 1. / fy * ( rho - Marmot::Constants::sqrt2_3 * betaP );
      return f;
    }

    /** \brief Yield function w.r.t. Mandel stress.
     *  \param mandelStress Mandel stress \c M.
     *  \param betaP        Hardening measure.
     */

    std::tuple< double, Tensor33d, Tensor3333d, double > yieldFunctionFromStress( const Tensor33d& mandelStress,
                                                                                  const double     betaP )
    {
      Tensor33d    dev = deviatoric( mandelStress );
      const double rho = std::max( sqrt( Fastor::inner( dev, dev ) ), 1e-15 );
      const double f   = 1. / fy * ( rho - Marmot::Constants::sqrt2_3 * betaP );

      Tensor33d dRho_dMandel = 1. / rho * dev;

      Tensor3333d d2Rho_dMandel_dMandel = -1.0 / std::pow( rho, 3 ) * outer( dev, dev ) +
                                          1. / rho *
                                            ( einsum< ik, jl, to_ijkl >( Spatial3D::I, Spatial3D::I ) -
                                              1. / 3 *
                                                einsum< ij, IK, IL, to_ijKL >( Spatial3D::I,
                                                                               Spatial3D::I,
                                                                               Spatial3D::I ) );

      Tensor3333d d2f_dMandel_dMandel = 1. / fy * d2Rho_dMandel_dMandel;
      Tensor33d   df_dMandel          = 1. / fy * dRho_dMandel;

      return { f, df_dMandel, d2f_dMandel_dMandel, -Constants::sqrt2_3 / fy };
    }

    /**
     * \brief First-order derivatives of \c f with respect to Mandel stress.
     * \tparam T Scalar type (double, autodiff dual, complex, ...).
     * \return \c { f, df_dMandel, df_dBetaP }.
     */

    template < typename T >
    std::tuple< T, Tensor33t< T >, T > yieldFunctionFromStressFirstOrderDerived( const Tensor33t< T >& mandelStress,
                                                                                 const T               betaP )
    {
      Tensor33t< T > dev = deviatoric( mandelStress );
      T              rho = sqrt( Fastor::inner( dev, dev ) );
      if ( Math::makeReal( rho ) == 0.0 )
        rho += 1e-15;
      const T f = 1. / fy * ( rho - Marmot::Constants::sqrt2_3 * betaP );

      const Tensor33t< T > dRho_dMandel = multiplyFastorTensorWithScalar( dev, 1. / rho );

      const Tensor33t< T > df_dMandel = multiplyFastorTensorWithScalar( dRho_dMandel, T( 1. / fy ) );

      return { f, df_dMandel, T( -Constants::sqrt2_3 / fy ) };
    }

    /** \brief Check yield condition \c f(Fe,betaP) > 0.
     */

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

    /** \brief Compute Mandel stress \c M = C_e : PK2(C_e) and its derivative w.r.t. \c Fe. */
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
      Tensor33d       PK2                      = 2.0 * dPsi_dCe;
      const Tensor33d mandel                   = Ce % PK2;
      /* const Tensor33d mandel = einsum< Ii, iJ >( Ce, PK2 ); */
      dMandel_dCe = einsum< Ii, iJKL, to_IJKL >( Ce, 2. * d2Psi_dCedCe ) +
                    einsum< IK, iL, iJ, to_IJKL >( Spatial3D::I, Spatial3D::I, PK2 );
      Tensor3333d dMandel_dFe = einsum< IJKL, KLMN >( dMandel_dCe, dCe_dFe );
      return { mandel, dMandel_dFe };
    }

    /** \brief Compute Mandel stress only (templated scalar type). */

    template < typename T >
    Tensor33t< T > computeMandelStressOnly( const Tensor33t< T >& Fe )
    {
      using namespace Marmot::ContinuumMechanics;
      Tensor33t< T > Ce = DeformationMeasures::CauchyGreen( Fe );

      T              psi_;
      Tensor33t< T > dPsi_dCe;

      std::tie( psi_, dPsi_dCe )  = EnergyDensityFunctions::FirstOrderDerived::PenceGouPotentialB( Ce, K, G );
      Tensor33t< T >       PK2    = multiplyFastorTensorWithScalar( dPsi_dCe, T( 2.0 ) );
      const Tensor33t< T > mandel = Ce % PK2;
      return mandel;
    }

    /** \brief Templated version of the hardening law returning \c beta only. */

    template < typename T >
    T computeBetaPOnly( const T alphaP )
    {
      const T beta = fyInf + ( fy - fyInf ) * exp( -alphaP * eta ) + alphaP * H;
      return beta;
    }

    /**  \brief Isotropic hardening law
     * \f$\beta(\alpha_p)=f_{y\infty}+(f_y-f_{y\infty})e^{-\eta\alpha_p}+H\,\alpha_p\f$. And its derivative w.r.t. \c
     * alphaP.
     */

    std::tuple< double, double > computeBetaP( const double alphaP )
    {
      const double beta           = fyInf + ( fy - fyInf ) * exp( -alphaP * eta ) + alphaP * H;
      const double dBetaP_dAlphaP = -eta * ( fy - fyInf ) * exp( -alphaP * eta ) + H;
      return { beta, dBetaP_dAlphaP };
    }

    /** \brief Scalar consistency function.
     */

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

    /** \brief Flow direction for return mapping, \c df/dS (templated).
     *  \tparam T   Scalar type (double, AD, complex).
     *  \param Fe   Elastic deformation gradient.
     *  \param betaP Hardening measure.
     */
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

    /** \brief Incremental plastic flow via exponential map.
     * \tparam T Scalar type (double, autodiff dual, complex, ...).
     * \param df_dS Flow direction in Mandel stress space.
     * \param dLambda Plastic multiplier increment.
     */

    template < typename T >
    std::tuple< Tensor33d, Tensor33d > computePlasticIncrement( Tensor33t< T > df_dS, T dLambda )
    {

      Tensor33d   dGp = multiplyFastorTensorWithScalar( df_dS, dLambda );
      Tensor33d   dFp;
      Tensor3333d ddFp_ddGp;
      std::tie( dFp, ddFp_ddGp ) = ContinuumMechanics::FiniteStrain::Plasticity::FlowIntegration::FirstOrderDerived::
        exponentialMap( dGp );

      return { dFp, einsum< IJKL, KL >( ddFp_ddGp, df_dS ) };
    }

    /** \brief Residual vector \c R for the full return mapping (templated).
     *  \details State vector \c X has 11 unknowns: 9 for \c Fe, 1 for \c alphaP, 1 for \c dLambda.
     *           Residual components:
     *            - \c R[0..8] : elastic part consistency (updated \c Fe vs trial),
     *            - \c R[9]    : \c alphaP update,
     *            - \c R[10]   : yield function \c f.
     *  \tparam T         Scalar type (double/dual/complex).
     *  \param X          Current iterate of unknowns \c [Fe(9), alphaP, dLambda]^T.
     *  \param FeTrial    Trial elastic deformation gradient.
     *  \param alphaPTrial Trial accumulated plastic strain.
     */

    template < typename T >
    VectorXt< T > computeResidualVector( const VectorXt< T >& X, const Tensor33d& FeTrial, const double alphaPTrial )
    {

      const int idxA = 9;
      const int idxF = 10;
      using namespace Eigen;
      using mV9t = Eigen::Map< Eigen::Matrix< T, 9, 1 > >;
      VectorXt< T > R( 11 );
      // initialize residual
      /* R.segment< 9 >( 0 ) = -mV9t( fastorTensorFromDoubleTensor< T >( FeTrial
       * ).data() ); */
      R( 9 ) = -alphaPTrial;

      const Tensor33t< T > Fe( X.segment( 0, 9 ).data() );

      const T dLambda = X( 10 );
      const T alphaP  = X( 9 );

      const T betaP = computeBetaPOnly( alphaP );

      // compute mandel stress
      Tensor33t< T > mandelStress = computeMandelStressOnly( Fe );

      T              f, df_dBetaP;
      Tensor33t< T > df_dMandel;
      std::tie( f, df_dMandel, df_dBetaP ) = yieldFunctionFromStressFirstOrderDerived( mandelStress, betaP );

      const Tensor33t< T > dGp = multiplyFastorTensorWithScalar( df_dMandel, dLambda );
      const Tensor33t< T > dFp = ContinuumMechanics::FiniteStrain::Plasticity::FlowIntegration::exponentialMap( dGp );

      VectorXt< T > aux = mV9t( Tensor33t< T >( einsum< iJ, JK >( Fe, dFp ) ).data() ) -
                          mV9t( fastorTensorFromDoubleTensor< T >( FeTrial ).data() );
      R( 0 )    = aux( 0 );
      R( 1 )    = aux( 1 );
      R( 2 )    = aux( 2 );
      R( 3 )    = aux( 3 );
      R( 4 )    = aux( 4 );
      R( 5 )    = aux( 5 );
      R( 6 )    = aux( 6 );
      R( 7 )    = aux( 7 );
      R( 8 )    = aux( 8 );
      R( idxA ) = ( alphaP + dLambda * df_dBetaP ) - alphaPTrial;
      R( idxF ) = f;

      return R;
    }

    /** \brief Residual and Jacobian for the full return mapping (double precision).
     *  \param X           State vector \c [Fe(9), alphaP, dLambda]^T.
     *  \param FeTrial     Trial elastic deformation gradient.
     *  \param alphaPTrial Trial accumulated plastic strain.
     */

    std::tuple< Eigen::VectorXd, Eigen::MatrixXd > computeResidualVectorAndTangent( const Eigen::VectorXd& X,
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
      std::tie( dFp, ddFp_ddGp ) = ContinuumMechanics::FiniteStrain::Plasticity::FlowIntegration::FirstOrderDerived::
        exponentialMap( dGp );

      Tensor3333d ddGp_dFe      = dLambda * einsum< ijmn, mnkL >( d2f_dMandel_dMandel, dMandel_dFe );
      Tensor33d   ddFp_ddLambda = einsum< IJKL, KL >( ddFp_ddGp, df_dMandel );

      Tensor3333d ddFp_dFe = einsum< iImn, mnkL >( ddFp_ddGp, ddGp_dFe );
      /* Tensor3333d dFeTrial_dFe = einsum< iIkL, iJ, to_IJkL >( ddFp_dFe, Fe ) +
       */
      /*                            einsum< iI, ik, JL, to_IJkL >( dFp,
       * Spatial3D::I, Spatial3D::I ); */
      Tensor3333d dFeTrial_dFe = einsum< iI, IJKL >( Fe, ddFp_dFe ) +
                                 einsum< IK, JL, to_IJKL >( Spatial3D::I, transpose( Spatial3D::I % dFp ) );

      /* std::cout << "ddFeTrial_dFe: " << dFeTrial_dFe << std::endl; */

      /* Tensor33d dFe_ddLambda = einsum< iI, iJ >( ddFp_ddLambda, Fe ); */
      Tensor33d dFe_ddLambda = einsum< Ii, iJ >( Fe, ddFp_ddLambda );

      df_dBetaP = -Constants::sqrt2_3 / fy;

      /* Tensor33d df_dFe                 = einsum< IJ, IJKL >( df_dMandel,
       * dMandel_dFe ); */
      Tensor33d df_dFe                 = einsum< IJ, IJKL >( df_dMandel, dMandel_dFe );
      std::tie( f, df_dFe, df_dBetaP ) = yieldFunction( Fe, betaP );
      df_dFe                           = einsum< IJ, IJKL >( df_dMandel, dMandel_dFe );

      /* R.segment< 9 >( 0 ) += mV9d( Tensor33d( einsum< iI, iJ >( dFp, Fe )
       * ).data() ); */
      R.segment< 9 >( 0 ) += mV9d( Tensor33d( einsum< iJ, JK >( Fe, dFp ) ).data() );
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
