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
   * \brief Finite-strain J2 plasticity with isotropic hardening. 
   *
   * Elastic law: hyperelastic, compressible Neo-Hookean (Pence–Gou potential, variant B).
   *
   *
   * \par Material parameters
   * - \b #K  — bulk modulus
   * - \b #G   — shear modulus
   * - \b #fy   — initial yield stress (reference)
   * - \b #fyInf  — saturated (asymptotic) yield stress 
   * - \b #eta — saturation rate parameter
   * - \b #H    — linear hardening modulus
   * - \b #implementationType  — algorithm selector (see below)
   * - \b #density (optional) — density
   *
   * \par State variables
   * - \b Fp  — plastic deformation gradient
   * - \b alphaP  — strain-like hardening variable
   *
   * \par Implementation variants (implementationType)
   * - \b 0: scalar return mapping (not implemented)
   * - \b 1: Full return mapping; derivatives computed analytically. 
   * - \b 2: FDAF — 
    Full return mapping; derivatives approximated via forward finite differences.
   * - \b 3: FDAC — Full return mapping; derivatives approximated via central finite differences.
   * - \b 4: CSDA — Full return mapping; derivatives via complex-step differentiation.
   * \ingroup materials_plasticity
   */

  class FiniteStrainJ2Plasticity : public MarmotMaterialFiniteStrain {
  public:
    using MarmotMaterialFiniteStrain::MarmotMaterialFiniteStrain;

    // elastic constants
    /** Bulk modulus (read from `materialProperties[0]`) */
    const double K;
    /** Shear modulus (read from `materialProperties[1]`) */ 
    const double G;
  
    // plasticity parameters
    /** Initial yield stress (read from `materialProperties[2]`) */
    const double fy;
    /** Saturated (asymptotic) yield stress (read from `materialProperties[3]`) */
    const double fyInf;
    /** Saturation rate parameter (read from `materialProperties[4]`) */
    const double eta; 
    /** Linear hardening modulus (read from `materialProperties[5]`) */
    const double H; 
    
    // implementation
    /** Algorithm variant selector (read from `materialProperties[6]`). */
    const int implementationType; 
    // mass properties;
/** Density (read from `materialProperties[7]`) (if provided). */
    const double density; 

    /**
     * \brief Construct the finite-strain J2 plasticity model.
     * \param materialProperties Array with parameters: #K, #G, #fy, #fyInf, #eta, #H, #implementationType, #density (optional).
     * \param nMaterialProperties Length of `materialProperties`.
     * \param materialLabel Material label.
     */
    FiniteStrainJ2Plasticity( const double* materialProperties, int nMaterialProperties, int materialLabel );

    /**
     * \brief Compute the Kirchhoff stress and the algorithmic tangent for the current step.  Performs an elastic trial; if yielding occurs, return mapping is performed. 
     *
     * @param[in,out] response  
     *   - `tau` - Kirchhoff stress tensor \f$\boldsymbol{\tau}\f$.
     *   - `elasticEnergyDensity` - elastic energy density  \f$\psi\f$.
     *   - `rho` - density (set to 1.0 here).  
     * @param[in,out] tangents
     *   - `dTau_dF` - algorithmic tangent \f$\partial\boldsymbol{\tau}/\partial\mathbf{F}\f$.
     * @param[in]  deformation
     *   - `F` - deformation gradient \f$\mathbf{F}\f$.  
     * @param[in]  timeIncrement  
     *   - `t` - old (pseudo-)time.
     *   - `dT`- (pseudo-)time increment.
     *   
     * Template parameter `<3>` indicates 3D.
     */
    void computeStress( ConstitutiveResponse< 3 >& response,
                        AlgorithmicModuli< 3 >&    tangents,
                        const Deformation< 3 >&    deformation,
                        const TimeIncrement&       timeIncrement );

    /** \brief Scalar-return mapping variant; would solve a 1D consistency equation for Δλ (not yet implemented). */

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

    /** \brief Get the number of required state variables (9 for `Fp` and 1 for `alphaP`). 
     */

    int getNumberOfRequiredStateVars() { return FiniteStrainJ2PlasticityStateVarManager::layout.nRequiredStateVars; }

    /** \brief Return the material density if provided in material parameters */

    double getDensity() { return density; }

 
/** \brief State variable manager for `Fp` and `alphaP`; provides named views and layout. */
    class FiniteStrainJ2PlasticityStateVarManager : public MarmotStateVarVectorManager {

    public:
/** \brief Memory layout of internal variables: 9 for \c Fp and 1 for \c alphaP (total 10). */
      inline const static auto layout = makeLayout( {
        { .name = "Fp", .length = 9 },
        { .name = "alphaP", .length = 1 },
      } );
/** \brief Plastic deformation gradient \f$\mathbf F_p\f$. */
      Fastor::TensorMap< double, 3, 3 > Fp;
/** \brief Strain-like hardening variable \f$\alpha_p\f$. */
      double&                           alphaP;

/** \brief Bind the manager to an external contiguous buffer holding (\c Fp,\c alphaP). */
      FiniteStrainJ2PlasticityStateVarManager( double* theStateVarVector )
        : MarmotStateVarVectorManager( theStateVarVector, layout ), Fp( &find( "Fp" ) ), alphaP( find( "alphaP" ) ) {};
    };
    std::unique_ptr< FiniteStrainJ2PlasticityStateVarManager > stateVars;

/** \brief Bind external storage for internal variables (`Fp`, `alphaP`).
 *  @param stateVars  Pointer to a contiguous array provided by the caller.
 *  @param nStateVars Number of entries in that array.
 *  @throws std::invalid_argument If nStateVars < getNumberOfRequiredStateVars().
 *
 *  Creates/updates the internal state manager that maps (Fp, alphaP) onto this buffer.
 */    void assignStateVars( double* stateVars, int nStateVars );

    /** \brief Expose a named view into the state vector.
     *  \param result One of: `Fp` (length 9), `alphaP` (length 1).
     */
    StateView getStateView( const std::string& result );

    /** \brief Initialize state (sets \f$\mathbf F_p = \mathbf I\f$) */
    void initializeYourself();

/** \brief Yield function \f$f(\mathbf F_e,\beta_p)\f$ and first derivatives.
 *
 *  Definition:
 *  \f[
 *    f(\mathbf M,\beta_p)=\frac{1}{f_y}\Big(\rho(\mathbf M)-\sqrt{\tfrac{2}{3}}\,\beta_p\Big),
 *    \qquad \rho(\mathbf M)=\sqrt{\operatorname{dev}\mathbf M:\operatorname{dev}\mathbf M},
 *  \f]
 *  with Mandel stress \f$\mathbf M=\mathbf C_e:\mathbf{PK2}(\mathbf C_e)\f$, \f$\mathbf C_e=\mathbf F_e^\mathrm T\mathbf F_e\f$.
 *
 *  Derivatives used:
 *  \f[
 *    \frac{\partial f}{\partial \mathbf M}=\frac{1}{f_y}\,\frac{\operatorname{dev}\mathbf M}{\rho},\qquad
 *    \frac{\partial f}{\partial \beta_p}=-\frac{\sqrt{2/3}}{f_y},\qquad
 *    \frac{\partial f}{\partial \mathbf F_e}
 *     = \frac{\partial f}{\partial \mathbf M}:\frac{\partial \mathbf M}{\partial \mathbf F_e}.
 *  \f]
 *
 *  @param Fe     Elastic deformation gradient \f$\mathbf F_e\f$.
 *  @param betaP  Stress-like hardening variable \f$\beta_p\f$.
 *  @return \f$\{\,f,\;\partial f/\partial \mathbf F_e,\;\partial f/\partial \beta_p\,\}\f$.
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

/** \brief Yield function \f$f(\mathbf M,\beta_p)\f$ w.r.t. Mandel stress (templated).
 *  @tparam T Scalar type (double, autodiff dual, complex).
 *  @param mandelStress Mandel stress \f$\mathbf M\f$.
 *  @param betaP        Stress-like hardening variable \f$\beta_p\f$.
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


/** \brief Yield function w.r.t. Mandel stress (double version).
 *  @param mandelStress Mandel stress \f$\mathbf M\f$.
 *  @param betaP        Stress-like hardening variable \f$\beta_p\f$.
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

/** \brief First-order derivatives of \f$f(\mathbf M,\beta_p)\f$ w.r.t. \f$\mathbf M \f$.
 *  @tparam T Scalar type (double, autodiff dual, complex).
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

    /** \brief Check yield condition \f$f(\mathbf F_e,\beta_p)>0\f$.
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

   
/** \brief Mandel stress \f$\mathbf M=\mathbf C_e:\mathbf S(\mathbf C_e)\f$ and \f$\partial \mathbf M/\partial \mathbf F_e\f$. */

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

 
/** \brief Isotropic hardening law \f$\beta(\alpha_p)=f_{y\infty}+(f_y-f_{y\infty})e^{-\eta\alpha_p}+H\,\alpha_p\f$ (templated) */

    template < typename T >
    T computeBetaPOnly( const T alphaP )
    {
      const T beta = fyInf + ( fy - fyInf ) * exp( -alphaP * eta ) + alphaP * H;
      return beta;
    }

    /**  \brief Isotropic hardening law
     * \f$\beta(\alpha_p)=f_{y\infty}+(f_y-f_{y\infty})e^{-\eta\alpha_p}+H\,\alpha_p\f$. 
     * And its derivative \f$\partial\beta/\partial\alpha_p\f$.
     */

    std::tuple< double, double > computeBetaP( const double alphaP )
    {
      const double beta           = fyInf + ( fy - fyInf ) * exp( -alphaP * eta ) + alphaP * H;
      const double dBetaP_dAlphaP = -eta * ( fy - fyInf ) * exp( -alphaP * eta ) + H;
      return { beta, dBetaP_dAlphaP };
    }


/** \brief Scalar consistency function \f$g(\Delta\lambda)\f$ and its derivative for the scalar-return mapping. */

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


/** \brief Flow direction for return mapping, \f$\partial f/\partial \mathbf S\f$ (templated).
 *  @tparam T   Scalar type (double, AD, complex).
 *  @param Fe   Elastic deformation gradient \f$\mathbf F_e\f$.
 *  @param betaP Stress-like hardening \f$\beta_p\f$.
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
 *  @tparam T    Scalar type (double, autodiff dual, complex).
 *  @param df_dS Flow direction in Mandel-stress space.
 *  @param dLambda Plastic multiplier increment \f$\Delta\lambda\f$.
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


/** \brief Residual vector \f$\mathbf R(\mathbf X)\f$ for the full return mapping (templated).
 *  \details State vector \f$\mathbf X\f$ has 11 unknowns: 9 for \f$\mathbf F_e\f$, 1 for \f$\alpha_p\f$, 1 for \f$\Delta\lambda\f$.
 *           Residual components:
 *            - \f$R_{0..8}\f$: elastic consistency (updated \f$\mathbf F_e\f$ vs trial),
 *            - \f$R_{9}\f$: \f$\alpha_p\f$ update,
 *            - \f$R_{10}\f$: yield \f$f\f$.
 *  @tparam T         Scalar type (double/dual/complex).
 *  @param X          Current iterate \f$[\mathrm{vec}(\mathbf F_e),\;\alpha_p,\;\Delta\lambda]^\mathrm T\f$.
 *  @param FeTrial    Trial elastic deformation gradient \f$\mathbf F_e^\text{trial}\f$.
 *  @param alphaPTrial Trial strain-like hardening variable \f$\alpha_p^\text{trial}\f$.
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
 *  @param X           \f$[\mathrm{vec}(\mathbf F_e),\;\alpha_p,\;\Delta\lambda]^\mathrm T\f$.
 *  @param FeTrial     Trial elastic gradient \f$\mathbf F_e^\text{trial}\f$.
 *  @param alphaPTrial Trial strain-like hardening variable \f$\alpha_p^\text{trial}\f$.
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
