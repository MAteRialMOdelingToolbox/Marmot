#include "Marmot/FiniteStrainJ2Plasticity.h"
#include "Marmot/MarmotDeformationMeasures.h"
#include "Marmot/MarmotEnergyDensityFunctions.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include "Marmot/MarmotNumericalDifferentiation.h"
#include "Marmot/MarmotNumericalDifferentiationForFastor.h"
#include "Marmot/MarmotStressMeasures.h"
#include <autodiff/forward/dual/dual.hpp>

namespace Marmot::Materials {

  using namespace Marmot;
  using namespace Fastor;
  using namespace FastorIndices;
  using namespace FastorStandardTensors;

  FiniteStrainJ2Plasticity::FiniteStrainJ2Plasticity( const double* materialProperties,
                                                      int           nMaterialProperties,
                                                      int           materialLabel )
    : MarmotMaterialFiniteStrain( materialProperties, nMaterialProperties, materialLabel ),
      K( materialProperties[0] ),
      G( materialProperties[1] ),
      fy( materialProperties[2] ),
      fyInf( materialProperties[3] ),
      eta( materialProperties[4] ),
      H( materialProperties[5] ),
      implementationType( materialProperties[6] ),
      density( nMaterialProperties > 7 ? materialProperties[7] : 0.0 ) // TODO: make mandatory material parameter
  {
    initializeStateLayout();
  }

  void FiniteStrainJ2Plasticity::computeStress( ConstitutiveResponse< 3 >& response,
                                                AlgorithmicModuli< 3 >&    tangents,
                                                const Deformation< 3 >&    deformation,
                                                const TimeIncrement&       timeIncrement )
  {
    switch ( implementationType ) {

    case 0: computeStressWithScalarReturnMapping( response, tangents, deformation, timeIncrement ); break;

    case 1: computeStressWithFullReturnMapping( response, tangents, deformation, timeIncrement ); break;

    case 2: computeStressFDAF( response, tangents, deformation, timeIncrement ); break;

    case 3: computeStressFDAC( response, tangents, deformation, timeIncrement ); break;

    case 4: computeStressCSDA( response, tangents, deformation, timeIncrement ); break;
    default: throw std::invalid_argument( "implementation type not supported" );
    };
  }
  void FiniteStrainJ2Plasticity::computeStressWithScalarReturnMapping( ConstitutiveResponse< 3 >& response,
                                                                       AlgorithmicModuli< 3 >&    tangents,
                                                                       const Deformation< 3 >&    deformation,
                                                                       const TimeIncrement&       timeIncrement )
  {
    throw std::invalid_argument( "not implemented yet" );
  }
  void FiniteStrainJ2Plasticity::computeStressWithFullReturnMapping( ConstitutiveResponse< 3 >& response,
                                                                     AlgorithmicModuli< 3 >&    tangents,
                                                                     const Deformation< 3 >&    deformation,
                                                                     const TimeIncrement&       timeIncrement )
  {

    Tensor33d&      Fp = stateLayout.getAs< Tensor33d& >( response.stateVars, "Fp" );
    const Tensor33d FpOld( Fp );
    double&         alphaP    = stateLayout.getAs< double& >( response.stateVars, "alphaP" );
    const double    alphaPOld = alphaP;

    using namespace Marmot;
    using namespace Fastor;
    using namespace Eigen;
    using namespace autodiff;
    using namespace FastorIndices;
    using namespace FastorStandardTensors;

    Tensor33d FeTrial = deformation.F % Fastor::inverse( FpOld );
    double    betaP, dBetaP_dAlphaP;
    std::tie( betaP, dBetaP_dAlphaP ) = computeBetaP( alphaPOld );

    Tensor33d dFp;
    dFp.eye();
    Tensor33d Fe = FeTrial;
    /* std::cout << "FeTrial: " << std::endl << FeTrial << std::endl; */
    if ( isYielding( FeTrial, betaP ) ) {

      size_t counter = 0;

      using mV9d = Eigen::Map< Eigen::Matrix< double, 9, 1 > >;
      VectorXd X( 11 );
      X.segment( 0, 9 )  = mV9d( FeTrial.data() );
      X( 9 )             = alphaPOld;
      X( 10 )            = 0;
      VectorXd        dX = VectorXd::Zero( 11 );
      VectorXd        R  = VectorXd::Zero( 11 );
      Eigen::MatrixXd dR_dX( 11, 11 );

      std::tie( R, dR_dX ) = computeResidualVectorAndTangent( X, FeTrial, alphaPOld );

      while ( R.norm() > 1e-12 || dX.norm() > 1e-12 ) {

        if ( counter > 10 )
          throw std::runtime_error( "inner newton not converged" );

        dX = -dR_dX.colPivHouseholderQr().solve( R );
        X += dX;
        std::tie( R, dR_dX ) = computeResidualVectorAndTangent( X, FeTrial, alphaPOld );
        counter += 1;
      }
      /* std::cout << "inner newton iters: " << counter << std::endl; */

      // update plastic deformation increment
      Fe              = X.segment( 0, 9 ).data();
      dFp             = Fastor::inverse( Fe ) % FeTrial;
      alphaP          = X( 9 );
      Tensor33d FpNew = dFp % Fp;
      memcpy( Fp.data(), FpNew.data(), 9 * sizeof( double ) );

      using namespace ContinuumMechanics;
      double      psi_;
      Tensor33d   Ce, dPsi_dCe;
      Tensor3333d dCe_dFe, d2Psi_dCedCe;
      std::tie( Ce, dCe_dFe ) = DeformationMeasures::FirstOrderDerived::rightCauchyGreen( Fe );

      // compute energy density, first and second partial derivatives wrt Cauchy
      // Green deformation
      std::tie( psi_, dPsi_dCe, d2Psi_dCedCe ) = EnergyDensityFunctions::SecondOrderDerived::PenceGouPotentialB( Ce,
                                                                                                                 K,
                                                                                                                 G );
      // compute Kirchhoff stress
      Tensor33d   PK2 = 2. * dPsi_dCe;
      Tensor3333d dTau_dPK2, dTau_dFe_partial;
      std::tie( response.tau,
                dTau_dPK2,
                dTau_dFe_partial )  = StressMeasures::FirstOrderDerived::KirchhoffStressFromPK2( PK2, Fe );
      response.rho                  = 1.0;
      response.elasticEnergyDensity = psi_;

      // compute tangent operator
      using mM9d = Eigen::Map< Eigen::Matrix< double, 9, 9 > >;

      MatrixXd dYdDeformation              = MatrixXd::Zero( 11, 11 );
      dYdDeformation.block< 9, 9 >( 0, 0 ) = mM9d( Tensor3333d( einsum< IK, JL, to_IJKL >( Spatial3D::I,
                                                                                           transpose( Fastor::inverse(
                                                                                             FpOld ) ) ) )
                                                     .data() )
                                               .transpose();
      MatrixXd dXdDeformation = dR_dX.colPivHouseholderQr().solve( dYdDeformation );

      Tensor3333d dFe_dF = Tensor3333d( Matrix9d( dXdDeformation.block< 9, 9 >( 0, 0 ).transpose() ).data() );

      Tensor3333d dPK2_dFe = einsum< ijKL, KLMN >( 2. * d2Psi_dCedCe, dCe_dFe );
      Tensor3333d dPK2_dF  = einsum< ijKL, KLMN >( dPK2_dFe, dFe_dF );

      /* tangents.dTau_dF = einsum< IJKL, KLMN >( dTau_dPK2, dPK2_dF ) +
       * dTau_dF_partial; */
      tangents.dTau_dF = einsum< IJKL, KLMN >( dTau_dPK2, dPK2_dF ) + einsum< ijKL, KLMN >( dTau_dFe_partial, dFe_dF );
    }
    else {
      using namespace Marmot::ContinuumMechanics;
      double      psi_;
      Tensor33d   Ce, dPsi_dCe;
      Tensor3333d dCe_dFe, d2Psi_dCedCe;
      std::tie( Ce, dCe_dFe ) = DeformationMeasures::FirstOrderDerived::rightCauchyGreen( Fe );

      // compute energy density, first and second partial derivatives wrt Cauchy
      // Green deformation
      std::tie( psi_, dPsi_dCe, d2Psi_dCedCe ) = EnergyDensityFunctions::SecondOrderDerived::PenceGouPotentialB( Ce,
                                                                                                                 K,
                                                                                                                 G );
      // compute Kirchhoff stress
      Tensor33d   PK2 = 2. * dPsi_dCe;
      Tensor3333d dTau_dPK2, dTau_dFe_partial;
      std::tie( response.tau,
                dTau_dPK2,
                dTau_dFe_partial )  = StressMeasures::FirstOrderDerived::KirchhoffStressFromPK2( PK2, Fe );
      response.rho                  = 1.0;
      response.elasticEnergyDensity = psi_;

      // compute tangent operator
      Tensor3333d dPK2_dFe = einsum< ijKL, KLMN >( 2. * d2Psi_dCedCe, dCe_dFe );
      Tensor3333d dFe_dF   = einsum< IK, JL, to_IJKL >( Spatial3D::I, transpose( Fastor::inverse( FpOld ) ) );
      Tensor3333d dPK2_dF  = einsum< ijKL, KLMN >( dPK2_dFe, dFe_dF );

      /* tangents.dTau_dF = einsum< IJKL, KLMN >( dTau_dPK2, dPK2_dF ) +
       * dTau_dF_partial; */
      tangents.dTau_dF = einsum< IJKL, KLMN >( dTau_dPK2, dPK2_dF ) + einsum< ijKL, KLMN >( dTau_dFe_partial, dFe_dF );
    }
  }

  void FiniteStrainJ2Plasticity::computeStressFDAF( ConstitutiveResponse< 3 >& response,
                                                    AlgorithmicModuli< 3 >&    tangents,
                                                    const Deformation< 3 >&    deformation,
                                                    const TimeIncrement&       timeIncrement )
  {

    Tensor33d&      Fp = stateLayout.getAs< Tensor33d& >( response.stateVars, "Fp" );
    const Tensor33d FpOld( Fp );
    double&         alphaP    = stateLayout.getAs< double& >( response.stateVars, "alphaP" );
    const double    alphaPOld = alphaP;

    using namespace Marmot;
    using namespace Fastor;
    using namespace Eigen;
    using namespace autodiff;
    using namespace FastorIndices;
    using namespace FastorStandardTensors;

    Tensor33d FeTrial = deformation.F % Fastor::inverse( FpOld );
    double    betaP, dBetaP_dAlphaP;
    std::tie( betaP, dBetaP_dAlphaP ) = computeBetaP( alphaPOld );

    Tensor33d dFp;
    dFp.eye();
    Tensor33d Fe = FeTrial;
    /* std::cout << "FeTrial: " << std::endl << FeTrial << std::endl; */
    if ( isYielding( FeTrial, betaP ) ) {

      size_t counter = 0;

      using mV9d = Eigen::Map< Eigen::Matrix< double, 9, 1 > >;
      VectorXd X( 11 );
      X.segment( 0, 9 )  = mV9d( FeTrial.data() );
      X( 9 )             = alphaPOld;
      X( 10 )            = 0;
      VectorXd        dX = VectorXd::Zero( 11 );
      VectorXd        R  = VectorXd::Zero( 11 );
      Eigen::MatrixXd dR_dX( 11, 11 );

      /* std::tie( R, dR_dX ) = computeResidualVectorAndTangent( X, FeTrial,
       * alphaPOld ); */
      R     = computeResidualVector( X, FeTrial, alphaP );
      dR_dX = NumericalAlgorithms::Differentiation::forwardDifference(
        [&]( const VectorXd& X_ ) {
          VectorXd Xout = computeResidualVector( X_, FeTrial, alphaP );
          return Xout;
        },
        X );
      try {
        while ( R.norm() > 1e-12 || dX.norm() > 1e-12 ) {

          if ( counter > 10 )
            throw std::runtime_error( "inner newton not converged" );

          dX = -dR_dX.colPivHouseholderQr().solve( R );
          X += dX;
          /* std::tie( R, dR_dX ) = computeResidualVectorAndTangent( X, FeTrial,
           * alphaPOld ); */
          R     = computeResidualVector( X, FeTrial, alphaP );
          dR_dX = NumericalAlgorithms::Differentiation::forwardDifference(
            [&]( const VectorXd& X_ ) {
              VectorXd Xout = computeResidualVector( X_, FeTrial, alphaP );
              return Xout;
            },
            X );
          counter += 1;
        }
      }
      catch ( std::exception& e ) {
        throw std::runtime_error( "return mapping failed: " + std::string( e.what() ) );
      }
      /* std::cout << "inner newton iters: " << counter << std::endl; */

      // update plastic deformation increment
      Fe              = X.segment( 0, 9 ).data();
      dFp             = Fastor::inverse( Fe ) % FeTrial;
      alphaP          = X( 9 );
      Tensor33d FpNew = dFp % Fp;
      memcpy( Fp.data(), FpNew.data(), 9 * sizeof( double ) );

      using namespace ContinuumMechanics;
      double      psi_;
      Tensor33d   Ce, dPsi_dCe;
      Tensor3333d dCe_dFe, d2Psi_dCedCe;
      std::tie( Ce, dCe_dFe ) = DeformationMeasures::FirstOrderDerived::rightCauchyGreen( Fe );

      // compute energy density, first and second partial derivatives wrt Cauchy
      // Green deformation
      using func_type    = std::function< Tensor33d( const Tensor33d& ) >;
      func_type computeS = [&]( const Tensor33d& Ce_ ) {
        const auto [_psi, _dPsi_dCe] = EnergyDensityFunctions::FirstOrderDerived::PenceGouPotentialB( Ce_, K, G );
        return _dPsi_dCe;
      };

      d2Psi_dCedCe = NumericalAlgorithms::Differentiation::TensorToTensor::forwardDifference( computeS, Ce );
      std::tie( psi_, dPsi_dCe ) = EnergyDensityFunctions::FirstOrderDerived::PenceGouPotentialB( Ce, K, G );
      // compute Kirchhoff stress
      Tensor33d   PK2 = 2. * dPsi_dCe;
      Tensor3333d dTau_dPK2, dTau_dFe_partial;
      std::tie( response.tau,
                dTau_dPK2,
                dTau_dFe_partial )  = StressMeasures::FirstOrderDerived::KirchhoffStressFromPK2( PK2, Fe );
      response.rho                  = 1.0;
      response.elasticEnergyDensity = psi_;

      // compute tangent operator
      using mM9d = Eigen::Map< Eigen::Matrix< double, 9, 9 > >;

      MatrixXd dYdDeformation              = MatrixXd::Zero( 11, 11 );
      dYdDeformation.block< 9, 9 >( 0, 0 ) = mM9d( Tensor3333d( einsum< IK, JL, to_IJKL >( Spatial3D::I,
                                                                                           transpose( Fastor::inverse(
                                                                                             FpOld ) ) ) )
                                                     .data() )
                                               .transpose();
      MatrixXd dXdDeformation = dR_dX.colPivHouseholderQr().solve( dYdDeformation );

      Tensor3333d dFe_dF = Tensor3333d( Matrix9d( dXdDeformation.block< 9, 9 >( 0, 0 ).transpose() ).data() );

      Tensor3333d dPK2_dFe = einsum< ijKL, KLMN >( 2. * d2Psi_dCedCe, dCe_dFe );
      Tensor3333d dPK2_dF  = einsum< ijKL, KLMN >( dPK2_dFe, dFe_dF );

      /* tangents.dTau_dF = einsum< IJKL, KLMN >( dTau_dPK2, dPK2_dF ) +
       * dTau_dF_partial; */
      tangents.dTau_dF = einsum< IJKL, KLMN >( dTau_dPK2, dPK2_dF ) + einsum< ijKL, KLMN >( dTau_dFe_partial, dFe_dF );
    }
    else {
      using namespace Marmot::ContinuumMechanics;
      double      psi_;
      Tensor33d   Ce, dPsi_dCe;
      Tensor3333d dCe_dFe, d2Psi_dCedCe;
      std::tie( Ce, dCe_dFe ) = DeformationMeasures::FirstOrderDerived::rightCauchyGreen( Fe );

      // compute energy density, first and second partial derivatives wrt Cauchy
      // Green deformation
      using func_type    = std::function< Tensor33d( const Tensor33d& ) >;
      func_type computeS = [&]( const Tensor33d& Ce_ ) {
        const auto [_psi, _dPsi_dCe] = EnergyDensityFunctions::FirstOrderDerived::PenceGouPotentialB( Ce_, K, G );
        return _dPsi_dCe;
      };

      d2Psi_dCedCe = NumericalAlgorithms::Differentiation::TensorToTensor::forwardDifference( computeS, Ce );

      std::tie( psi_, dPsi_dCe ) = EnergyDensityFunctions::FirstOrderDerived::PenceGouPotentialB( Ce, K, G );
      // compute Kirchhoff stress
      Tensor33d   PK2 = 2. * dPsi_dCe;
      Tensor3333d dTau_dPK2, dTau_dFe_partial;
      std::tie( response.tau,
                dTau_dPK2,
                dTau_dFe_partial )  = StressMeasures::FirstOrderDerived::KirchhoffStressFromPK2( PK2, Fe );
      response.rho                  = 1.0;
      response.elasticEnergyDensity = psi_;

      // compute tangent operator
      Tensor3333d dPK2_dFe = einsum< ijKL, KLMN >( 2. * d2Psi_dCedCe, dCe_dFe );
      Tensor3333d dFe_dF   = einsum< IK, JL, to_IJKL >( Spatial3D::I, transpose( Fastor::inverse( FpOld ) ) );
      Tensor3333d dPK2_dF  = einsum< ijKL, KLMN >( dPK2_dFe, dFe_dF );

      /* tangents.dTau_dF = einsum< IJKL, KLMN >( dTau_dPK2, dPK2_dF ) +
       * dTau_dF_partial; */
      tangents.dTau_dF = einsum< IJKL, KLMN >( dTau_dPK2, dPK2_dF ) + einsum< ijKL, KLMN >( dTau_dFe_partial, dFe_dF );
    }
  }
  void FiniteStrainJ2Plasticity::computeStressFDAC( ConstitutiveResponse< 3 >& response,
                                                    AlgorithmicModuli< 3 >&    tangents,
                                                    const Deformation< 3 >&    deformation,
                                                    const TimeIncrement&       timeIncrement )
  {

    Tensor33d&      Fp = stateLayout.getAs< Tensor33d& >( response.stateVars, "Fp" );
    const Tensor33d FpOld( Fp );
    double&         alphaP    = stateLayout.getAs< double& >( response.stateVars, "alphaP" );
    const double    alphaPOld = alphaP;

    using namespace Marmot;
    using namespace Fastor;
    using namespace Eigen;
    using namespace autodiff;
    using namespace FastorIndices;
    using namespace FastorStandardTensors;

    Tensor33d FeTrial = deformation.F % Fastor::inverse( FpOld );
    double    betaP, dBetaP_dAlphaP;
    std::tie( betaP, dBetaP_dAlphaP ) = computeBetaP( alphaPOld );

    Tensor33d dFp;
    dFp.eye();
    Tensor33d Fe = FeTrial;
    /* std::cout << "FeTrial: " << std::endl << FeTrial << std::endl; */
    if ( isYielding( FeTrial, betaP ) ) {

      size_t counter = 0;

      using mV9d = Eigen::Map< Eigen::Matrix< double, 9, 1 > >;
      VectorXd X( 11 );
      X.segment( 0, 9 )  = mV9d( FeTrial.data() );
      X( 9 )             = alphaPOld;
      X( 10 )            = 0;
      VectorXd        dX = VectorXd::Zero( 11 );
      VectorXd        R  = VectorXd::Zero( 11 );
      Eigen::MatrixXd dR_dX( 11, 11 );

      /* std::tie( R, dR_dX ) = computeResidualVectorAndTangent( X, FeTrial,
       * alphaPOld ); */
      R     = computeResidualVector( X, FeTrial, alphaP );
      dR_dX = NumericalAlgorithms::Differentiation::centralDifference(
        [&]( const VectorXd& X_ ) {
          VectorXd Xout = computeResidualVector( X_, FeTrial, alphaP );
          return Xout;
        },
        X );
      try {
        while ( R.norm() > 1e-12 || dX.norm() > 1e-12 ) {

          if ( counter > 10 )
            throw std::runtime_error( "inner newton not converged" );

          dX = -dR_dX.colPivHouseholderQr().solve( R );
          X += dX;
          /* std::tie( R, dR_dX ) = computeResidualVectorAndTangent( X, FeTrial,
           * alphaPOld ); */
          R     = computeResidualVector( X, FeTrial, alphaP );
          dR_dX = NumericalAlgorithms::Differentiation::centralDifference(
            [&]( const VectorXd& X_ ) {
              VectorXd Xout = computeResidualVector( X_, FeTrial, alphaP );
              return Xout;
            },
            X );
          counter += 1;
        }
      }
      catch ( std::exception& e ) {
        throw std::runtime_error( "return mapping failed: " + std::string( e.what() ) );
      }
      /* std::cout << "inner newton iters: " << counter << std::endl; */

      // update plastic deformation increment
      Fe              = X.segment( 0, 9 ).data();
      dFp             = Fastor::inverse( Fe ) % FeTrial;
      alphaP          = X( 9 );
      Tensor33d FpNew = dFp % Fp;
      memcpy( Fp.data(), FpNew.data(), 9 * sizeof( double ) );

      using namespace ContinuumMechanics;
      double      psi_;
      Tensor33d   Ce, dPsi_dCe;
      Tensor3333d dCe_dFe, d2Psi_dCedCe;
      std::tie( Ce, dCe_dFe ) = DeformationMeasures::FirstOrderDerived::rightCauchyGreen( Fe );

      // compute energy density, first and second partial derivatives wrt Cauchy
      // Green deformation
      using func_type    = std::function< Tensor33d( const Tensor33d& ) >;
      func_type computeS = [&]( const Tensor33d& Ce_ ) {
        const auto [_psi, _dPsi_dCe] = EnergyDensityFunctions::FirstOrderDerived::PenceGouPotentialB( Ce_, K, G );
        return _dPsi_dCe;
      };

      d2Psi_dCedCe = NumericalAlgorithms::Differentiation::TensorToTensor::centralDifference( computeS, Ce );
      std::tie( psi_, dPsi_dCe ) = EnergyDensityFunctions::FirstOrderDerived::PenceGouPotentialB( Ce, K, G );
      // compute Kirchhoff stress
      Tensor33d   PK2 = 2. * dPsi_dCe;
      Tensor3333d dTau_dPK2, dTau_dFe_partial;
      std::tie( response.tau,
                dTau_dPK2,
                dTau_dFe_partial )  = StressMeasures::FirstOrderDerived::KirchhoffStressFromPK2( PK2, Fe );
      response.rho                  = 1.0;
      response.elasticEnergyDensity = psi_;

      // compute tangent operator
      using mM9d = Eigen::Map< Eigen::Matrix< double, 9, 9 > >;

      MatrixXd dYdDeformation              = MatrixXd::Zero( 11, 11 );
      dYdDeformation.block< 9, 9 >( 0, 0 ) = mM9d( Tensor3333d( einsum< IK, JL, to_IJKL >( Spatial3D::I,
                                                                                           transpose( Fastor::inverse(
                                                                                             FpOld ) ) ) )
                                                     .data() )
                                               .transpose();
      MatrixXd dXdDeformation = dR_dX.colPivHouseholderQr().solve( dYdDeformation );

      Tensor3333d dFe_dF = Tensor3333d( Matrix9d( dXdDeformation.block< 9, 9 >( 0, 0 ).transpose() ).data() );

      Tensor3333d dPK2_dFe = einsum< ijKL, KLMN >( 2. * d2Psi_dCedCe, dCe_dFe );
      Tensor3333d dPK2_dF  = einsum< ijKL, KLMN >( dPK2_dFe, dFe_dF );

      /* tangents.dTau_dF = einsum< IJKL, KLMN >( dTau_dPK2, dPK2_dF ) +
       * dTau_dF_partial; */
      tangents.dTau_dF = einsum< IJKL, KLMN >( dTau_dPK2, dPK2_dF ) + einsum< ijKL, KLMN >( dTau_dFe_partial, dFe_dF );
    }
    else {
      using namespace Marmot::ContinuumMechanics;
      double      psi_;
      Tensor33d   Ce, dPsi_dCe;
      Tensor3333d dCe_dFe, d2Psi_dCedCe;
      std::tie( Ce, dCe_dFe ) = DeformationMeasures::FirstOrderDerived::rightCauchyGreen( Fe );

      // compute energy density, first and second partial derivatives wrt Cauchy
      // Green deformation
      using func_type    = std::function< Tensor33d( const Tensor33d& ) >;
      func_type computeS = [&]( const Tensor33d& Ce_ ) {
        const auto [_psi, _dPsi_dCe] = EnergyDensityFunctions::FirstOrderDerived::PenceGouPotentialB( Ce_, K, G );
        return _dPsi_dCe;
      };

      d2Psi_dCedCe = NumericalAlgorithms::Differentiation::TensorToTensor::centralDifference( computeS, Ce );

      std::tie( psi_, dPsi_dCe ) = EnergyDensityFunctions::FirstOrderDerived::PenceGouPotentialB( Ce, K, G );
      // compute Kirchhoff stress
      Tensor33d   PK2 = 2. * dPsi_dCe;
      Tensor3333d dTau_dPK2, dTau_dFe_partial;
      std::tie( response.tau,
                dTau_dPK2,
                dTau_dFe_partial )  = StressMeasures::FirstOrderDerived::KirchhoffStressFromPK2( PK2, Fe );
      response.rho                  = 1.0;
      response.elasticEnergyDensity = psi_;

      // compute tangent operator
      Tensor3333d dPK2_dFe = einsum< ijKL, KLMN >( 2. * d2Psi_dCedCe, dCe_dFe );
      Tensor3333d dFe_dF   = einsum< IK, JL, to_IJKL >( Spatial3D::I, transpose( Fastor::inverse( FpOld ) ) );
      Tensor3333d dPK2_dF  = einsum< ijKL, KLMN >( dPK2_dFe, dFe_dF );

      /* tangents.dTau_dF = einsum< IJKL, KLMN >( dTau_dPK2, dPK2_dF ) +
       * dTau_dF_partial; */
      tangents.dTau_dF = einsum< IJKL, KLMN >( dTau_dPK2, dPK2_dF ) + einsum< ijKL, KLMN >( dTau_dFe_partial, dFe_dF );
    }
  }

  void FiniteStrainJ2Plasticity::computeStressCSDA( ConstitutiveResponse< 3 >& response,
                                                    AlgorithmicModuli< 3 >&    tangents,
                                                    const Deformation< 3 >&    deformation,
                                                    const TimeIncrement&       timeIncrement )
  {

    Tensor33d&      Fp = stateLayout.getAs< Tensor33d& >( response.stateVars, "Fp" );
    const Tensor33d FpOld( Fp );
    double&         alphaP    = stateLayout.getAs< double& >( response.stateVars, "alphaP" );
    const double    alphaPOld = alphaP;

    using namespace Marmot;
    using namespace Fastor;
    using namespace Eigen;
    using namespace autodiff;
    using namespace FastorIndices;
    using namespace FastorStandardTensors;

    Tensor33d FeTrial = deformation.F % Fastor::inverse( FpOld );
    double    betaP, dBetaP_dAlphaP;
    std::tie( betaP, dBetaP_dAlphaP ) = computeBetaP( alphaPOld );

    Tensor33d dFp;
    dFp.eye();
    Tensor33d Fe = FeTrial;
    /* std::cout << "FeTrial: " << std::endl << FeTrial << std::endl; */
    if ( isYielding( FeTrial, betaP ) ) {

      size_t counter = 0;

      using mV9d = Eigen::Map< Eigen::Matrix< double, 9, 1 > >;
      VectorXd X( 11 );
      X.segment( 0, 9 )  = mV9d( FeTrial.data() );
      X( 9 )             = alphaPOld;
      X( 10 )            = 0;
      VectorXd        dX = VectorXd::Zero( 11 );
      VectorXd        R  = VectorXd::Zero( 11 );
      Eigen::MatrixXd dR_dX( 11, 11 );

      /* std::tie( R, dR_dX ) = computeResidualVectorAndTangent( X, FeTrial,
       * alphaPOld ); */
      /* R     = computeResidualVector( X, FeTrial, alphaP ); */
      std::tie( R, dR_dX ) = NumericalAlgorithms::Differentiation::Complex::forwardDifference(
        [&]( const VectorXcd& X_ ) {
          VectorXcd Xout = computeResidualVector( X_, FeTrial, alphaP );
          return Xout;
        },
        X );
      R = computeResidualVector( X, FeTrial, alphaP );
      try {
        while ( R.norm() > 1e-12 || dX.norm() > 1e-12 ) {

          if ( counter > 10 )
            throw std::runtime_error( "inner newton not converged" );

          dX = -dR_dX.colPivHouseholderQr().solve( R );
          X += dX;
          std::tie( R, dR_dX ) = NumericalAlgorithms::Differentiation::Complex::forwardDifference(
            [&]( const VectorXcd& X_ ) {
              VectorXcd Xout = computeResidualVector( X_, FeTrial, alphaP );
              return Xout;
            },
            X );
          R = computeResidualVector( X, FeTrial, alphaP );
          counter += 1;
        }
      }
      catch ( std::exception& e ) {
        throw std::runtime_error( "return mapping failed: " + std::string( e.what() ) );
      }
      /* std::cout << "inner newton iters: " << counter << std::endl; */

      // update plastic deformation increment
      Fe              = X.segment( 0, 9 ).data();
      dFp             = Fastor::inverse( Fe ) % FeTrial;
      alphaP          = X( 9 );
      Tensor33d FpNew = dFp % Fp;
      memcpy( Fp.data(), FpNew.data(), 9 * sizeof( double ) );

      using namespace ContinuumMechanics;
      double      psi_;
      Tensor33d   Ce, dPsi_dCe;
      Tensor3333d dCe_dFe, d2Psi_dCedCe;
      std::tie( Ce, dCe_dFe ) = DeformationMeasures::FirstOrderDerived::rightCauchyGreen( Fe );

      // compute energy density, first and second partial derivatives wrt Cauchy
      // Green deformation
      using func_type    = std::function< Tensor33t< complexDouble >( const Tensor33t< complexDouble >& ) >;
      func_type computeS = [&]( const Tensor33t< complexDouble >& Ce_ ) {
        const auto [_psi, _dPsi_dCe] = EnergyDensityFunctions::FirstOrderDerived::PenceGouPotentialB( Ce_, K, G );
        return _dPsi_dCe;
      };

      d2Psi_dCedCe = NumericalAlgorithms::Differentiation::Complex::TensorToTensor::forwardDifference( computeS, Ce );
      std::tie( psi_, dPsi_dCe ) = EnergyDensityFunctions::FirstOrderDerived::PenceGouPotentialB( Ce, K, G );
      // compute Kirchhoff stress
      Tensor33d   PK2 = 2. * dPsi_dCe;
      Tensor3333d dTau_dPK2, dTau_dFe_partial;
      std::tie( response.tau,
                dTau_dPK2,
                dTau_dFe_partial )  = StressMeasures::FirstOrderDerived::KirchhoffStressFromPK2( PK2, Fe );
      response.rho                  = 1.0;
      response.elasticEnergyDensity = psi_;

      // compute tangent operator
      using mM9d = Eigen::Map< Eigen::Matrix< double, 9, 9 > >;

      MatrixXd dYdDeformation              = MatrixXd::Zero( 11, 11 );
      dYdDeformation.block< 9, 9 >( 0, 0 ) = mM9d( Tensor3333d( einsum< IK, JL, to_IJKL >( Spatial3D::I,
                                                                                           transpose( Fastor::inverse(
                                                                                             FpOld ) ) ) )
                                                     .data() )
                                               .transpose();
      MatrixXd dXdDeformation = dR_dX.colPivHouseholderQr().solve( dYdDeformation );

      Tensor3333d dFe_dF = Tensor3333d( Matrix9d( dXdDeformation.block< 9, 9 >( 0, 0 ).transpose() ).data() );

      Tensor3333d dPK2_dFe = einsum< ijKL, KLMN >( 2. * d2Psi_dCedCe, dCe_dFe );
      Tensor3333d dPK2_dF  = einsum< ijKL, KLMN >( dPK2_dFe, dFe_dF );

      /* tangents.dTau_dF = einsum< IJKL, KLMN >( dTau_dPK2, dPK2_dF ) +
       * dTau_dF_partial; */
      tangents.dTau_dF = einsum< IJKL, KLMN >( dTau_dPK2, dPK2_dF ) + einsum< ijKL, KLMN >( dTau_dFe_partial, dFe_dF );
    }
    else {
      using namespace Marmot::ContinuumMechanics;
      double      psi_;
      Tensor33d   Ce, dPsi_dCe;
      Tensor3333d dCe_dFe, d2Psi_dCedCe;
      std::tie( Ce, dCe_dFe ) = DeformationMeasures::FirstOrderDerived::rightCauchyGreen( Fe );

      // compute energy density, first and second partial derivatives wrt Cauchy
      // Green deformation
      using func_type    = std::function< Tensor33t< complexDouble >( const Tensor33t< complexDouble >& ) >;
      func_type computeS = [&]( const Tensor33t< complexDouble >& Ce_ ) {
        const auto [_psi, _dPsi_dCe] = EnergyDensityFunctions::FirstOrderDerived::PenceGouPotentialB( Ce_, K, G );
        return _dPsi_dCe;
      };

      d2Psi_dCedCe = NumericalAlgorithms::Differentiation::Complex::TensorToTensor::forwardDifference( computeS, Ce );

      std::tie( psi_, dPsi_dCe ) = EnergyDensityFunctions::FirstOrderDerived::PenceGouPotentialB( Ce, K, G );
      // compute Kirchhoff stress
      Tensor33d   PK2 = 2. * dPsi_dCe;
      Tensor3333d dTau_dPK2, dTau_dFe_partial;
      std::tie( response.tau,
                dTau_dPK2,
                dTau_dFe_partial )  = StressMeasures::FirstOrderDerived::KirchhoffStressFromPK2( PK2, Fe );
      response.rho                  = 1.0;
      response.elasticEnergyDensity = psi_;

      // compute tangent operator
      Tensor3333d dPK2_dFe = einsum< ijKL, KLMN >( 2. * d2Psi_dCedCe, dCe_dFe );
      Tensor3333d dFe_dF   = einsum< IK, JL, to_IJKL >( Spatial3D::I, transpose( Fastor::inverse( FpOld ) ) );
      Tensor3333d dPK2_dF  = einsum< ijKL, KLMN >( dPK2_dFe, dFe_dF );

      /* tangents.dTau_dF = einsum< IJKL, KLMN >( dTau_dPK2, dPK2_dF ) +
       * dTau_dF_partial; */
      tangents.dTau_dF = einsum< IJKL, KLMN >( dTau_dPK2, dPK2_dF ) + einsum< ijKL, KLMN >( dTau_dFe_partial, dFe_dF );
    }
  }

  void FiniteStrainJ2Plasticity::initializeYourself( double* stateVars, int nStateVars )
  {
    // set all state variables to zero
    for ( int i = 0; i < nStateVars; ++i ) {
      stateVars[i] = 0.0;
    }

    Tensor33d& Fp = stateLayout.getAs< Tensor33d& >( stateVars, "Fp" );
    Fp.eye();
  }
} // namespace Marmot::Materials
