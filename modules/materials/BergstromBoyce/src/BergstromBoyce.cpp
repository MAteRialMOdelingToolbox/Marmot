#include "Marmot/BergstromBoyce.h"
#include "Marmot/MarmotDeformationMeasures.h"
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include "Marmot/MarmotNumericalDifferentiation.h"
#include "Marmot/MarmotNumericalDifferentiationForFastor.h"
#include "Marmot/MarmotStressMeasures.h"
#include <complex>

namespace Marmot::Materials {

  using namespace Marmot;
  using namespace Fastor;
  using namespace FastorIndices;
  using namespace FastorStandardTensors;

  BergstromBoyce::BergstromBoyce( const double* materialProperties, int nMaterialProperties, int materialLabel )
    : MarmotMaterialFiniteStrain( materialProperties, nMaterialProperties, materialLabel ),
      hyperelasticBase( materialProperties[0] ),
      kappaA( materialProperties[1] ),
      kappaB( materialProperties[2] ),
      A1( materialProperties[3] ),
      A2( materialProperties[4] ),
      A3( materialProperties[5] ),
      B1( materialProperties[6] ),
      B2( materialProperties[7] ),
      B3( materialProperties[8] ),
      c1( materialProperties[9] ),
      c2( materialProperties[10] ),
      c3( materialProperties[11] ),
      implementationType( materialProperties[12] ),
      density( nMaterialProperties > 13 ? materialProperties[13] : 0.0 )
  {
    stateLayout.add( "Fv", 9 ); // viscous deformation gradient of network B
    stateLayout.finalize();
  }

  void BergstromBoyce::computeStress( ConstitutiveResponse< 3 >& response,
                                      AlgorithmicModuli< 3 >&    tangents,
                                      const Deformation< 3 >&    deformation,
                                      const TimeIncrement&       timeIncrement ) const
  {
    switch ( implementationType ) {
    case 0: computeStressCSDA( response, tangents, deformation, timeIncrement ); break;
    case 1: computeStressWithFullReturnMapping( response, tangents, deformation, timeIncrement ); break;
    default: throw std::invalid_argument( "implementation type not supported" );
    };
  }

  void BergstromBoyce::computeStressWithFullReturnMapping( ConstitutiveResponse< 3 >& response,
                                                           AlgorithmicModuli< 3 >&    tangents,
                                                           const Deformation< 3 >&    deformation,
                                                           const TimeIncrement&       timeIncrement ) const
  {
    throw std::invalid_argument( "not implemented yet -- use CSDA (implementationType = 0)" );
  }

  void BergstromBoyce::computeStressCSDA( ConstitutiveResponse< 3 >& response,
                                          AlgorithmicModuli< 3 >&    tangents,
                                          const Deformation< 3 >&    deformation,
                                          const TimeIncrement&       timeIncrement ) const
  {
    using namespace Marmot;
    using namespace Fastor;
    using namespace Eigen;
    using namespace FastorIndices;
    using namespace FastorStandardTensors;
    using complexDouble = std::complex< double >;

    TensorMap33d    Fv = stateLayout.getAs< TensorMap33d >( response.stateVars, "Fv" );
    const Tensor33d FvOld( Fv );

    const double dT = timeIncrement.dT;

    Tensor33d FeTrial = deformation.F % Fastor::inverse( FvOld );

    // --- solve the local 10x10 Newton problem for (Fe, dGamma) ---
    using mV9d = Eigen::Map< Eigen::Matrix< double, 9, 1 > >;
    using mM9d = Eigen::Map< Eigen::Matrix< double, 9, 9 > >;

    VectorXd X( 10 );
    X.segment( 0, 9 ) = mV9d( FeTrial.data() );
    X( 9 )             = 0.0;
    VectorXd        dX = VectorXd::Zero( 10 );
    VectorXd        R  = VectorXd::Zero( 10 );
    Eigen::MatrixXd dR_dX( 10, 10 );

    std::tie( R, dR_dX ) = NumericalAlgorithms::Differentiation::Complex::forwardDifference(
      [&]( const VectorXcd& X_ ) { return computeResidualVector( X_, FeTrial, dT ); },
      X );
    R = computeResidualVector( X, FeTrial, dT );

    size_t counter = 0;
    try {
      while ( R.norm() > 1e-12 || dX.norm() > 1e-12 ) {

        if ( counter > 20 )
          throw StressUpdateFailed( "BergstromBoyce: inner newton not converged" );

        dX = -dR_dX.colPivHouseholderQr().solve( R );
        X += dX;
        std::tie( R, dR_dX ) = NumericalAlgorithms::Differentiation::Complex::forwardDifference(
          [&]( const VectorXcd& X_ ) { return computeResidualVector( X_, FeTrial, dT ); },
          X );
        R = computeResidualVector( X, FeTrial, dT );
        counter += 1;
      }
    }
    catch ( std::exception& e ) {
      throw StressUpdateFailed( "BergstromBoyce: return mapping failed: " + std::string( e.what() ) );
    }

    Tensor33d Fe( X.segment( 0, 9 ).data() );

    // --- update the viscous deformation gradient ---
    Tensor33d dFv    = Fastor::inverse( Fe ) % FeTrial;
    Tensor33d FvNew  = dFv % FvOld;
    memcpy( Fv.data(), FvNew.data(), 9 * sizeof( double ) );

    // --- network A: direct evaluation on total C ---
    using namespace ContinuumMechanics;
    Tensor33d   C, dPsiA_dC;
    Tensor3333d dC_dF;
    std::tie( C, dC_dF ) = DeformationMeasures::FirstOrderDerived::rightCauchyGreen( deformation.F );

    double psiA;
    std::tie( psiA, dPsiA_dC ) = hyperelasticPotential( C, hyperelasticBase, A1, A2, A3, kappaA );

    using func_type_A    = std::function< Tensor33t< complexDouble >( const Tensor33t< complexDouble >& ) >;
    func_type_A computeSA = [&]( const Tensor33t< complexDouble >& C_ ) {
      const auto [_psi, _dPsi_dC] = hyperelasticPotential( C_, hyperelasticBase, A1, A2, A3, kappaA );
      return _dPsi_dC;
    };
    Tensor3333d d2PsiA_dCdC = NumericalAlgorithms::Differentiation::Complex::TensorToTensor::forwardDifference(
      computeSA,
      C );

    Tensor33d   PK2_A = 2. * dPsiA_dC;
    Tensor3333d dTauA_dPK2, dTauA_dF_partial;
    Tensor33d   tauA;
    std::tie( tauA, dTauA_dPK2, dTauA_dF_partial ) = StressMeasures::FirstOrderDerived::KirchhoffStressFromPK2(
      PK2_A,
      deformation.F );

    Tensor3333d dPK2A_dF = einsum< ijKL, KLMN >( 2. * d2PsiA_dCdC, dC_dF );
    Tensor3333d dTauA_dF = einsum< IJKL, KLMN >( dTauA_dPK2, dPK2A_dF ) +
                           einsum< ijKL, KLMN >( dTauA_dF_partial, dC_dF );

    // --- network B: evaluation on the elastic Ce, pushed forward through Fe ---
    Tensor33d   Ce, dPsiB_dCe;
    Tensor3333d dCe_dFe;
    std::tie( Ce, dCe_dFe ) = DeformationMeasures::FirstOrderDerived::rightCauchyGreen( Fe );

    double psiB;
    std::tie( psiB, dPsiB_dCe ) = hyperelasticPotential( Ce, hyperelasticBase, B1, B2, B3, kappaB );

    using func_type_B    = std::function< Tensor33t< complexDouble >( const Tensor33t< complexDouble >& ) >;
    func_type_B computeSB = [&]( const Tensor33t< complexDouble >& Ce_ ) {
      const auto [_psi, _dPsi_dCe] = hyperelasticPotential( Ce_, hyperelasticBase, B1, B2, B3, kappaB );
      return _dPsi_dCe;
    };
    Tensor3333d d2PsiB_dCedCe = NumericalAlgorithms::Differentiation::Complex::TensorToTensor::forwardDifference(
      computeSB,
      Ce );

    Tensor33d   PK2_B = 2. * dPsiB_dCe;
    Tensor3333d dTauB_dPK2, dTauB_dFe_partial;
    Tensor33d   tauB;
    std::tie( tauB, dTauB_dPK2, dTauB_dFe_partial ) = StressMeasures::FirstOrderDerived::KirchhoffStressFromPK2(
      PK2_B,
      Fe );

    // --- implicit-function-theorem sensitivity dFe/dF ---
    MatrixXd dYdDeformation               = MatrixXd::Zero( 10, 10 );
    dYdDeformation.block< 9, 9 >( 0, 0 ) = mM9d( Tensor3333d( einsum< IK, JL, to_IJKL >( Spatial3D::I,
                                                                                         transpose( Fastor::inverse(
                                                                                           FvOld ) ) ) )
                                                   .data() )
                                             .transpose();
    MatrixXd dXdDeformation = dR_dX.colPivHouseholderQr().solve( dYdDeformation );

    Tensor3333d dFe_dF = Tensor3333d( Matrix9d( dXdDeformation.block< 9, 9 >( 0, 0 ).transpose() ).data() );

    Tensor3333d dPK2B_dFe = einsum< ijKL, KLMN >( 2. * d2PsiB_dCedCe, dCe_dFe );
    Tensor3333d dPK2B_dF  = einsum< ijKL, KLMN >( dPK2B_dFe, dFe_dF );

    Tensor3333d dTauB_dF = einsum< IJKL, KLMN >( dTauB_dPK2, dPK2B_dF ) +
                           einsum< ijKL, KLMN >( dTauB_dFe_partial, dFe_dF );

    // --- assemble total response ---
    response.tau                  = tauA + tauB;
    response.elasticEnergyDensity = psiA + psiB;

    Tensor33d N;
    double    rho, gammaDot;
    std::tie( N, rho, gammaDot ) = computeFlowQuantities( Fe );
    const double dGamma          = X( 9 );
    response.dissipation += dGamma * rho;

    tangents.dTau_dF = dTauA_dF + dTauB_dF;
  }

  void BergstromBoyce::initializeYourself( double* stateVars, int nStateVars )
  {
    for ( int i = 0; i < nStateVars; ++i ) {
      stateVars[i] = 0.0;
    }

    TensorMap33d Fv = stateLayout.getAs< TensorMap33d >( stateVars, "Fv" );
    memcpy( Fv.data(), Spatial3D::I.data(), 9 * sizeof( double ) );
  }
} // namespace Marmot::Materials
