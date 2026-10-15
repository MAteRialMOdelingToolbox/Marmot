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
      hasNetworkC( nMaterialProperties >= 20 ),
      kappaC( hasNetworkC ? materialProperties[13] : 0.0 ),
      C1( hasNetworkC ? materialProperties[14] : 0.0 ),
      C2( hasNetworkC ? materialProperties[15] : 0.0 ),
      C3( hasNetworkC ? materialProperties[16] : 0.0 ),
      c4( hasNetworkC ? materialProperties[17] : 0.0 ),
      c5( hasNetworkC ? materialProperties[18] : 0.0 ),
      c6( hasNetworkC ? materialProperties[19] : 0.0 ),
      density( hasNetworkC ? ( nMaterialProperties > 20 ? materialProperties[20] : 0.0 )
                           : ( nMaterialProperties > 13 ? materialProperties[13] : 0.0 ) )
  {
    stateLayout.add( "Fv", 9 );    // viscous deformation gradient of network B
    if ( hasNetworkC ) {
      stateLayout.add( "FvC", 9 ); // viscous deformation gradient of network C
    }
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

  BergstromBoyce::ViscousNetworkResponse BergstromBoyce::solveViscousNetwork( const Tensor33d& F,
                                                                              Tensor33d&       Fv,
                                                                              const double     dT,
                                                                              const double     kappa,
                                                                              const double     p1,
                                                                              const double     p2,
                                                                              const double     p3,
                                                                              const double     c1_,
                                                                              const double     c2_,
                                                                              const double     c3_ ) const
  {
    using namespace Marmot;
    using namespace Fastor;
    using namespace Eigen;
    using namespace FastorIndices;
    using namespace FastorStandardTensors;
    using complexDouble = std::complex< double >;

    const Tensor33d FvOld( Fv );
    Tensor33d       FeTrial = F % Fastor::inverse( FvOld );

    // --- solve the local 10x10 Newton problem for (Fe, dGamma) ---
    using mV9d = Eigen::Map< Eigen::Matrix< double, 9, 1 > >;
    using mM9d = Eigen::Map< Eigen::Matrix< double, 9, 9 > >;

    VectorXd X( 10 );
    X.segment( 0, 9 )  = mV9d( FeTrial.data() );
    X( 9 )             = 0.0;
    VectorXd        dX = VectorXd::Zero( 10 );
    VectorXd        R  = VectorXd::Zero( 10 );
    Eigen::MatrixXd dR_dX( 10, 10 );

    auto evalResidual = [&]( const VectorXcd& X_ ) {
      return computeResidualVector( X_, FeTrial, dT, kappa, p1, p2, p3, c1_, c2_, c3_ );
    };

    std::tie( R, dR_dX ) = NumericalAlgorithms::Differentiation::Complex::forwardDifference( evalResidual, X );
    R                    = computeResidualVector( X, FeTrial, dT, kappa, p1, p2, p3, c1_, c2_, c3_ );

    size_t counter = 0;
    try {
      while ( R.norm() > 1e-12 || dX.norm() > 1e-12 ) {

        if ( counter > 20 )
          throw StressUpdateFailed( "BergstromBoyce: inner newton not converged" );

        dX = -dR_dX.colPivHouseholderQr().solve( R );
        X += dX;
        std::tie( R, dR_dX ) = NumericalAlgorithms::Differentiation::Complex::forwardDifference( evalResidual, X );
        R                    = computeResidualVector( X, FeTrial, dT, kappa, p1, p2, p3, c1_, c2_, c3_ );
        counter += 1;
      }
    }
    catch ( std::exception& e ) {
      throw StressUpdateFailed( "BergstromBoyce: return mapping failed: " + std::string( e.what() ) );
    }

    Tensor33d Fe( X.segment( 0, 9 ).data() );

    // --- update the viscous deformation gradient ---
    Tensor33d dFv   = Fastor::inverse( Fe ) % FeTrial;
    Tensor33d FvNew = dFv % FvOld;
    Fv              = FvNew;

    // --- evaluation on the elastic Ce, pushed forward through Fe ---
    using namespace ContinuumMechanics;
    Tensor33d   Ce, dPsi_dCe;
    Tensor3333d dCe_dFe;
    std::tie( Ce, dCe_dFe ) = DeformationMeasures::FirstOrderDerived::rightCauchyGreen( Fe );

    double psi;
    std::tie( psi, dPsi_dCe ) = hyperelasticPotential( Ce, hyperelasticBase, p1, p2, p3, kappa );

    using func_type    = std::function< Tensor33t< complexDouble >( const Tensor33t< complexDouble >& ) >;
    func_type computeS = [&]( const Tensor33t< complexDouble >& Ce_ ) {
      const auto [_psi, _dPsi_dCe] = hyperelasticPotential( Ce_, hyperelasticBase, p1, p2, p3, kappa );
      return _dPsi_dCe;
    };
    Tensor3333d
      d2Psi_dCedCe = NumericalAlgorithms::Differentiation::Complex::TensorToTensor::forwardDifference( computeS, Ce );

    Tensor33d   PK2 = 2. * dPsi_dCe;
    Tensor3333d dTau_dPK2, dTau_dFe_partial;
    Tensor33d   tau;
    std::tie( tau, dTau_dPK2, dTau_dFe_partial ) = StressMeasures::FirstOrderDerived::KirchhoffStressFromPK2( PK2, Fe );

    // --- implicit-function-theorem sensitivity dFe/dF ---
    MatrixXd dYdDeformation              = MatrixXd::Zero( 10, 10 );
    dYdDeformation.block< 9, 9 >( 0, 0 ) = mM9d( Tensor3333d(
                                                   einsum< IK, JL, to_IJKL >( Spatial3D::I,
                                                                              transpose( Fastor::inverse( FvOld ) ) ) )
                                                   .data() )
                                             .transpose();
    MatrixXd dXdDeformation = dR_dX.colPivHouseholderQr().solve( dYdDeformation );

    Tensor3333d dFe_dF = Tensor3333d( Matrix9d( dXdDeformation.block< 9, 9 >( 0, 0 ).transpose() ).data() );

    Tensor3333d dPK2_dFe = einsum< ijKL, KLMN >( 2. * d2Psi_dCedCe, dCe_dFe );
    Tensor3333d dPK2_dF  = einsum< ijKL, KLMN >( dPK2_dFe, dFe_dF );

    Tensor3333d dTau_dF = einsum< IJKL, KLMN >( dTau_dPK2, dPK2_dF ) + einsum< ijKL, KLMN >( dTau_dFe_partial, dFe_dF );

    Tensor33d N;
    double    rho, gammaDot;
    std::tie( N, rho, gammaDot ) = computeFlowQuantities( Fe, kappa, p1, p2, p3, c1_, c2_, c3_ );
    const double dGamma          = X( 9 );
    const double dissipation     = dGamma * rho;

    return { tau, dTau_dF, psi, dissipation };
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

    const double dT = timeIncrement.dT;

    // --- network A: direct evaluation on total C ---
    using namespace ContinuumMechanics;
    Tensor33d   C, dPsiA_dC;
    Tensor3333d dC_dF;
    std::tie( C, dC_dF ) = DeformationMeasures::FirstOrderDerived::rightCauchyGreen( deformation.F );

    double psiA;
    std::tie( psiA, dPsiA_dC ) = hyperelasticPotential( C, hyperelasticBase, A1, A2, A3, kappaA );

    using func_type_A     = std::function< Tensor33t< complexDouble >( const Tensor33t< complexDouble >& ) >;
    func_type_A computeSA = [&]( const Tensor33t< complexDouble >& C_ ) {
      const auto [_psi, _dPsi_dC] = hyperelasticPotential( C_, hyperelasticBase, A1, A2, A3, kappaA );
      return _dPsi_dC;
    };
    Tensor3333d
      d2PsiA_dCdC = NumericalAlgorithms::Differentiation::Complex::TensorToTensor::forwardDifference( computeSA, C );

    Tensor33d   PK2_A = 2. * dPsiA_dC;
    Tensor3333d dTauA_dPK2, dTauA_dF_partial;
    Tensor33d   tauA;
    std::tie( tauA,
              dTauA_dPK2,
              dTauA_dF_partial ) = StressMeasures::FirstOrderDerived::KirchhoffStressFromPK2( PK2_A, deformation.F );

    Tensor3333d dPK2A_dF = einsum< ijKL, KLMN >( 2. * d2PsiA_dCdC, dC_dF );
    // dTauA_dF_partial is already the direct partial derivative w.r.t. the total F (the second
    // argument KirchhoffStressFromPK2 was called with above), so it is added directly -- it must
    // NOT be chain-ruled through dC_dF again (that tensor maps derivatives w.r.t. C, not F).
    Tensor3333d dTauA_dF = einsum< IJKL, KLMN >( dTauA_dPK2, dPK2A_dF ) + dTauA_dF_partial;

    // --- network B: evaluation on elastic Ce_B, pushed forward through Fe_B ---
    TensorMap33d FvB_map = stateLayout.getAs< TensorMap33d >( response.stateVars, "Fv" );
    Tensor33d    FvB( FvB_map );
    auto         resB = solveViscousNetwork( deformation.F, FvB, dT, kappaB, B1, B2, B3, c1, c2, c3 );
    memcpy( FvB_map.data(), FvB.data(), 9 * sizeof( double ) );

    // --- network C (optional): evaluation on elastic Ce_C, pushed forward through Fe_C ---
    ViscousNetworkResponse resC{};
    if ( hasNetworkC ) {
      TensorMap33d FvC_map = stateLayout.getAs< TensorMap33d >( response.stateVars, "FvC" );
      Tensor33d    FvC( FvC_map );
      resC = solveViscousNetwork( deformation.F, FvC, dT, kappaC, C1, C2, C3, c4, c5, c6 );
      memcpy( FvC_map.data(), FvC.data(), 9 * sizeof( double ) );
    }

    // --- assemble total response ---
    response.tau                  = tauA + resB.tau + resC.tau;
    response.elasticEnergyDensity = psiA + resB.psi + resC.psi;

    // Assign, not accumulate: the solver reuses the same response object across trial Newton
    // iterations within one increment and does not reset dissipation between them (unlike tau
    // and elasticEnergyDensity, which it does reset), so += would make the reported dissipation
    // depend on iteration count rather than this constitutive update alone.
    response.dissipation = resB.dissipation + resC.dissipation;

    tangents.dTau_dF = dTauA_dF + resB.dTau_dF + resC.dTau_dF;
  }

  void BergstromBoyce::initializeYourself( double* stateVars, int nStateVars )
  {
    for ( int i = 0; i < nStateVars; ++i ) {
      stateVars[i] = 0.0;
    }

    TensorMap33d Fv = stateLayout.getAs< TensorMap33d >( stateVars, "Fv" );
    memcpy( Fv.data(), Spatial3D::I.data(), 9 * sizeof( double ) );

    if ( hasNetworkC ) {
      TensorMap33d FvC = stateLayout.getAs< TensorMap33d >( stateVars, "FvC" );
      memcpy( FvC.data(), Spatial3D::I.data(), 9 * sizeof( double ) );
    }
  }
} // namespace Marmot::Materials
