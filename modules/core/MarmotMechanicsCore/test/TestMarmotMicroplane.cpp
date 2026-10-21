#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotMicroplane.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotVoigt.h"
#include <cmath>

using namespace Marmot::Testing;
using namespace Marmot::Microplane;
using namespace Marmot::ContinuumMechanics::VoigtNotation;
using namespace Marmot::ContinuumMechanics::Elasticity::Isotropic;

namespace {

  // Generic checks that must hold for every scheme, independent of the literal reference values:
  // correct row count, unit-length normals, and a total weight sum of 3 (the normalization used
  // by all 5 hard-coded schemes). Cheap, and localizes a broken table (wrong row count, a
  // non-normalized direction, a mis-scaled weight) far more precisely than the round-trip checks
  // below, which would just fail with "does not recover the original tensor".
  void checkGenericInvariants( QuadratureScheme scheme, const std::string& name )
  {
    const Discretization::Data data = Discretization::get( scheme );
    const size_t               nqp  = nqpFromScheme( scheme );

    throwExceptionOnFailure( static_cast< size_t >( data.normals.rows() ) == nqp, name + ": row count" );
    throwExceptionOnFailure( static_cast< size_t >( data.weights.rows() ) == nqp, name + ": weight count" );

    for ( size_t i = 0; i < nqp; i++ ) {
      const double normSq = data.normals.row( i ).squaredNorm();
      throwExceptionOnFailure( checkIfEqual( normSq, 1.0, 1e-6 ),
                               name + ": unit normal at row " + std::to_string( i ) );
    }

    throwExceptionOnFailure( checkIfEqual( data.weights.sum(), 3.0, 1e-6 ), name + ": weight sum" );
  }

  // Projects the macroscopic strain tensor eps onto every microplane's deviatoric strain vector
  // via VDSplit::pTensors.Dev and checks it against the independent, closed-form microplane theory formula
  // (eps_D^(m))_i = eps_ij n_j^(m) - eps_V n_i^(m). Then homogenizes the very same per-microplane volumetric/deviatoric
  // values back to a macroscopic tensor and checks that the original eps is recovered exactly: this is the completeness
  // property of the microplane quadrature.
  template < QuadratureScheme scheme >
  void checkVDSplitProjectionAndHomogenization( const Matrix3d& eps, const std::string& name )
  {
    VDSplit< scheme > vdSplit;
    constexpr size_t  nqp = VDSplit< scheme >::nqp;

    const double epsV = eps.trace() / 3.;

    // Independent, closed-form kinematic projection: row m = eps * n^(m) - eps_V * n^(m).
    const auto epsD_analytic = ( vdSplit.pTensors.nMatrix * eps - epsV * vdSplit.pTensors.nMatrix ).eval();

    // Actual kinematic projection via the Dev tensor, exactly as done in MDP.cpp/GMDP.cpp.
    TensorMap< const Tensor< double, 2 > > epsT( eps.data(), 3, 3 );
    Eigen::Tensor< double, 2 > epsD_actualT = vdSplit.pTensors.Dev.contract( epsT, contractionDims< 2, 0, 3, 1 >() );
    Eigen::Map< const Eigen::Matrix< double, Eigen::Dynamic, 3 > > epsD_actual( epsD_actualT.data(), nqp, 3 );

    throwExceptionOnFailure( checkIfEqual( Eigen::MatrixXd( epsD_actual ), Eigen::MatrixXd( epsD_analytic ), 1e-12 ),
                             name + ": Dev kinematic projection vs. eps*n - eps_V*n" );

    // Homogenize the per-microplane strains and check that the original macroscopic strain tensor
    // is recovered exactly.
    auto volValue = vdSplit.pTensors.w;
    volValue.setConstant( epsV );
    const Matrix3d recovered = vdSplit.homogenize( volValue, epsD_analytic );

    if constexpr ( scheme != QuadratureScheme::N61 ) {
      throwExceptionOnFailure( checkIfEqual( Eigen::MatrixXd( recovered ), Eigen::MatrixXd( eps ), 1e-8 ),
                               name + ": homogenize(project(eps)) recovers eps" );
    }

    // The weighted sum of the single-microplane contributions (applyProjectionTensors) must agree
    // with the batched homogenize() call above, regardless of whether the completeness identity
    // holds for this scheme.
    Matrix3d accumulated = Matrix3d::Zero();
    for ( size_t m = 0; m < nqp; m++ )
      accumulated += vdSplit.pTensors.w( m ) *
                     vdSplit.applyProjectionTensors( epsV, Vector3d( epsD_analytic.row( m ).transpose() ), m );

    throwExceptionOnFailure( checkIfEqual( Eigen::MatrixXd( accumulated ), Eigen::MatrixXd( recovered ), 1e-10 ),
                             name + ": weighted sum of applyProjectionTensors matches homogenize" );
  }

  // Projects the macroscopic strain tensor eps onto every microplane's normal/shear strains via
  // NTSplit::pTensors.N/M/L, exactly as MicroplaneM7.cpp/GM7.cpp do (dotting the engineering-Voigt
  // strain with each projection tensor column), then homogenizes those projections back via
  // homogenizeStress(). Since N_ij N_kl + M_ij M_kl + L_ij L_kl integrated over all microplanes is
  // the symmetric identity tensor, this round trip must recover the original tensor exactly - the
  // N-T split analogue of the V-D split completeness check above.
  template < QuadratureScheme scheme >
  void checkNTSplitProjectionAndHomogenization( const Matrix3d& eps, const std::string& name )
  {
    NTSplit< scheme > ntSplit;

    const Marmot::Vector6d epsVoigt = strainToVoigt( eps );

    const auto eps_N = ( epsVoigt.transpose() * ntSplit.pTensors.N ).eval();
    const auto eps_M = ( epsVoigt.transpose() * ntSplit.pTensors.M ).eval();
    const auto eps_L = ( epsVoigt.transpose() * ntSplit.pTensors.L ).eval();

    const Marmot::Vector6d recovered = ntSplit.homogenizeStress( eps_N, eps_M, eps_L );

    throwExceptionOnFailure( checkIfEqual( Eigen::MatrixXd( recovered ),
                                           Eigen::MatrixXd( stressToVoigt( eps ) ),
                                           1e-8 ),
                             name + ": homogenizeStress(project(eps)) recovers eps" );
  }

  // Verifies that a linear-elastic V-D split microplane constitutive law - sigma_V^(m) = 3K
  // deps_V and sigma_D^(m) = 2G deps_D^(m), the classical microplane elastic moduli reproducing
  // isotropic elasticity - homogenizes back to exactly the macroscopic isotropic elastic stress
  // stiffnessTensor(E,nu) * depsVoigt for the same strain increment. This is the actual
  // constitutive-model consistency check a real V-D split material relies on, not just a
  // kinematic identity.
  template < QuadratureScheme scheme >
  void checkVDSplitLinearElasticity( const Marmot::Vector6d& depsVoigt,
                                     const double            E,
                                     const double            nu,
                                     const std::string&      name )
  {
    VDSplit< scheme > vdSplit;

    const double G = shearModulus( E, nu );
    const double K = lameParameter( E, nu ) + 2. / 3. * G;

    const Matrix3d deps  = voigtToStrain( depsVoigt );
    const double   depsV = deps.trace() / 3.;
    const auto     depsD = ( vdSplit.pTensors.nMatrix * deps - depsV * vdSplit.pTensors.nMatrix ).eval();

    auto sigV = vdSplit.pTensors.w;
    sigV.setConstant( 3. * K * depsV );
    const auto sigD = ( 2. * G * depsD ).eval();

    const Matrix3d         stress_microplane = vdSplit.homogenize( sigV, sigD );
    const Marmot::Vector6d stress_macro      = stiffnessTensor( E, nu ) * depsVoigt;

    throwExceptionOnFailure( checkIfEqual( Eigen::MatrixXd( stressToVoigt( stress_microplane ) ),
                                           Eigen::MatrixXd( stress_macro ),
                                           1e-4 ),
                             name + ": VDSplit linear elasticity matches macroscopic isotropic stiffness" );
  }

  // N-T split analogue of checkVDSplitLinearElasticity: sigma_N^(m) = E_N0 deps_N^(m),
  // sigma_M^(m) = E_T deps_M^(m), sigma_L^(m) = E_T deps_L^(m) with the classical N-T split
  // elastic moduli (as used e.g. by MicroplaneM7's elastic response), homogenized via
  // homogenizeStress() and compared to the macroscopic isotropic elastic stress.
  template < QuadratureScheme scheme >
  void checkNTSplitLinearElasticity( const Marmot::Vector6d& depsVoigt,
                                     const double            E,
                                     const double            nu,
                                     const std::string&      name )
  {
    NTSplit< scheme > ntSplit;

    const double E_N0 = E / ( 1.0 - 2.0 * nu );
    const double E_T  = E / ( 1.0 - 2.0 * nu ) * ( 1.0 - 4.0 * nu ) / ( 1.0 + nu );

    const auto deps_N = ( depsVoigt.transpose() * ntSplit.pTensors.N ).eval();
    const auto deps_M = ( depsVoigt.transpose() * ntSplit.pTensors.M ).eval();
    const auto deps_L = ( depsVoigt.transpose() * ntSplit.pTensors.L ).eval();

    const auto sig_N = ( E_N0 * deps_N ).eval();
    const auto sig_M = ( E_T * deps_M ).eval();
    const auto sig_L = ( E_T * deps_L ).eval();

    const Marmot::Vector6d stress_microplane = ntSplit.homogenizeStress( sig_N, sig_M, sig_L );
    const Marmot::Vector6d stress_macro      = stiffnessTensor( E, nu ) * depsVoigt;

    // See checkVDSplitLinearElasticity for why this tolerance is looser than machine precision.
    throwExceptionOnFailure( checkIfEqual( Eigen::MatrixXd( stress_microplane ),
                                           Eigen::MatrixXd( stress_macro ),
                                           1e-4 ),
                             name + ": NTSplit linear elasticity matches macroscopic isotropic stiffness" );
  }

} // namespace

void testDiscretizationSchemes()
{
  checkGenericInvariants( QuadratureScheme::N21Symmetric, "N21Symmetric" );
  checkGenericInvariants( QuadratureScheme::N21Asymmetric, "N21Asymmetric" );
  checkGenericInvariants( QuadratureScheme::N28, "N28" );
  checkGenericInvariants( QuadratureScheme::N37, "N37" );
  checkGenericInvariants( QuadratureScheme::N61, "N61" );
}

void testVDSplitProjectionAndHomogenization()
{
  Matrix3d eps;
  // clang-format off
  eps << 0.020,  0.010, -0.015,
         0.010, -0.030,  0.008,
        -0.015,  0.008,  0.015;
  // clang-format on

  checkVDSplitProjectionAndHomogenization< QuadratureScheme::N21Symmetric >( eps, "N21Symmetric" );
  checkVDSplitProjectionAndHomogenization< QuadratureScheme::N21Asymmetric >( eps, "N21Asymmetric" );
  checkVDSplitProjectionAndHomogenization< QuadratureScheme::N28 >( eps, "N28" );
  checkVDSplitProjectionAndHomogenization< QuadratureScheme::N37 >( eps, "N37" );
  checkVDSplitProjectionAndHomogenization< QuadratureScheme::N61 >( eps, "N61" );
}

void testNTSplitProjectionAndHomogenization()
{
  Matrix3d eps;
  // clang-format off
  eps << 0.020,  0.010, -0.015,
         0.010, -0.030,  0.008,
        -0.015,  0.008,  0.015;
  // clang-format on

  checkNTSplitProjectionAndHomogenization< QuadratureScheme::N21Symmetric >( eps, "N21Symmetric" );
  checkNTSplitProjectionAndHomogenization< QuadratureScheme::N21Asymmetric >( eps, "N21Asymmetric" );
  checkNTSplitProjectionAndHomogenization< QuadratureScheme::N28 >( eps, "N28" );
  checkNTSplitProjectionAndHomogenization< QuadratureScheme::N37 >( eps, "N37" );
  checkNTSplitProjectionAndHomogenization< QuadratureScheme::N61 >( eps, "N61" );
}

void testVDSplitLinearElasticity()
{
  const double E  = 30000.;
  const double nu = 0.2;

  Marmot::Vector6d depsVoigt;
  depsVoigt << 0.0010, 0.0004, -0.0006, 0.0016, -0.0008, 0.0012;

  checkVDSplitLinearElasticity< QuadratureScheme::N21Symmetric >( depsVoigt, E, nu, "N21Symmetric" );
  checkVDSplitLinearElasticity< QuadratureScheme::N21Asymmetric >( depsVoigt, E, nu, "N21Asymmetric" );
  checkVDSplitLinearElasticity< QuadratureScheme::N28 >( depsVoigt, E, nu, "N28" );
  checkVDSplitLinearElasticity< QuadratureScheme::N37 >( depsVoigt, E, nu, "N37" );
  checkVDSplitLinearElasticity< QuadratureScheme::N61 >( depsVoigt, E, nu, "N61" );
}

void testNTSplitLinearElasticity()
{
  const double E  = 30000.;
  const double nu = 0.2;

  Marmot::Vector6d depsVoigt;
  depsVoigt << 0.0010, 0.0004, -0.0006, 0.0016, -0.0008, 0.0012;

  checkNTSplitLinearElasticity< QuadratureScheme::N21Symmetric >( depsVoigt, E, nu, "N21Symmetric" );
  checkNTSplitLinearElasticity< QuadratureScheme::N21Asymmetric >( depsVoigt, E, nu, "N21Asymmetric" );
  checkNTSplitLinearElasticity< QuadratureScheme::N28 >( depsVoigt, E, nu, "N28" );
  checkNTSplitLinearElasticity< QuadratureScheme::N37 >( depsVoigt, E, nu, "N37" );
  checkNTSplitLinearElasticity< QuadratureScheme::N61 >( depsVoigt, E, nu, "N61" );
}

int main()
{
  executeTestsAndCollectExceptions( { testDiscretizationSchemes,
                                      testVDSplitProjectionAndHomogenization,
                                      testNTSplitProjectionAndHomogenization,
                                      testVDSplitLinearElasticity,
                                      testNTSplitLinearElasticity } );

  return 0;
}
