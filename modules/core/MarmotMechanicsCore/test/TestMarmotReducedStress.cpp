#include "Marmot/LinearElastic.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotLowerDimensionalStress.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/RankineSofteningPlasticity.h"
#include "Marmot/VonMises.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

/* Tests of Marmot::ContinuumMechanics::ReducedStress (condense(), computeReducedStress()) and of the reduced stress
 * states of MaterialHypoElastic built on it. All stiffnesses and stresses are in Pa. */

using namespace Marmot;
using namespace Marmot::Testing;
using namespace Marmot::ContinuumMechanics::ReducedStress;
using Marmot::Math::Matrix6d;
using Marmot::Math::Vector6d;

namespace {

  constexpr double E  = 30e9; // Pa
  constexpr double nu = 0.2;

  double maxAbs( const Eigen::MatrixXd& m )
  {
    return m.cwiseAbs().maxCoeff();
  }

  /// relative comparison: max |a - b| <= tol * max |b|
  bool checkRelative( const Eigen::MatrixXd& a, const Eigen::MatrixXd& b, double tol )
  {
    return a.rows() == b.rows() && a.cols() == b.cols() && maxAbs( a - b ) <= tol * maxAbs( b );
  }

  std::string here( const char* function, const std::string& what )
  {
    return std::string( function ) + ": " + what;
  }

  /// an anisotropic stiffness: isotropic plus couplings between normal and shear components, among them the
  /// out-of-plane shear components 13 and 23
  Matrix6d anisotropicStiffness()
  {
    Matrix6d C      = Marmot::ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( E, nu );
    auto     couple = [&]( int i, int j, double value ) {
      C( i, j ) = value;
      C( j, i ) = value;
    };
    couple( 0, 4, 4e9 );   // 11 - 13
    couple( 1, 5, -3e9 );  // 22 - 23
    couple( 2, 4, 2e9 );   // 33 - 13
    couple( 0, 3, 1.5e9 ); // 11 - 12
    couple( 3, 5, 1e9 );   // 12 - 23
    return C;
  }

  /// the condensed tangent from the compliance, independent of condense(): for p and f covering all six components
  /// it is the inverse of the compliance block S_pp
  template < std::size_t nP_ >
  auto condensedFromCompliance( const Matrix6d& C, const std::array< int, nP_ >& p )
  {
    const Matrix6d                  S  = C.inverse();
    constexpr int                   nP = static_cast< int >( nP_ );
    Eigen::Matrix< double, nP, nP > Spp;
    for ( int i = 0; i < nP; i++ )
      for ( int j = 0; j < nP; j++ )
        Spp( i, j ) = S( p[i], p[j] );
    return Eigen::Matrix< double, nP, nP >( Spp.inverse() );
  }

  /// a linear 3D stress update from a zero stress, which records the strain of its last call
  struct LinearStressUpdate {
    Matrix6d  C;
    Vector6d* lastStrain = nullptr;

    void operator()( const Vector6d& dStrain, Vector6d& stress, Matrix6d& dStress_dStrain ) const
    {
      stress          = C * dStrain;
      dStress_dStrain = C;
      if ( lastStrain )
        *lastStrain = dStrain;
    }
  };

  /// a nonlinear elastic 3D stress update: sigma = C eps + k E eps^3 (componentwise), so that Newton iterates
  struct CubicStressUpdate {
    Matrix6d  C;
    double    k          = 2e3;
    Vector6d* lastStrain = nullptr;

    void operator()( const Vector6d& dStrain, Vector6d& stress, Matrix6d& dStress_dStrain ) const
    {
      stress          = C * dStrain + k * E * dStrain.array().cube().matrix();
      dStress_dStrain = C;
      dStress_dStrain.diagonal() += 3 * k * E * dStrain.array().square().matrix();
      if ( lastStrain )
        *lastStrain = dStrain;
    }
  };

  /// the converged stress of a reduced stress state of a stress update
  template < typename StressUpdate3D >
  Vector6d stressAt( const StressUpdate3D& update, const Vector6d& dStrain )
  {
    Vector6d stress;
    Matrix6d C;
    update( dStrain, stress, C );
    return stress;
  }

  template < std::size_t nF >
  double maxAbsAt( const Vector6d& v, const std::array< int, nF >& f )
  {
    double m = 0;
    for ( int i : f )
      m = std::max( m, std::abs( v( i ) ) );
    return m;
  }

} // namespace

// plane stress with normal - out-of-plane shear coupling: sigma_33 = sigma_13 = sigma_23 = 0, the condensed tangent
// is the exact Schur complement and differs from the one with gamma13 = gamma23 = 0 held (the assumption of master)
void testPlaneStressAnisotropicLinear()
{
  const Matrix6d     C = anisotropicStiffness();
  Vector6d           lastStrain;
  LinearStressUpdate update{ C, &lastStrain };

  Vector6d dStrain          = Vector6d::Zero();
  dStrain( 0 )              = 1e-4;
  dStrain( 1 )              = -3e-5;
  dStrain( 3 )              = 5e-5;
  dStrain( 2 )              = -dStrain( 0 ) - dStrain( 1 );
  const Vector6d prescribed = dStrain;

  const auto condensation = computeReducedStress( update, dStrain, planeStressPrescribed, planeStressFree );

  const Vector6d stress = stressAt( update, dStrain );
  throwExceptionOnFailure( maxAbsAt( stress, planeStressFree ) <= 1e-12 * maxAbs( stress ),
                           here( __PRETTY_FUNCTION__, "sigma_f does not vanish" ) );
  for ( int i : planeStressPrescribed )
    throwExceptionOnFailure( dStrain( i ) == prescribed( i ),
                             here( __PRETTY_FUNCTION__, "a prescribed strain was changed" ) );
  throwExceptionOnFailure( lastStrain == dStrain,
                           here( __PRETTY_FUNCTION__, "the last evaluation is not at the converged strain" ) );
  throwExceptionOnFailure( condensation.nIterations == 2,
                           here( __PRETTY_FUNCTION__, "a linear update must converge after one correction" ) );

  const Eigen::Matrix3d exact = condensedFromCompliance( C, planeStressPrescribed );
  throwExceptionOnFailure( checkRelative( condensation.dStressP_dStrainP, exact, 1e-12 ),
                           here( __PRETTY_FUNCTION__, "condensed tangent != Schur complement" ) );

  // the in-plane stress is the condensed tangent times the in-plane strain
  const Eigen::Vector3d strainP( prescribed( 0 ), prescribed( 1 ), prescribed( 3 ) );
  const Eigen::Vector3d stressP( stress( 0 ), stress( 1 ), stress( 3 ) );
  throwExceptionOnFailure( checkRelative( stressP, exact * strainP, 1e-12 ),
                           here( __PRETTY_FUNCTION__, "in-plane stress != condensed tangent * strain" ) );

  // the sensitivity: C_f: [ I ; deps_f/deps_p ] = 0
  Eigen::Matrix< double, 6, 3 > dStrain_dStrainP = Eigen::Matrix< double, 6, 3 >::Zero();
  for ( int i = 0; i < 3; i++ ) {
    dStrain_dStrainP( planeStressPrescribed[i], i ) = 1;
    dStrain_dStrainP.row( planeStressFree[i] )      = condensation.dStrainF_dStrainP.row( i );
  }
  const Eigen::Matrix< double, 6, 3 > dStress_dStrainP = C * dStrain_dStrainP;
  for ( int i : planeStressFree )
    throwExceptionOnFailure( dStress_dStrainP.row( i ).cwiseAbs().maxCoeff() <= 1e-12 * maxAbs( C ),
                             here( __PRETTY_FUNCTION__, "dStrainF_dStrainP does not keep sigma_f = 0" ) );

  // master held gamma13 = gamma23 = 0 and condensed only eps33: wrong for this material
  const auto onlyEps33 = condense( C, planeStressPrescribed, std::array< int, 1 >{ 2 } );
  throwExceptionOnFailure( !checkRelative( onlyEps33.dStressP_dStrainP, exact, 1e-3 ),
                           here( __PRETTY_FUNCTION__, "the test stiffness does not couple normal and shear" ) );
}

// uniaxial stress with the same anisotropic stiffness: all stresses but sigma_11 vanish, tangent 1 / S_11
void testUniaxialStressAnisotropicLinear()
{
  const Matrix6d     C = anisotropicStiffness();
  LinearStressUpdate update{ C };

  Vector6d dStrain = Vector6d::Zero();
  dStrain( 0 )     = 2e-4;

  const auto condensation = computeReducedStress( update, dStrain, uniaxialStressPrescribed, uniaxialStressFree );

  const Vector6d stress = stressAt( update, dStrain );
  throwExceptionOnFailure( maxAbsAt( stress, uniaxialStressFree ) <= 1e-12 * maxAbs( stress ),
                           here( __PRETTY_FUNCTION__, "sigma_f does not vanish" ) );

  const double exact = 1. / C.inverse()( 0, 0 );
  throwExceptionOnFailure( std::abs( condensation.dStressP_dStrainP( 0, 0 ) - exact ) <= 1e-12 * exact,
                           here( __PRETTY_FUNCTION__, "tangent != 1 / S_11" ) );
  throwExceptionOnFailure( std::abs( stress( 0 ) - exact * 2e-4 ) <= 1e-12 * std::abs( stress( 0 ) ),
                           here( __PRETTY_FUNCTION__, "sigma_11 != tangent * eps_11" ) );

  // with the shear strains held zero (the iteration of master), the state would differ
  const auto lateralOnly = condense( C, uniaxialStressPrescribed, std::array< int, 2 >{ 1, 2 } );
  throwExceptionOnFailure( std::abs( lateralOnly.dStressP_dStrainP( 0, 0 ) - exact ) > 1e-3 * exact,
                           here( __PRETTY_FUNCTION__, "the test stiffness does not couple normal and shear" ) );
}

// index sets that are neither plane nor uniaxial stress: biaxial stress (sigma_12 free as well), and a set with a
// component in neither set, whose given strain is kept
void testOtherIndexSets()
{
  const Matrix6d     C = anisotropicStiffness();
  LinearStressUpdate update{ C };

  {
    constexpr std::array< int, 2 > p       = { 0, 1 };
    constexpr std::array< int, 4 > f       = { 2, 3, 4, 5 };
    Vector6d                       dStrain = Vector6d::Zero();
    dStrain( 0 )                           = 1e-4;
    dStrain( 1 )                           = 4e-5;

    const auto     condensation = computeReducedStress( update, dStrain, p, f );
    const Vector6d stress       = stressAt( update, dStrain );
    throwExceptionOnFailure( maxAbsAt( stress, f ) <= 1e-12 * maxAbs( stress ),
                             here( __PRETTY_FUNCTION__, "biaxial: sigma_f does not vanish" ) );
    throwExceptionOnFailure( checkRelative( condensation.dStressP_dStrainP, condensedFromCompliance( C, p ), 1e-12 ),
                             here( __PRETTY_FUNCTION__, "biaxial: condensed tangent != Schur complement" ) );
  }
  {
    // p = {11, 13}, f = {22, 33, 23}, gamma12 is given and kept
    constexpr std::array< int, 2 > p       = { 0, 4 };
    constexpr std::array< int, 3 > f       = { 1, 2, 5 };
    Vector6d                       dStrain = Vector6d::Zero();
    dStrain( 0 )                           = 1e-4;
    dStrain( 4 )                           = -2e-4;
    dStrain( 3 )                           = 3e-5;

    const auto     condensation = computeReducedStress( update, dStrain, p, f );
    const Vector6d stress       = stressAt( update, dStrain );
    throwExceptionOnFailure( dStrain( 3 ) == 3e-5, here( __PRETTY_FUNCTION__, "the given gamma12 was changed" ) );
    throwExceptionOnFailure( maxAbsAt( stress, f ) <= 1e-12 * maxAbs( stress ),
                             here( __PRETTY_FUNCTION__, "partial set: sigma_f does not vanish" ) );

    // the Schur complement of the 5x5 stiffness without gamma12
    constexpr std::array< int, 5 > pf = { 0, 4, 1, 2, 5 };
    Eigen::Matrix< double, 5, 5 >  C5;
    for ( int i = 0; i < 5; i++ )
      for ( int j = 0; j < 5; j++ )
        C5( i, j ) = C( pf[i], pf[j] );
    const Eigen::Matrix2d exact = C5.inverse().topLeftCorner< 2, 2 >().inverse();
    throwExceptionOnFailure( checkRelative( condensation.dStressP_dStrainP, exact, 1e-12 ),
                             here( __PRETTY_FUNCTION__, "partial set: condensed tangent != Schur complement" ) );
  }
}

// a nonlinear update, so that Newton iterates: converged sigma_f, last evaluation at the converged strain, and the
// condensed tangent and the sensitivity against central finite differences of the reduced stress state
void testNonlinearPlaneStress()
{
  const Matrix6d    C = anisotropicStiffness();
  Vector6d          lastStrain;
  CubicStressUpdate update{ C, 2e3, &lastStrain };

  auto solve = [&]( const Eigen::Vector3d& strainP, Vector6d& dStrain ) {
    dStrain.setZero();
    for ( int i = 0; i < 3; i++ )
      dStrain( planeStressPrescribed[i] ) = strainP( i );
    dStrain( 2 ) = -strainP( 0 ) - strainP( 1 );
    return computeReducedStress( update, dStrain, planeStressPrescribed, planeStressFree );
  };

  const Eigen::Vector3d strainP( 1.2e-2, -4e-3, 6e-3 );
  Vector6d              dStrain;
  const auto            condensation = solve( strainP, dStrain );

  const Vector6d stress = stressAt( update, dStrain );
  throwExceptionOnFailure( maxAbsAt( stress, planeStressFree ) <= 1e-12 * maxAbs( stress ),
                           here( __PRETTY_FUNCTION__, "sigma_f does not vanish" ) );
  throwExceptionOnFailure( condensation.nIterations > 2,
                           here( __PRETTY_FUNCTION__, "the nonlinear update did not need iterations" ) );
  throwExceptionOnFailure( lastStrain == dStrain,
                           here( __PRETTY_FUNCTION__, "the last evaluation is not at the converged strain" ) );

  Eigen::Matrix3d dStressP_dStrainP_FD;
  Eigen::Matrix3d dStrainF_dStrainP_FD;
  for ( int j = 0; j < 3; j++ ) {
    const double    h    = 1e-7;
    Eigen::Vector3d plus = strainP, minus = strainP;
    plus( j ) += h;
    minus( j ) -= h;
    Vector6d dStrainPlus, dStrainMinus;
    solve( plus, dStrainPlus );
    solve( minus, dStrainMinus );
    const Vector6d stressPlus = stressAt( update, dStrainPlus ), stressMinus = stressAt( update, dStrainMinus );
    for ( int i = 0; i < 3; i++ ) {
      dStressP_dStrainP_FD( i,
                            j ) = ( stressPlus( planeStressPrescribed[i] ) - stressMinus( planeStressPrescribed[i] ) ) /
                                  ( 2 * h );
      dStrainF_dStrainP_FD( i, j ) = ( dStrainPlus( planeStressFree[i] ) - dStrainMinus( planeStressFree[i] ) ) /
                                     ( 2 * h );
    }
  }
  throwExceptionOnFailure( checkRelative( condensation.dStressP_dStrainP, dStressP_dStrainP_FD, 1e-6 ),
                           here( __PRETTY_FUNCTION__, "condensed tangent != finite differences" ) );
  throwExceptionOnFailure( checkRelative( condensation.dStrainF_dStrainP, dStrainF_dStrainP_FD, 1e-6 ),
                           here( __PRETTY_FUNCTION__, "dStrainF_dStrainP != finite differences" ) );
}

// tiny and zero stresses: the tolerance is relative, a zero increment converges at once
void testScaleFreeTolerance()
{
  const Matrix6d     C = anisotropicStiffness();
  LinearStressUpdate update{ C };

  Vector6d dStrain = Vector6d::Zero();
  dStrain( 0 )     = 1e-15; // stresses of about 3e-5 Pa
  computeReducedStress( update, dStrain, planeStressPrescribed, planeStressFree );
  const Vector6d stress = stressAt( update, dStrain );
  throwExceptionOnFailure( maxAbsAt( stress, planeStressFree ) <= 1e-12 * maxAbs( stress ),
                           here( __PRETTY_FUNCTION__, "tiny stresses: sigma_f does not vanish relatively" ) );

  Vector6d   zero         = Vector6d::Zero();
  const auto condensation = computeReducedStress( update, zero, planeStressPrescribed, planeStressFree );
  throwExceptionOnFailure( condensation.nIterations == 1 && zero.isZero( 0 ),
                           here( __PRETTY_FUNCTION__, "a zero increment must converge with one evaluation" ) );
}

// a consistent singular C_ff: no shear stiffness in 13 and 23 and nothing driving them. The minimum norm solution keeps
// the initial guess of the undetermined shear strains, and the tangent is the plane stress tangent of the isotropic
// stiffness.
void testConsistentSingularCff()
{
  Matrix6d C = Marmot::ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( E, nu );
  C( 4, 4 )  = 0;
  C( 5, 5 )  = 0;
  LinearStressUpdate update{ C };

  Vector6d dStrain = Vector6d::Zero();
  dStrain( 0 )     = 1e-4;
  dStrain( 1 )     = -2e-5;
  dStrain( 3 )     = 3e-5;
  dStrain( 2 )     = -dStrain( 0 ) - dStrain( 1 );
  dStrain( 4 )     = 2e-5; // initial guesses of the undetermined shear strains
  dStrain( 5 )     = -1e-5;

  const auto     condensation = computeReducedStress( update, dStrain, planeStressPrescribed, planeStressFree );
  const Vector6d stress       = stressAt( update, dStrain );
  throwExceptionOnFailure( maxAbsAt( stress, planeStressFree ) <= 1e-12 * maxAbs( stress ),
                           here( __PRETTY_FUNCTION__, "sigma_f does not vanish" ) );
  throwExceptionOnFailure( dStrain( 4 ) == 2e-5 && dStrain( 5 ) == -1e-5,
                           here( __PRETTY_FUNCTION__, "the undetermined shear strains were changed" ) );
  throwExceptionOnFailure( condensation.CffCOD.rank() == 1,
                           here( __PRETTY_FUNCTION__, "C_ff is not recognized as singular" ) );

  const double    G = E / ( 2 * ( 1 + nu ) );
  Eigen::Matrix3d exact;
  exact << E / ( 1 - nu * nu ), nu * E / ( 1 - nu * nu ), 0, nu * E / ( 1 - nu * nu ), E / ( 1 - nu * nu ), 0, 0, 0, G;
  throwExceptionOnFailure( checkRelative( condensation.dStressP_dStrainP, exact, 1e-12 ),
                           here( __PRETTY_FUNCTION__, "condensed tangent != plane stress tangent" ) );
  throwExceptionOnFailure( condensation.dStrainF_dStrainP.bottomRows< 2 >().isZero( 1e-14 ),
                           here( __PRETTY_FUNCTION__, "the undetermined shear strains have a sensitivity" ) );
  throwExceptionOnFailure( std::abs( condensation.dStrainF_dStrainP( 0, 0 ) + nu / ( 1 - nu ) ) <= 1e-12,
                           here( __PRETTY_FUNCTION__, "deps33 / deps11 != -nu / (1 - nu)" ) );
}

// a singular C_ff (vanishing out-of-plane shear stiffness) throws StressUpdateFailed if the reduced problem is
// inconsistent: a stress-free component is driven (coupling or a stress without stiffness); a diverging iteration
// throws
void testFailures()
{
  auto throwsStressUpdateFailed = []( const std::function< void() >& f ) {
    try {
      f();
    }
    catch ( const Marmot::Utilities::StressUpdateFailed& ) {
      return true;
    }
    return false;
  };

  // no shear stiffness in 13 and 23, but coupling of 11 and 13: sigma_11 depends on the undetermined gamma13
  {
    Matrix6d C = Marmot::ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( E, nu );
    C( 4, 4 )  = 0;
    C( 5, 5 )  = 0;
    C( 0, 4 )  = 4e9;
    C( 4, 0 )  = 4e9;
    throwExceptionOnFailure( throwsStressUpdateFailed(
                               [&]() { condense( C, planeStressPrescribed, planeStressFree ); } ),
                             here( __PRETTY_FUNCTION__, "condense() with an inconsistent singular C_ff must throw" ) );

    LinearStressUpdate coupled{ C };
    throwExceptionOnFailure( throwsStressUpdateFailed( [&]() {
                               Vector6d dStrain = Vector6d::Zero();
                               dStrain( 0 )     = 1e-4;
                               computeReducedStress( coupled, dStrain, planeStressPrescribed, planeStressFree );
                             } ),
                             here( __PRETTY_FUNCTION__, "a coupled singular C_ff must throw" ) );
  }

  // no shear stiffness in 13 and 23, but a shear stress sigma_13 that no strain can remove
  {
    struct DrivenShear {
      Matrix6d C;
      void     operator()( const Vector6d& dStrain, Vector6d& stress, Matrix6d& dStress_dStrain ) const
      {
        stress = C * dStrain;
        stress( 4 ) += 1e6; // Pa
        dStress_dStrain = C;
      }
    };
    DrivenShear driven{ Marmot::ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( E, nu ) };
    driven.C( 4, 4 ) = 0;
    driven.C( 5, 5 ) = 0;
    throwExceptionOnFailure( throwsStressUpdateFailed( [&]() {
                               Vector6d dStrain = Vector6d::Zero();
                               dStrain( 0 )     = 1e-4;
                               computeReducedStress( driven, dStrain, planeStressPrescribed, planeStressFree );
                             } ),
                             here( __PRETTY_FUNCTION__, "a stress outside the range of a singular C_ff must throw" ) );
  }

  // a tangent 10 times too soft: every correction overshoots by a factor of 9, Newton diverges
  struct WrongTangent {
    Matrix6d C;
    void     operator()( const Vector6d& dStrain, Vector6d& stress, Matrix6d& dStress_dStrain ) const
    {
      stress          = C * dStrain;
      dStress_dStrain = 0.1 * C;
    }
  };
  WrongTangent diverging{ anisotropicStiffness() };
  throwExceptionOnFailure( throwsStressUpdateFailed( [&]() {
                             Vector6d dStrain = Vector6d::Zero();
                             dStrain( 0 )     = 1e-4;
                             computeReducedStress( diverging, dStrain, planeStressPrescribed, planeStressFree );
                           } ),
                           here( __PRETTY_FUNCTION__, "a diverging iteration must throw" ) );

  bool invalid = false;
  try {
    condense( anisotropicStiffness(), std::array< int, 2 >{ 0, 1 }, std::array< int, 2 >{ 1, 2 } );
  }
  catch ( const std::invalid_argument& ) {
    invalid = true;
  }
  throwExceptionOnFailure( invalid, here( __PRETTY_FUNCTION__, "overlapping index sets must be rejected" ) );

  // an index out of range is rejected before the stress update or the strain is accessed at it
  bool rejected = false, evaluated = false;
  try {
    Vector6d dStrain = Vector6d::Zero();
    computeReducedStress( [&]( const Vector6d&, Vector6d&, Matrix6d& ) { evaluated = true; },
                          dStrain,
                          std::array< int, 1 >{ 0 },
                          std::array< int, 1 >{ 6 } );
  }
  catch ( const std::invalid_argument& ) {
    rejected = true;
  }
  throwExceptionOnFailure( rejected && !evaluated,
                           here( __PRETTY_FUNCTION__, "an index out of range must be rejected before the update" ) );
}

// freeStrainSensitivity(): the free strains keep the free stresses zero when a further variable x changes; with a
// singular C_ff, a dStressF_dX in its range passes, one with a component along its null space throws
void testFreeStrainSensitivity()
{
  using Matrix32 = Eigen::Matrix< double, 3, 2 >;

  // regular C_ff: dStrainF_dX = -C_ff^-1 dStressF_dX
  {
    const auto condensation = condense( anisotropicStiffness(), planeStressPrescribed, planeStressFree );
    Matrix32   dStressF_dX;
    dStressF_dX << 1e9, -2e8, 3e8, 5e8, -4e8, 1e8;
    const Matrix32 dStrainF_dX = freeStrainSensitivity( condensation, dStressF_dX );
    throwExceptionOnFailure( checkRelative( condensation.Cff * dStrainF_dX, -dStressF_dX, 1e-12 ),
                             here( __PRETTY_FUNCTION__, "regular: C_ff dStrainF_dX != -dStressF_dX" ) );
  }

  // singular C_ff: no shear stiffness in 13 and 23
  Matrix6d C              = Marmot::ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( E, nu );
  C( 4, 4 )               = 0;
  C( 5, 5 )               = 0;
  const auto condensation = condense( C, planeStressPrescribed, planeStressFree );

  Matrix32 inRange           = Matrix32::Zero();
  inRange( 0, 0 )            = 1e9; // only sigma_33 depends on x
  inRange( 0, 1 )            = -3e8;
  const Matrix32 dStrainF_dX = freeStrainSensitivity( condensation, inRange );
  throwExceptionOnFailure( checkRelative( condensation.Cff * dStrainF_dX, -inRange, 1e-12 ),
                           here( __PRETTY_FUNCTION__, "singular: C_ff dStrainF_dX != -dStressF_dX" ) );
  throwExceptionOnFailure( maxAbs( dStrainF_dX.bottomRows< 2 >() ) <= 1e-14 * maxAbs( dStrainF_dX ),
                           here( __PRETTY_FUNCTION__, "singular: the undetermined shear strains have a sensitivity" ) );

  Matrix32 outsideRange = inRange;
  outsideRange( 1, 1 )  = 1e6; // sigma_13 depends on x, but no shear strain can compensate
  bool thrown           = false;
  try {
    freeStrainSensitivity( condensation, outsideRange );
  }
  catch ( const Marmot::Utilities::StressUpdateFailed& ) {
    thrown = true;
  }
  throwExceptionOnFailure( thrown, here( __PRETTY_FUNCTION__, "singular: dStressF_dX outside the range must throw" ) );
}

// MaterialHypoElastic::computePlaneStress / computeUniaxialStress with LinearElastic, orthotropic with axes rotated
// out of the plane (normal - out-of-plane shear coupling), in Pa
void testMaterialHypoElasticAnisotropic()
{
  // E1, E2, E3, nu12, nu23, nu13, G12, G23, G13, direction 1, direction 2
  const double                     s     = std::sqrt( 0.5 );
  const std::array< double, 15 >   props = { 40e9, 10e9, 12e9, 0.25, 0.3, 0.2, 5e9, 4e9, 6e9, s, 0, s, 0, 1, 0 };
  Marmot::Materials::LinearElastic material( props.data(), static_cast< int >( props.size() ), 1 );

  const Marmot::Materials::MaterialHypoElastic::timeInfo timeInfo{ 0., 1. };

  Matrix6d C;
  {
    Marmot::Materials::MaterialHypoElastic::state3D state;
    const Marmot::Materials::MaterialHypoElastic&   base = material; // computeStress is protected in LinearElastic
    base.computeStress( state, C, Vector6d::Zero(), timeInfo );
  }
  throwExceptionOnFailure( std::abs( C( 0, 4 ) ) > 1e-3 * C( 0, 0 ),
                           here( __PRETTY_FUNCTION__, "the rotated material does not couple 11 and 13" ) );

  // plane stress
  {
    Marmot::Materials::MaterialHypoElastic::state2D state;
    Marmot::Math::Matrix3d                          tangent;
    const Marmot::Math::Vector3d                    dStrain( 1e-4, -2e-5, 3e-5 );
    material.computePlaneStress( state, tangent, dStrain, timeInfo );

    const Eigen::Matrix3d exact = condensedFromCompliance( C, planeStressPrescribed );
    throwExceptionOnFailure( checkRelative( tangent, exact, 1e-11 ),
                             here( __PRETTY_FUNCTION__, "plane stress tangent != Schur complement" ) );
    throwExceptionOnFailure( checkRelative( state.stress, exact * dStrain, 1e-11 ),
                             here( __PRETTY_FUNCTION__, "plane stress != Schur complement * strain" ) );
  }
  // uniaxial stress
  {
    Marmot::Materials::MaterialHypoElastic::state1D state;
    double                                          tangent;
    material.computeUniaxialStress( state, tangent, 1e-4, timeInfo );

    const double exact = 1. / C.inverse()( 0, 0 );
    throwExceptionOnFailure( std::abs( tangent - exact ) <= 1e-11 * exact,
                             here( __PRETTY_FUNCTION__, "uniaxial tangent != 1 / S_11" ) );
    throwExceptionOnFailure( std::abs( state.stress - exact * 1e-4 ) <= 1e-11 * exact * 1e-4,
                             here( __PRETTY_FUNCTION__, "uniaxial stress != tangent * strain" ) );
  }
}

// MaterialHypoElastic::computePlaneStress / computeUniaxialStress with VonMises in Pa along a plastic path: the
// condensed tangent against central finite differences of the reduced stress update
void testMaterialHypoElasticPlastic()
{
  // E, nu, fy, HLin, deltaYieldStress, delta, density
  const std::array< double, 7 >    props = { 210e9, 0.3, 235e6, 2e9, 100e6, 50, 7850 };
  Marmot::Materials::VonMisesModel material( props.data(), static_cast< int >( props.size() ), 1 );
  const int                        nStateVars = material.getNumberOfRequiredStateVars();
  const Marmot::Materials::MaterialHypoElastic::timeInfo timeInfo{ 0., 1. };

  // plane stress
  {
    std::vector< double >                           stateVars( nStateVars, 0.0 );
    Marmot::Materials::MaterialHypoElastic::state2D state;
    state.stateVars = stateVars.data();
    Marmot::Math::Matrix3d       tangent;
    const Marmot::Math::Vector3d dStrain( 8e-4, -2e-4, 6e-4 );
    for ( int n = 0; n < 5; n++ ) // into plasticity
      material.computePlaneStress( state, tangent, dStrain, timeInfo );
    throwExceptionOnFailure( stateVars[0] > 0, here( __PRETTY_FUNCTION__, "the path did not yield" ) );

    auto update = [&]( const Marmot::Math::Vector3d& dE ) {
      std::vector< double >                           sv = stateVars;
      Marmot::Materials::MaterialHypoElastic::state2D st = state;
      st.stateVars                                       = sv.data();
      Marmot::Math::Matrix3d t;
      material.computePlaneStress( st, t, dE, timeInfo );
      return std::make_pair( Marmot::Math::Vector3d( st.stress ), t );
    };
    const auto [stress, analytic] = update( dStrain );
    Eigen::Matrix3d fd;
    for ( int j = 0; j < 3; j++ ) {
      const double           h  = 1e-9;
      Marmot::Math::Vector3d dP = dStrain, dM = dStrain;
      dP( j ) += h;
      dM( j ) -= h;
      fd.col( j ) = ( update( dP ).first - update( dM ).first ) / ( 2 * h );
    }
    throwExceptionOnFailure( checkRelative( analytic, fd, 1e-5 ),
                             here( __PRETTY_FUNCTION__, "plastic plane stress tangent != finite differences" ) );
  }
  // uniaxial stress
  {
    std::vector< double >                           stateVars( nStateVars, 0.0 );
    Marmot::Materials::MaterialHypoElastic::state1D state;
    state.stateVars = stateVars.data();
    double tangent;
    for ( int n = 0; n < 5; n++ )
      material.computeUniaxialStress( state, tangent, 8e-4, timeInfo );
    throwExceptionOnFailure( stateVars[0] > 0, here( __PRETTY_FUNCTION__, "the uniaxial path did not yield" ) );

    auto update = [&]( double dE ) {
      std::vector< double >                           sv = stateVars;
      Marmot::Materials::MaterialHypoElastic::state1D st = state;
      st.stateVars                                       = sv.data();
      double t;
      material.computeUniaxialStress( st, t, dE, timeInfo );
      return std::make_pair( st.stress, t );
    };
    const double h        = 1e-9;
    const double analytic = update( 8e-4 ).second;
    const double fd       = ( update( 8e-4 + h ).first - update( 8e-4 - h ).first ) / ( 2 * h );
    throwExceptionOnFailure( std::abs( analytic - fd ) <= 1e-5 * std::abs( fd ),
                             here( __PRETTY_FUNCTION__, "plastic uniaxial tangent != finite differences" ) );
  }
}

// regression: RankineSofteningPlasticity (Pa) in plane stress through MaterialHypoElastic::computePlaneStress, driven
// to full softening. The fully softened crack has no out-of-plane shear stiffness, so C_ff is singular but
// consistent; the plane stress update must not throw.
void testRankineFullSofteningPlaneStress()
{
  // E, nu, fy0, HLin, deltaYieldStress, delta, density
  const std::array< double, 7 >                   props = { 30e9, 0.2, 3e6, -2e9, 0., 0., 2400. };
  Marmot::Materials::RankineSofteningPlasticity   material( props.data(), static_cast< int >( props.size() ), 1 );
  std::vector< double >                           stateVars( material.getNumberOfRequiredStateVars(), 0.0 );
  Marmot::Materials::MaterialHypoElastic::state2D state;
  state.stateVars = stateVars.data();
  Marmot::Math::Matrix3d                                 tangent;
  const Marmot::Materials::MaterialHypoElastic::timeInfo timeInfo{ 0., 1. };
  const Marmot::Math::Vector3d                           dStrain( 1e-4, 2e-5, 3e-5 );

  double peak = 0;
  for ( int n = 0; n < 30; n++ ) {
    try {
      material.computePlaneStress( state, tangent, dStrain, timeInfo );
    }
    catch ( const Marmot::Utilities::StressUpdateFailed& e ) {
      throw std::runtime_error(
        here( __PRETTY_FUNCTION__, "increment " + std::to_string( n + 1 ) + " threw: " + e.what() ) );
    }
    peak = std::max( peak, state.stress.cwiseAbs().maxCoeff() );
    throwExceptionOnFailure( tangent.allFinite(), here( __PRETTY_FUNCTION__, "the tangent is not finite" ) );
  }
  throwExceptionOnFailure( peak >= 1e6 && state.stress.cwiseAbs().maxCoeff() <= 1e-9 * peak,
                           here( __PRETTY_FUNCTION__, "the path did not reach full softening" ) );
}

int main()
{
  const std::vector< std::function< void() > > tests = {
    testPlaneStressAnisotropicLinear,
    testUniaxialStressAnisotropicLinear,
    testOtherIndexSets,
    testNonlinearPlaneStress,
    testScaleFreeTolerance,
    testConsistentSingularCff,
    testFailures,
    testFreeStrainSensitivity,
    testMaterialHypoElasticAnisotropic,
    testMaterialHypoElasticPlastic,
    testRankineFullSofteningPlaneStress,
  };
  executeTestsAndCollectExceptions( tests );
  return 0;
}
