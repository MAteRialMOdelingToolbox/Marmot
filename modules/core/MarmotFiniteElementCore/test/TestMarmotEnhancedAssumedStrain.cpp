#include "Marmot/MarmotEnhancedAssumedStrain.h"
#include "Marmot/MarmotTesting.h"
#include <Eigen/Dense>

using namespace Marmot;
using namespace Marmot::Testing;
using namespace Marmot::FiniteElement::EAS;

// ---------------------------------------------------------------------------------------------
// F(J): the EAS transformation matrix built from the element Jacobian.
// ---------------------------------------------------------------------------------------------

void testFReducesToIdentityForIdentityJacobian2D()
{
  const Eigen::MatrixXd F2D = F( Eigen::Matrix2d::Identity() );
  throwExceptionOnFailure( checkIfEqual( F2D, Eigen::MatrixXd( Eigen::Matrix3d::Identity() ), 1e-12 ),
                           "F(Identity) must reduce to the identity for a 2D Jacobian." );
}

void testFReducesToIdentityForIdentityJacobian3D()
{
  const Eigen::MatrixXd F3D = F( Eigen::Matrix3d::Identity() );
  throwExceptionOnFailure( checkIfEqual( F3D, Eigen::MatrixXd( Eigen::Matrix< double, 6, 6 >::Identity() ), 1e-12 ),
                           "F(Identity) must reduce to the identity for a 3D Jacobian." );
}

void testFMatchesHandDerivedValuesForDiagonalJacobian2D()
{
  Eigen::Matrix2d J;
  J << 2.0, 0.0, 0.0, 3.0;

  const Eigen::MatrixXd F2D = F( J );

  Eigen::Matrix3d expected;
  // clang-format off
  expected << 4.0, 0.0, 0.0,
              0.0, 9.0, 0.0,
              0.0, 0.0, 6.0;
  // clang-format on

  throwExceptionOnFailure( checkIfEqual( F2D, Eigen::MatrixXd( expected ), 1e-12 ),
                           "F() does not match the hand-derived value for a diagonal 2D Jacobian." );
}

void testFThrowsForUnsupportedDimension()
{
  const Eigen::MatrixXd J = Eigen::MatrixXd::Identity( 4, 4 );

  bool threw = false;
  try {
    F( J );
  }
  catch ( const std::invalid_argument& ) {
    threw = true;
  }
  throwExceptionOnFailure( threw, "F() must throw for an unsupported (non-2D/3D) Jacobian." );
}

// ---------------------------------------------------------------------------------------------
// EASInterpolation(): each variant's M-matrix, checked against hand-computed values at a fixed
// natural coordinate.
// ---------------------------------------------------------------------------------------------

void testEASInterpolationDeBorstEAS2()
{
  Eigen::Vector2d       xi( 2.0, 3.0 );
  const Eigen::MatrixXd E = EASInterpolation( DeBorstEAS2, xi );

  Eigen::MatrixXd expected( 3, 2 );
  // clang-format off
  expected << 3.0, 0.0,
              0.0, 2.0,
              0.0, 0.0;
  // clang-format on
  throwExceptionOnFailure( checkIfEqual( E, expected, 1e-12 ), "DeBorstEAS2 interpolation matrix mismatch." );
}

void testEASInterpolationEAS3()
{
  Eigen::Vector3d       xi( 2.0, 3.0, 5.0 );
  const Eigen::MatrixXd E = EASInterpolation( EAS3, xi );

  Eigen::MatrixXd expected = Eigen::MatrixXd::Zero( 6, 3 );
  expected( 0, 0 )         = 2.0;
  expected( 1, 1 )         = 3.0;
  expected( 2, 2 )         = 5.0;

  throwExceptionOnFailure( checkIfEqual( E, expected, 1e-12 ), "EAS3 interpolation matrix mismatch." );
}

void testEASInterpolationDeBorstEAS9()
{
  Eigen::Vector3d       xi( 2.0, 3.0, 5.0 );
  const Eigen::MatrixXd E = EASInterpolation( DeBorstEAS9, xi );

  Eigen::MatrixXd expected = Eigen::MatrixXd::Zero( 6, 9 );
  expected( 0, 0 )         = 2.0;
  expected( 1, 1 )         = 3.0;
  expected( 2, 2 )         = 5.0;
  expected( 0, 3 )         = 6.0;  // xi0*xi1
  expected( 0, 4 )         = 10.0; // xi0*xi2
  expected( 1, 5 )         = 6.0;  // xi1*xi0
  expected( 1, 6 )         = 15.0; // xi1*xi2
  expected( 2, 7 )         = 10.0; // xi2*xi0
  expected( 2, 8 )         = 15.0; // xi2*xi1

  throwExceptionOnFailure( checkIfEqual( E, expected, 1e-12 ), "DeBorstEAS9 interpolation matrix mismatch." );
}

void testEASInterpolationDeBorstEAS6b()
{
  Eigen::Vector3d       xi( 2.0, 3.0, 5.0 );
  const Eigen::MatrixXd E = EASInterpolation( DeBorstEAS6b, xi );

  Eigen::MatrixXd expected = Eigen::MatrixXd::Zero( 6, 6 );
  expected( 0, 0 )         = 2.0;
  expected( 1, 1 )         = 3.0;
  expected( 2, 2 )         = 5.0;
  expected( 0, 3 )         = 6.0;  // xi1*xi0
  expected( 0, 4 )         = 10.0; // xi2*xi0
  expected( 1, 3 )         = 6.0;  // xi0*xi1
  expected( 1, 5 )         = 15.0; // xi2*xi1
  expected( 2, 4 )         = 10.0; // xi0*xi2
  expected( 2, 5 )         = 15.0; // xi1*xi2

  throwExceptionOnFailure( checkIfEqual( E, expected, 1e-12 ), "DeBorstEAS6b interpolation matrix mismatch." );
}

void testEASInterpolationSimoRifaiEAS5()
{
  Eigen::Vector2d       xi( 2.0, 3.0 );
  const Eigen::MatrixXd E = EASInterpolation( SimoRifaiEAS5, xi );

  Eigen::MatrixXd expected( 3, 5 );
  // clang-format off
  expected << 2.0, 0.0, 0.0, 0.0,  6.0,
              0.0, 3.0, 0.0, 0.0, -6.0,
              0.0, 0.0, 2.0, 3.0, -5.0; // xi0^2 - xi1^2 = 4 - 9 = -5
  // clang-format on
  throwExceptionOnFailure( checkIfEqual( E, expected, 1e-12 ), "SimoRifaiEAS5 interpolation matrix mismatch." );
}

void testEASInterpolationSimoRifaiEAS4()
{
  Eigen::Vector2d       xi( 2.0, 3.0 );
  const Eigen::MatrixXd E = EASInterpolation( SimoRifaiEAS4, xi );

  Eigen::MatrixXd expected( 3, 4 );
  // clang-format off
  expected << 2.0, 0.0, 0.0, 0.0,
              0.0, 3.0, 0.0, 0.0,
              0.0, 0.0, 2.0, 3.0;
  // clang-format on
  throwExceptionOnFailure( checkIfEqual( E, expected, 1e-12 ), "SimoRifaiEAS4 interpolation matrix mismatch." );
}

void testEASInterpolationThrowsForUnimplementedType()
{
  // DeBorstEAS2_P2 is declared in the EASType enum but has no case in EASInterpolation()'s switch.
  Eigen::Vector2d xi( 1.0, 1.0 );

  bool threw = false;
  try {
    EASInterpolation( DeBorstEAS2_P2, xi );
  }
  catch ( const std::invalid_argument& ) {
    threw = true;
  }
  throwExceptionOnFailure( threw, "EASInterpolation() must throw for an unimplemented EAS type." );
}

int main()
{
  auto tests = std::vector< std::function< void() > >{
    testFReducesToIdentityForIdentityJacobian2D,
    testFReducesToIdentityForIdentityJacobian3D,
    testFMatchesHandDerivedValuesForDiagonalJacobian2D,
    testFThrowsForUnsupportedDimension,
    testEASInterpolationDeBorstEAS2,
    testEASInterpolationEAS3,
    testEASInterpolationDeBorstEAS9,
    testEASInterpolationDeBorstEAS6b,
    testEASInterpolationSimoRifaiEAS5,
    testEASInterpolationSimoRifaiEAS4,
    testEASInterpolationThrowsForUnimplementedType,
  };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
