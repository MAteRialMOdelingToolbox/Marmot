#include "Marmot/HaighWestergaard.h"
#include "Marmot/MarmotTesting.h"

using namespace Marmot::Testing;
using namespace Marmot::ContinuumMechanics::HaighWestergaard;
using namespace Marmot::ContinuumMechanics::VoigtNotation::Invariants;
using namespace autodiff;

void testHaighWestergaardDouble()
{
  Eigen::Matrix< double, 6, 1 >         stress, smean, sdev;
  HaighWestergaardCoordinates< double > hwoutput, hwref;

  // first test: hydrostatic stress state
  stress << -10., -10., -10., 0., 0., 0.;

  hwoutput    = haighWestergaard( stress );
  hwref.xi    = -10. * sqrt( 3. );
  hwref.rho   = 0.;
  hwref.theta = 0.;
  throwExceptionOnFailure( checkIfEqual( hwoutput.xi, hwref.xi, 1e-13 ),
                           MakeString() << __PRETTY_FUNCTION__ << " test 1 failed: error in xi" );
  throwExceptionOnFailure( checkIfEqual( hwoutput.rho, hwref.rho ),
                           MakeString() << __PRETTY_FUNCTION__ << " test 1 failed: error in rho" );
  throwExceptionOnFailure( checkIfEqual( hwoutput.theta, hwref.theta ),
                           MakeString() << __PRETTY_FUNCTION__ << " test 1 failed: error in theta" );

  // second test: compression meridian
  stress << -10., 0., 0., 0, 0., 0.;
  hwoutput    = haighWestergaard( stress );
  hwref.xi    = -10. / 3. * sqrt( 3. );
  hwref.rho   = 10. * sqrt( 2. / 3. );
  hwref.theta = Marmot::Constants::Pi / 3.;
  throwExceptionOnFailure( checkIfEqual( hwoutput.xi, hwref.xi, 1e-13 ),
                           MakeString() << __PRETTY_FUNCTION__ << " test 2 failed: error in xi" );
  throwExceptionOnFailure( checkIfEqual( hwoutput.rho, hwref.rho, 1e-13 ),
                           MakeString() << __PRETTY_FUNCTION__ << " test 2 failed: error in rho" );
  throwExceptionOnFailure( checkIfEqual( hwoutput.theta, hwref.theta ),
                           MakeString() << __PRETTY_FUNCTION__ << " test 2 failed: error in theta" );

  // third test: pure shear
  stress << -10., 10., 0., 0, 0., 0.;
  hwoutput    = haighWestergaard( stress );
  hwref.xi    = 0.;
  hwref.rho   = sqrt( 2. ) * 10.;
  hwref.theta = Marmot::Constants::Pi / 6.;
  throwExceptionOnFailure( checkIfEqual( hwoutput.xi, hwref.xi, 1e-13 ),
                           MakeString() << __PRETTY_FUNCTION__ << " test 3 failed: error in xi" );
  throwExceptionOnFailure( checkIfEqual( hwoutput.rho, hwref.rho, 1e-13 ),
                           MakeString() << __PRETTY_FUNCTION__ << " test 3 failed: error in rho" );
  throwExceptionOnFailure( checkIfEqual( hwoutput.theta, hwref.theta ),
                           MakeString() << __PRETTY_FUNCTION__ << " test 3 failed: error in theta" );

  // fourth test: general stress state
  stress << 10., -12., -23., 1.3, -6.8, 10.3;
  hwoutput = haighWestergaard( stress );
  double p = stress.head< 3 >().sum() / 3.;
  smean.setZero();
  smean.head< 3 >().setConstant( p );
  sdev = stress - smean;
  sdev.tail< 3 >() *= sqrt( 2. );
  const double J2_       = J2( stress );
  const double J3_       = J3( stress );
  double       cos3theta = 3. * sqrt( 3. ) / 2. * J3_ / pow( J2_, 3. / 2. );

  hwref.xi    = sqrt( 3. ) * p;
  hwref.rho   = sdev.norm();
  hwref.theta = 1. / 3. * acos( cos3theta );

  throwExceptionOnFailure( checkIfEqual( hwoutput.xi, hwref.xi, 1e-13 ),
                           MakeString() << __PRETTY_FUNCTION__ << " test 4  failed: error in xi" );
  throwExceptionOnFailure( checkIfEqual( hwoutput.rho, hwref.rho ),
                           MakeString() << __PRETTY_FUNCTION__ << " test 4  failed: error in rho" );
  throwExceptionOnFailure( checkIfEqual( hwoutput.theta, hwref.theta ),
                           MakeString() << __PRETTY_FUNCTION__ << " test 4  failed: error in theta" );
}

void testHaighWestergaardDual()
{
  Eigen::Matrix< dual, 6, 1 >         stress, smean, sdev;
  HaighWestergaardCoordinates< dual > hwoutput, hwref;

  // first test: hydrostatic stress state
  stress << dual( -10. ), dual( -10. ), dual( -10. ), dual( 0. ), dual( 0. ), dual( 0. );
  stress( 0 ).grad = 1.0;
  stress( 1 ).grad = 1.0;
  stress( 2 ).grad = 1.0;

  hwoutput    = haighWestergaard( stress );
  hwref.xi    = stress( 0 ) * sqrt( 3. );
  hwref.rho   = 0.;
  hwref.theta = 0.;

  throwExceptionOnFailure( checkIfEqual( hwoutput.xi, hwref.xi, 1e-13 ),
                           MakeString() << __PRETTY_FUNCTION__ << " test 1  failed: error in xi" );
  // throwExceptionOnFailure( checkIfEqual( hwoutput.rho, hwref.rho ),
  //                          MakeString() << __PRETTY_FUNCTION__ << " test 1  failed: error in rho" );
  throwExceptionOnFailure( checkIfEqual( hwoutput.theta, hwref.theta ),
                           MakeString() << __PRETTY_FUNCTION__ << " test 1  failed: error in theta" );

  // second test: compression meridian
  stress << dual( -10. ), dual( 0. ), dual( 0. ), dual( 0. ), dual( 0. ), dual( 0. );
  hwoutput    = haighWestergaard( stress );
  hwref.xi    = -10. / 3. * sqrt( 3. );
  hwref.rho   = 10. * sqrt( 2. / 3. );
  hwref.theta = Marmot::Constants::Pi / 3.;
  throwExceptionOnFailure( checkIfEqual( hwoutput.xi, hwref.xi, 1e-13 ),
                           MakeString() << __PRETTY_FUNCTION__ << " test 2 failed: error in xi" );
  throwExceptionOnFailure( checkIfEqual( hwoutput.rho, hwref.rho, 1e-13 ),
                           MakeString() << __PRETTY_FUNCTION__ << " test 2 failed: error in rho" );
  throwExceptionOnFailure( checkIfEqual( hwoutput.theta, hwref.theta ),
                           MakeString() << __PRETTY_FUNCTION__ << " test 2 failed: error in theta" );

  // third test: pure shear
  stress << dual( -10. ), dual( 10. ), dual( 0. ), dual( 0. ), dual( 0. ), dual( 0. );
  hwoutput    = haighWestergaard( stress );
  hwref.xi    = 0.;
  hwref.rho   = sqrt( 2. ) * 10.;
  hwref.theta = Marmot::Constants::Pi / 6.;
  throwExceptionOnFailure( checkIfEqual( hwoutput.xi, hwref.xi, 1e-13 ),
                           MakeString() << __PRETTY_FUNCTION__ << " test 3 failed: error in xi" );
  throwExceptionOnFailure( checkIfEqual( hwoutput.rho, hwref.rho, 1e-13 ),
                           MakeString() << __PRETTY_FUNCTION__ << " test 3 failed: error in rho" );
  throwExceptionOnFailure( checkIfEqual( hwoutput.theta, hwref.theta ),
                           MakeString() << __PRETTY_FUNCTION__ << " test 3 failed: error in theta" );

  // fourth test: general stress state
  stress << dual( 10. ), dual( -12. ), dual( -23. ), dual( 1.3 ), dual( -6.8 ), dual( 10.3 );
  stress( 1 ).grad = 1.0; // set imaginary part of the seconde component of the stress tensor to 1.0 (to calculate
                          // derivative with regard to that component)
  hwoutput = haighWestergaard( stress );
  dual p   = stress.head< 3 >().sum() / 3.;
  smean.setZero();
  smean.head< 3 >().setConstant( p );
  sdev = stress - smean;
  sdev.tail< 3 >() *= sqrt( 2. );
  const dual J2_       = J2( stress );
  const dual J3_       = J3( stress );
  dual       cos3theta = 3. * sqrt( 3. ) / 2. * J3_ / pow( J2_, 3. / 2. );

  hwref.xi    = sqrt( 3. ) * p;
  hwref.rho   = sdev.norm();
  hwref.theta = 1. / 3. * acos( cos3theta );

  throwExceptionOnFailure( checkIfEqual( hwoutput.xi, hwref.xi, 1e-13 ),
                           MakeString() << __PRETTY_FUNCTION__ << " test 4  failed: error in xi" );
  throwExceptionOnFailure( checkIfEqual( hwoutput.rho, hwref.rho ),
                           MakeString() << __PRETTY_FUNCTION__ << " test 4  failed: error in rho" );
  throwExceptionOnFailure( checkIfEqual( hwoutput.theta, hwref.theta ),
                           MakeString() << __PRETTY_FUNCTION__ << " test 4  failed: error in theta" );
}

void testHaighWestergaardComplexDouble()
{
  Eigen::Matrix< std::complex< double >, 6, 1 >         stress, smean, sdev;
  HaighWestergaardCoordinates< std::complex< double > > hwoutput, hwref;

  // general stress state
  stress << std::complex< double >( 10., 0. ), std::complex< double >( -12., 0. ), std::complex< double >( -23., 0. ),
    std::complex< double >( 1.3, 0. ), std::complex< double >( -6.8, 0. ), std::complex< double >( 10.3, 0. );

  hwoutput                 = haighWestergaard( stress );
  std::complex< double > p = stress.head< 3 >().sum() / 3.;
  smean.setZero();
  smean.head< 3 >().setConstant( p );
  sdev = stress - smean;
  sdev.tail< 3 >() *= sqrt( 2. );
  const std::complex< double > J2_       = J2( stress );
  const std::complex< double > J3_       = J3( stress );
  const std::complex< double > cos3theta = 3. * sqrt( 3. ) / 2. * J3_ / pow( J2_, 3. / 2. );

  hwref.xi    = sqrt( 3. ) * p;
  hwref.rho   = sdev.norm();
  hwref.theta = 1. / 3. * acos( cos3theta );

  throwExceptionOnFailure( checkIfEqual( hwoutput.xi, hwref.xi, 1e-13 ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed: error in xi" );
  throwExceptionOnFailure( checkIfEqual( hwoutput.rho, hwref.rho ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed: error in rho" );
  throwExceptionOnFailure( checkIfEqual( hwoutput.theta, hwref.theta ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed: error in theta" );
}

void testHaighWestergaardFromStrain()
{
  Eigen::Matrix< double, 6, 1 > strain, epsmean, epsdev;
  strain << 10e-3, -1.2e-3, 7e-4, 1.3e-3, 6.8e-3, 10.3e-3;
  HaighWestergaardCoordinates< double > hwoutput, hwref;
  hwoutput = haighWestergaardFromStrain( strain );
  double p = strain.head< 3 >().sum() / 3.;
  epsmean.setZero();
  epsmean.head< 3 >().setConstant( p );
  epsdev = strain - epsmean;
  epsdev.tail< 3 >() /= sqrt( 2. );
  const double J2_       = J2Strain( strain );
  const double J3_       = J3Strain( strain );
  double       cos3theta = 3. * sqrt( 3. ) / 2. * J3_ / pow( J2_, 3. / 2. );

  hwref.xi    = sqrt( 3. ) * p;
  hwref.rho   = epsdev.norm();
  hwref.theta = 1. / 3. * acos( cos3theta );

  throwExceptionOnFailure( checkIfEqual( hwoutput.xi, hwref.xi, 1e-13 ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed: error in xi" );
  throwExceptionOnFailure( checkIfEqual( hwoutput.rho, hwref.rho ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed: error in rho" );
  throwExceptionOnFailure( checkIfEqual( hwoutput.theta, hwref.theta ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed: error in theta" );
}

int main()
{
  auto tests = std::vector< std::function< void() > >{
    testHaighWestergaardDouble,
    testHaighWestergaardDual,
    testHaighWestergaardComplexDouble,
    testHaighWestergaardFromStrain,
  };

  executeTestsAndCollectExceptions( tests );
}
