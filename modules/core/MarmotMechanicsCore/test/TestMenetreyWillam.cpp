#include "Marmot/HaighWestergaard.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MenetreyWillam.h"
#include "autodiff/forward/dual.hpp"
#include <iomanip>

using namespace Marmot::Testing;
using namespace Marmot::ContinuumMechanics::CommonConstitutiveModels;
using namespace Marmot::ContinuumMechanics::HaighWestergaard;

void testYieldFunctions()
{
  const double                                ft = 2.1;
  const double                                fc = 22.;
  const HaighWestergaardCoordinates< double > hw = { -1.0, 10., 0.4 };

  auto mwVM = MenetreyWillam( ft, MenetreyWillam::MenetreyWillamType::Mises );
  throwExceptionOnFailure( checkIfEqual( mwVM.yieldFunction( hw ), 4.8321184351980424 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " failed: Von Mises yield criterion doen't yield the right result" );

  auto mwR = MenetreyWillam( ft, MenetreyWillam::MenetreyWillamType::Rankine );
  throwExceptionOnFailure( checkIfEqual( mwR.yieldFunction( hw ), 2.2381914955371371 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " failed: Rankine yield criterion doen't yield the right result" );

  auto mwDP = MenetreyWillam( ft, MenetreyWillam::MenetreyWillamType::DruckerPrager, fc );
  throwExceptionOnFailure( checkIfEqual( mwDP.yieldFunction( hw ), 1.5483064286295773 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " failed: Drucker-Prager yield criterion doen't yield the right result" );

  auto mwMC = MenetreyWillam( ft, MenetreyWillam::MenetreyWillamType::MohrCoulomb, fc );
  throwExceptionOnFailure( checkIfEqual( mwMC.yieldFunction( hw ), 2.5218527589353812 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " failed: Mohr-Coulomb yield criterion doen't yield the right result" );
}

void testSmoothingWithFillet()
{
  const double                                ft     = 2.1;
  const double                                fc     = 22.;
  const HaighWestergaardCoordinates< double > hw     = { 2.0, 0., 0. };
  double                                      varEps = 0.2;

  auto mwDP = MenetreyWillam( ft, MenetreyWillam::MenetreyWillamType::DruckerPrager, fc );
  throwExceptionOnFailure( checkIfEqual( mwDP.yieldFunction( hw, varEps ), 0.4922077922077921 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " failed: Smoothing with fillet doesn't work for Drucker-Prager" );

  auto mwMC = MenetreyWillam( ft, MenetreyWillam::MenetreyWillamType::MohrCoulomb, fc );
  throwExceptionOnFailure( checkIfEqual( mwMC.yieldFunction( hw, varEps ), 0. - 0.3026289888799327 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " failed: Smoothing with fillet doesn't work for Mohr-Coulomb" );
}

void testInlineFunctions()
{
  const double ft  = 2.1;
  const double fc  = 22.;
  const double c   = MenetreyWillam::c( fc, ft );
  const double phi = MenetreyWillam::phi( fc, ft );
  const double ft_ = MenetreyWillam::ft( c, phi );
  const double fc_ = MenetreyWillam::fc( c, phi );

  throwExceptionOnFailure( checkIfEqual( ft, ft_, 1e-13 ) && checkIfEqual( fc, fc_, 1e-13 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " failed: Bug in inline conversion functions (ft, fc, phi, c)" );

  auto mw = MenetreyWillam( ft, MenetreyWillam::MenetreyWillamType::MohrCoulomb, fc );
  const HaighWestergaardCoordinates< double > hw = { -1.0, 10., 0.4 };
  throwExceptionOnFailure( checkIfEqual( mw.yieldFunction( hw ), 2.5218527589353812 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " failed: Mohr-Coulomb yield criterion doen't yield the right result" );
}

void testDerivatives()
{
  const double                                ft = 2.1;
  const double                                fc = 22.;
  const HaighWestergaardCoordinates< double > hw = { -1.0, 10., 0.4 };

  auto mwMC  = MenetreyWillam( ft, MenetreyWillam::MenetreyWillamType::MohrCoulomb, fc );
  auto dYdHW = mwMC.dYieldFunction_dHaighWestergaard( hw );

  autodiff::dual                                xi( hw.xi ), rho( hw.rho ), theta( hw.theta );
  HaighWestergaardCoordinates< autodiff::dual > hw_dual = { xi, rho, theta };
  hw_dual.xi.grad                                       = 1.0;
  const double dYddXi                                   = mwMC.yieldFunction( hw_dual ).grad;

  hw_dual.xi.grad      = 0.0;
  hw_dual.rho.grad     = 1.0;
  const double dYddRho = mwMC.yieldFunction( hw_dual ).grad;

  hw_dual.rho.grad       = 0.0;
  hw_dual.theta.grad     = 1.0;
  const double dYddTheta = mwMC.yieldFunction( hw_dual ).grad;

  throwExceptionOnFailure( checkIfEqual( std::get< 0 >( dYdHW ), dYddXi ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " failed: Derivative of yield function with respect to xi is wrong" );
  throwExceptionOnFailure( checkIfEqual( std::get< 1 >( dYdHW ), dYddRho ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " failed: Derivative of yield function with respect to rho is wrong" );
  throwExceptionOnFailure( checkIfEqual( std::get< 2 >( dYdHW ), dYddTheta ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " failed: Derivative of yield function with respect to theta is wrong" );
}

int main()
{
  testYieldFunctions();
  testSmoothingWithFillet();
  testInlineFunctions();
  testDerivatives();
  return 0;
}
