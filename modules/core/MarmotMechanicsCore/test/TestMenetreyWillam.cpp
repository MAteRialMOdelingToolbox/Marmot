#include "Marmot/HaighWestergaard.h"
#include "Marmot/MarmotConstants.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MenetreyWillam.h"
#include "autodiff/forward/dual.hpp"
#include <iomanip>

using namespace Marmot;
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

void testdPolarRadius_dTheta()
{
  const double e                        = 0.6;
  const double thetaExtensionMeridian   = 0.;
  const double thetaCompressionMeridian = Constants::Pi / 3;
  const double thetaShearMeridian       = Constants::Pi / 6;

  const auto [r_Ext, drdThetaExt] = MenetreyWillam::dPolarRadius_dTheta( thetaExtensionMeridian, e );
  const double r_ExtTarget        = 1. / e;
  const double drdThetaExtTarget  = 0.;

  const auto [r_Comp, drdThetaComp] = MenetreyWillam::dPolarRadius_dTheta( thetaCompressionMeridian, e );
  const double r_CompTarget         = 1.;
  const double drdThetaCompTarget   = 0.;

  const auto [r_Shear, drdThetaShear] = MenetreyWillam::dPolarRadius_dTheta( thetaShearMeridian, e );

  autodiff::dual thetaShearDual( thetaShearMeridian );
  thetaShearDual.grad                    = 1.0;
  const autodiff::dual rDualShear        = MenetreyWillam::dPolarRadius_dTheta( thetaShearDual, e ).first;
  const double         drdThetaDualShear = rDualShear.grad;

  throwExceptionOnFailure( checkIfEqual( r_Ext, r_ExtTarget ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " failed: Polar radius for extension meridian is wrong" );

  throwExceptionOnFailure( checkIfEqual( r_Comp, r_CompTarget ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " failed: Polar radius for compression meridian is wrong" );

  throwExceptionOnFailure( checkIfEqual( drdThetaExt, drdThetaExtTarget ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " failed: Derivative of polar radius for extension meridian is wrong" );

  throwExceptionOnFailure( checkIfEqual( drdThetaComp, drdThetaCompTarget, 1e-13 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " failed: Derivative of polar radius for compression meridian is wrong" );

  throwExceptionOnFailure( checkIfEqual( r_Shear, val( rDualShear ) ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " failed: Polar radius for shear meridian doesn't yield the same result for "
                                           "autodiff::dual" );

  throwExceptionOnFailure( checkIfEqual( drdThetaShear, drdThetaDualShear ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " failed: Derivative of polar radius for shear meridian is wrong" );
}

void testd2PolarRadius_dTheta2()
{
  const double e                        = 0.6;
  const double thetaExtensionMeridian   = 0.;
  const double thetaCompressionMeridian = Constants::Pi / 3;
  const double thetaShearMeridian       = Constants::Pi / 6;

  const auto [r_Ext, drdThetaExt, d2rdTheta2Ext] = MenetreyWillam::d2PolarRadius_dTheta2( thetaExtensionMeridian, e );
  const double r_ExtTarget                       = 1. / e;
  const double drdThetaExtTarget                 = 0.;

  const auto [r_Comp, drdThetaComp, d2rdTheta2Comp] = MenetreyWillam::d2PolarRadius_dTheta2( thetaCompressionMeridian,
                                                                                             e );
  const double r_CompTarget                         = 1.;
  const double drdThetaCompTarget                   = 0.;

  const auto [r_Shear, drdThetaShear, d2rdTheta2Shear] = MenetreyWillam::d2PolarRadius_dTheta2( thetaShearMeridian, e );

  autodiff::dual2nd thetaShearDual( thetaShearMeridian );
  auto              polarRadiusOnly = []( const autodiff::dual2nd& theta, const double& e ) {
    auto [r, dr, d2r] = MenetreyWillam::d2PolarRadius_dTheta2( theta, e );
    return r;
  };
  auto [rvalDual, drDual, d2rDual] = autodiff::derivatives( polarRadiusOnly,
                                                            autodiff::wrt( thetaShearDual ),
                                                            autodiff::at( thetaShearDual, e ) );

  throwExceptionOnFailure( checkIfEqual( r_Ext, r_ExtTarget ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " failed: Polar radius for extension meridian is wrong" );

  throwExceptionOnFailure( checkIfEqual( r_Comp, r_CompTarget ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " failed: Polar radius for compression meridian is wrong" );

  throwExceptionOnFailure( checkIfEqual( drdThetaExt, drdThetaExtTarget ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " failed: Derivative of polar radius for extension meridian is wrong" );

  throwExceptionOnFailure( checkIfEqual( drdThetaComp, drdThetaCompTarget, 1e-13 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " failed: Derivative of polar radius for compression meridian is wrong" );

  throwExceptionOnFailure( checkIfEqual( r_Shear, rvalDual ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " failed: Polar radius for shear meridian doesn't yield the same result for "
                                           "autodiff::dual" );

  throwExceptionOnFailure( checkIfEqual( drdThetaShear, drDual ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " failed: 1st derivative of polar radius for shear meridian is wrong" );

  throwExceptionOnFailure( checkIfEqual( d2rdTheta2Shear, d2rDual ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " failed: 2nd derivative of polar radius for shear meridian is wrong" );
}

int main()
{
  auto tests = std::vector< std::function< void() > >{ testYieldFunctions,
                                                       testSmoothingWithFillet,
                                                       testInlineFunctions,
                                                       testDerivatives,
                                                       testdPolarRadius_dTheta,
                                                       testd2PolarRadius_dTheta2 };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
