#include "Marmot/MarmotConstants.h"
#include "Marmot/MenetreyWillam.h"
#include <cmath>
#include <sstream>

namespace Marmot {
  namespace ContinuumMechanics::CommonConstitutiveModels {
    using namespace Constants;
    using namespace ContinuumMechanics::HaighWestergaard;

    MenetreyWillam::MenetreyWillam( const double ft, const MenetreyWillamType& type, const double fc )
    {
      setParameters( ft, fc, type );
    }

    void MenetreyWillam::setParameters( const double ft, const double fc, const MenetreyWillamType& type )
    {
      switch ( type ) {
      case MenetreyWillamType::Mises:
        param.Af = 0;
        param.Bf = sqrt3_2 / ft;
        param.Cf = 0;
        param.m  = 1;
        param.e  = 1;
        break;
      case MenetreyWillamType::Rankine:
        param.Af = 0;
        param.Bf = 1. / ( sqrt6 * ft );
        param.Cf = 1. / ( sqrt3 * ft );
        param.m  = 1;
        param.e  = 0.51;
        break;
      case MenetreyWillamType::DruckerPrager:
        param.Af = 0;
        param.Bf = sqrt3_8 * ( fc + ft ) / ( fc * ft );
        param.Cf = 3. / 2 * ( fc - ft ) / ( fc * ft );
        param.m  = 1;
        param.e  = 1;
        break;
      case MenetreyWillamType::MohrCoulomb:
        param.Af = 0;
        param.Bf = 1. / sqrt6 * ( fc + 2. * ft ) / ( fc * ft );
        param.Cf = 1. / sqrt3 * ( fc - ft ) / ( fc * ft );
        param.m  = 1.;
        param.e  = ( fc + 2 * ft ) / ( 2 * fc + ft );
        break;
      default: throw std::invalid_argument( "Requested MenetreyWillamType not found." );
      }
    }

    double MenetreyWillam::polarRadius( const double& theta ) const { return polarRadius( theta, param.e ); }

    double MenetreyWillam::polarRadius( const double& theta, const double& e )
    {
      // computes the deviatoric shape roundness for a given eccentricity (e) at a certain
      // position (theta) the numerator and denominator are stored for performance reasons, as
      // they are also needed for the derivative dRdTheta
      if ( e >= 1.0 )
        return 1;

      const double e2        = e * e;
      const double cosTheta  = std::cos( theta );
      const double cos2Theta = cosTheta * cosTheta;

      const double numerator   = ( 4 * ( 1 - e2 ) * cos2Theta + ( 2 * e - 1 ) * ( 2 * e - 1 ) );
      const double denominator = ( 2 * ( 1 - e2 ) * cosTheta +
                                   ( 2 * e - 1 ) * std::sqrt( 4 * ( 1 - e2 ) * cos2Theta + 5 * e2 - 4 * e ) );

      return numerator / denominator;
    }

    std::pair< double, double > MenetreyWillam::dPolarRadius_dTheta( const double& theta ) const
    {
      return dPolarRadius_dTheta( theta, param.e );
    }

    std::pair< double, double > MenetreyWillam::dPolarRadius_dTheta( const double& theta, const double& e )
    {
      //
      // computes the deviatoric shape roundness for a given eccentricity (e) at a certain
      // position (theta) the numerator and denominator are stored for performance reasons, as
      // they are also needed for the derivative dRdTheta
      double r, dRdTheta;

      if ( e >= 1.0 ) {
        r        = 1;
        dRdTheta = 0;
        return std::make_pair( r, dRdTheta );
      }

      const double e2        = e * e;
      const double cosTheta  = std::cos( theta );
      const double cos2Theta = cosTheta * cosTheta;

      double numerator   = ( 4 * ( 1 - e2 ) * cos2Theta + ( 2 * e - 1 ) * ( 2 * e - 1 ) );
      double denominator = ( 2 * ( 1 - e2 ) * cosTheta +
                             ( 2 * e - 1 ) * std::sqrt( 4 * ( 1 - e2 ) * cos2Theta + 5 * e2 - 4 * e ) );

      r = numerator / denominator;

      // compute the derivative
      const double sinTheta = std::sin( theta );

      const double a          = numerator;
      const double b          = 1. / denominator;
      const double aux        = 4 * ( 1 - e2 ) * cos2Theta + 5 * e2 - 4 * e;
      const double dAuxdTheta = 4 * ( 1 - e2 ) * 2 * cosTheta * -sinTheta;
      const double dadTheta   = 2 * cosTheta * ( -sinTheta ) * 4 * ( 1 - e2 );
      const double dbdTheta   = -pow( denominator, -2 ) * ( 2 * ( 1 - e2 ) * -sinTheta +
                                                          ( 2 * e - 1 ) * 1. / 2 * pow( aux, -1. / 2 ) * dAuxdTheta );
      dRdTheta                = b * dadTheta + a * dbdTheta;

      return { r, dRdTheta };
    }

    double MenetreyWillam::yieldFunction( const HaighWestergaardCoordinates<>& hw, const double varEps ) const
    {
      const double r_ = polarRadius( hw.theta, param.e );
      if ( varEps == 0 )
        return ( param.Af * hw.rho ) * ( param.Af * hw.rho ) + param.m * ( param.Bf * hw.rho * r_ + param.Cf * hw.xi ) -
               1;
      else
        return param.Af * param.Af * hw.rho * hw.rho +
               param.m *
                 ( std::sqrt( param.Bf * hw.rho * r_ * param.Bf * hw.rho * r_ + varEps * varEps ) + param.Cf * hw.xi ) -
               1;
    }

    std::tuple< double, double, double > MenetreyWillam::dYieldFunction_dHaighWestergaard(
      const HaighWestergaardCoordinates<>& hw,
      const double                         varEps ) const
    {
      const auto [r_, dRdTheta_] = dPolarRadius_dTheta( hw.theta, param.e );

      double dFdXi, dFdRho, dFdTheta;
      dFdXi = param.m * param.Cf;

      if ( varEps == 0.0 ) {
        dFdRho   = 2 * param.Af * hw.rho * param.Af + param.m * param.Bf * r_;
        dFdTheta = param.m * param.Bf * hw.rho * dRdTheta_;
      }
      else {
        const double auxTerm1 = param.m * 0.5 *
                                std::pow( param.Bf * param.Bf * hw.rho * hw.rho * r_ * r_ + varEps * varEps, -0.5 ) *
                                2 * param.Bf * hw.rho * r_ * param.Bf;

        dFdRho   = param.Af * param.Af * 2 * hw.rho + auxTerm1 * r_;
        dFdTheta = auxTerm1 * hw.rho * dRdTheta_;
      }

      return { dFdXi, dFdRho, dFdTheta };
    }

  } // namespace ContinuumMechanics::CommonConstitutiveModels
} // namespace Marmot
