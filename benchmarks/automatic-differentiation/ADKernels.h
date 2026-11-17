#pragma once
// Backend-generic AD drivers (ports of Marmot/MarmotAutomaticDifferentiation*.h) and the benchmark kernels.
// Every kernel mirrors a real Marmot code path and, where possible, calls the actual Marmot templates.
#include "Marmot/MarmotDeformationMeasures.h"
#include "Marmot/MarmotEnergyDensityFunctions.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotStressMeasures.h"
#include "Marmot/MarmotViscoelasticity.h"
#include <Eigen/Core>
#include <array>
#include <cmath>
#include <functional>
#include <tuple>
#include <vector>

namespace ADBench {

  using Fastor::Tensor;
  using Tensor33d     = Tensor< double, 3, 3 >;
  using Tensor3333d   = Tensor< double, 3, 3, 3, 3 >;
  using Tensor333333d = Tensor< double, 3, 3, 3, 3, 3, 3 >;
  using Vector6d      = Eigen::Matrix< double, 6, 1 >;
  using Matrix6d      = Eigen::Matrix< double, 6, 6 >;

  // -------------------------------------------------------------------------------------------------------------
  // Drivers
  // -------------------------------------------------------------------------------------------------------------

  template < typename S, std::size_t... Rest >
  Tensor< S, Rest... > lift( const Tensor< double, Rest... >& in )
  {
    Tensor< S, Rest... > out;
    for ( Fastor::FASTOR_INDEX i = 0; i < in.size(); ++i )
      out.data()[out.get_mem_index( i )] = S( in.data()[in.get_mem_index( i )] );
    return out;
  }

  /// Marmot::AutomaticDifferentiation::dF_dT (tensor -> tensor, first order; MarmotMaterialFiniteStrainAD hot path)
  template < typename B >
  std::pair< Tensor33d, Tensor3333d > dF_dT(
    const std::function< Tensor< typename B::template HyperDual< 1 >, 3, 3 >(
      const Tensor< typename B::template HyperDual< 1 >, 3, 3 >& ) >& F,
    const Tensor33d&                                                  T )
  {
    using S = typename B::template HyperDual< 1 >;
    Tensor33d            F_;
    Tensor3333d          dF_dT_;
    Tensor< S, 3, 3 >    T_ = lift< S >( T );
    Tensor< S, 3, 3 >    FT;
    for ( std::size_t i = 0; i < 9; ++i ) {
      B::template seed< 1 >( T_.data()[i], 1.0 );
      FT = F( T_ );
      for ( std::size_t j = 0; j < 9; ++j )
        dF_dT_.data()[j * 9 + i] = B::template derivative< 1 >( FT.data()[j] );
      B::template seed< 1 >( T_.data()[i], 0.0 );
    }
    for ( std::size_t j = 0; j < 9; ++j )
      F_.data()[j] = B::value( FT.data()[j] );
    return { F_, dF_dT_ };
  }

  /// Marmot::AutomaticDifferentiation::SecondOrder::d2f_dT2 (tensor -> scalar, hyper-dual 2nd order)
  template < typename B >
  std::tuple< double, Tensor33d, Tensor3333d > d2f_dT2(
    const std::function< typename B::template HyperDual< 2 >( const Tensor< typename B::template HyperDual< 2 >, 3, 3 >& ) >&
                     F,
    const Tensor33d& T )
  {
    using S = typename B::template HyperDual< 2 >;
    double            F_;
    S                 Fr;
    Tensor33d         dF;
    Tensor3333d       d2F;
    Tensor< S, 3, 3 > T_ = lift< S >( T );
    for ( std::size_t i = 0; i < 3; i++ )
      for ( std::size_t j = 0; j < 3; j++ ) {
        B::template seed< 1 >( T_( i, j ), 1.0 );
        for ( std::size_t k = 0; k < 3; k++ )
          for ( std::size_t l = 0; l < 3; l++ ) {
            B::template seed< 2 >( T_( k, l ), 1.0 );
            Fr                = F( T_ );
            d2F( i, j, k, l ) = B::template derivative< 2 >( Fr );
            B::template seed< 2 >( T_( k, l ), 0.0 );
          }
        dF( i, j ) = B::template derivative< 1 >( Fr );
        F_         = B::value( Fr );
        B::template seed< 1 >( T_( i, j ), 0.0 );
      }
    return { F_, dF, d2F };
  }

  /// Marmot::AutomaticDifferentiation::ThirdOrder::d3f_dT3 (tensor -> scalar, hyper-dual 3rd order, 729 passes)
  template < typename B >
  std::tuple< double, Tensor33d, Tensor3333d, Tensor333333d > d3f_dT3(
    const std::function< typename B::template HyperDual< 3 >( const Tensor< typename B::template HyperDual< 3 >, 3, 3 >& ) >&
                     F,
    const Tensor33d& T )
  {
    using S = typename B::template HyperDual< 3 >;
    double            F_;
    S                 Fr;
    Tensor33d         dF;
    Tensor3333d       d2F;
    Tensor333333d     d3F;
    Tensor< S, 3, 3 > T_ = lift< S >( T );
    for ( std::size_t i = 0; i < 3; i++ )
      for ( std::size_t j = 0; j < 3; j++ ) {
        B::template seed< 1 >( T_( i, j ), 1.0 );
        for ( std::size_t k = 0; k < 3; k++ )
          for ( std::size_t l = 0; l < 3; l++ ) {
            B::template seed< 2 >( T_( k, l ), 1.0 );
            for ( std::size_t m = 0; m < 3; m++ )
              for ( std::size_t n = 0; n < 3; n++ ) {
                B::template seed< 3 >( T_( m, n ), 1.0 );
                Fr                      = F( T_ );
                d3F( i, j, k, l, m, n ) = B::template derivative< 3 >( Fr );
                B::template seed< 3 >( T_( m, n ), 0.0 );
              }
            d2F( i, j, k, l ) = B::template derivative< 2 >( Fr );
            B::template seed< 2 >( T_( k, l ), 0.0 );
          }
        dF( i, j ) = B::template derivative< 1 >( Fr );
        F_         = B::value( Fr );
        B::template seed< 1 >( T_( i, j ), 0.0 );
      }
    return { F_, dF, d2F, d3F };
  }

  /// Marmot::AutomaticDifferentiation::dF_dX (Eigen vector -> vector; MarmotMaterialHypoElasticAD hot path)
  template < typename B >
  std::pair< Vector6d, Matrix6d > dF_dX(
    const std::function< Eigen::Matrix< typename B::template HyperDual< 1 >, 6, 1 >(
      const Eigen::Matrix< typename B::template HyperDual< 1 >, 6, 1 >& ) >& F,
    const Vector6d&                                                         X )
  {
    using S = typename B::template HyperDual< 1 >;
    Eigen::Matrix< S, 6, 1 > X_ = X.template cast< S >();
    Eigen::Matrix< S, 6, 1 > Fr;
    Matrix6d                 J;
    Vector6d                 F_;
    for ( int j = 0; j < 6; j++ ) {
      B::template seed< 1 >( X_( j ), 1.0 );
      Fr = F( X_ );
      for ( int i = 0; i < 6; i++ )
        J( i, j ) = B::template derivative< 1 >( Fr( i ) );
      B::template seed< 1 >( X_( j ), 0.0 );
      F_( j ) = B::value( Fr( j ) );
    }
    return { F_, J };
  }

  // -------------------------------------------------------------------------------------------------------------
  // Kernels
  // -------------------------------------------------------------------------------------------------------------

  inline constexpr double K = 3500., G = 1500.;

  /// ADCompressibleNeoHooke::computeStressAD + MarmotMaterialFiniteStrainAD: Kirchhoff stress and dTau/dF
  template < typename B >
  std::pair< Tensor33d, Tensor3333d > kernelFiniteStrainTangent( const Tensor33d& F )
  {
    using S = typename B::template HyperDual< 1 >;
    using namespace Marmot::ContinuumMechanics;
    std::function< Tensor< S, 3, 3 >( const Tensor< S, 3, 3 >& ) > tau = [&]( const Tensor< S, 3, 3 >& F_ ) {
      const auto C                = DeformationMeasures::rightCauchyGreen( F_ );
      const auto [psi_, dPsi_dC]  = EnergyDensityFunctions::FirstOrderDerived::PenceGouPotentialB( C, K, G );
      Tensor< S, 3, 3 >      PK2 = Marmot::multiplyFastorTensorWithScalar( dPsi_dC, S( 2.0 ) );
      return StressMeasures::KirchhoffStressFromPK2( PK2, F_ );
    };
    return dF_dT< B >( tau, F );
  }

  /// SecondOrder::d2f_dT2 on the Pence-Gou energy (energy, PK2/2, material tangent)
  template < typename B >
  std::tuple< double, Tensor33d, Tensor3333d > kernelEnergyHessian( const Tensor33d& C )
  {
    using S = typename B::template HyperDual< 2 >;
    std::function< S( const Tensor< S, 3, 3 >& ) > psi = []( const Tensor< S, 3, 3 >& C_ ) {
      return Marmot::ContinuumMechanics::EnergyDensityFunctions::PenceGouPotentialB< S >( C_, K, G );
    };
    return d2f_dT2< B >( psi, C );
  }

  /// CompressibleFiniteStrainLinearViscoelasticity: ThirdOrder::d3f_dT3 on the Yeoh energy (729 passes)
  template < typename B >
  std::tuple< double, Tensor33d, Tensor3333d, Tensor333333d > kernelEnergyThirdDerivative( const Tensor33d& C )
  {
    using S = typename B::template HyperDual< 3 >;
    std::function< S( const Tensor< S, 3, 3 >& ) > psi = []( const Tensor< S, 3, 3 >& C_ ) {
      return Marmot::ContinuumMechanics::EnergyDensityFunctions::YeohPotential< S >( C_, 500., -10., 1., 1e-3 );
    };
    return d3f_dT3< B >( psi, C );
  }

  /// MarmotMaterialHypoElasticAD + ADVonMises-like radial return (Eigen Vector6 of duals, 6 passes)
  template < typename B >
  std::pair< Vector6d, Matrix6d > kernelHypoElasticJ2( const Vector6d& stressOld, const Vector6d& dStrain )
  {
    using S  = typename B::template HyperDual< 1 >;
    using V6 = Eigen::Matrix< S, 6, 1 >;
    Matrix6d Cel                                  = Matrix6d::Zero();
    const double E = 30000., nu = 0.2, fy = 30.;
    const double lam = E * nu / ( 1 + nu ) / ( 1 - 2 * nu ), mu = E / 2 / ( 1 + nu );
    Cel.topLeftCorner< 3, 3 >().setConstant( lam );
    Cel.diagonal().head< 3 >().array() += 2 * mu;
    Cel.diagonal().tail< 3 >().setConstant( mu );

    std::function< V6( const V6& ) > stress = [&]( const V6& dE ) {
      V6      S_  = stressOld.template cast< S >() + Cel * dE;
      const S p   = ( S_( 0 ) + S_( 1 ) + S_( 2 ) ) / 3.;
      V6      dev = S_;
      for ( int i = 0; i < 3; i++ )
        dev( i ) -= p;
      S J2 = 0.5 * ( dev( 0 ) * dev( 0 ) + dev( 1 ) * dev( 1 ) + dev( 2 ) * dev( 2 ) ) +
             ( dev( 3 ) * dev( 3 ) + dev( 4 ) * dev( 4 ) + dev( 5 ) * dev( 5 ) );
      const S q = sqrt( 3. * J2 );
      if ( B::value( q ) > fy ) {
        const S scale = fy / q;
        for ( int i = 0; i < 6; i++ )
          dev( i ) *= scale;
        for ( int i = 0; i < 3; i++ )
          dev( i ) += p;
        return V6( dev );
      }
      return S_;
    };
    return dF_dX< B >( stress, dStrain );
  }

  /// MarmotMaterialPhaseFieldEnergyDegradation::generic (scalar, dual2nd with both levels seeded)
  template < typename B >
  std::array< double, 3 > kernelPhaseFieldDegradation( double pf )
  {
    using S = typename B::template HyperDual< 2 >;
    const double p = 2., a1 = 2.5, a2 = -0.5, a3 = 0.0;
    S            x( pf );
    B::template seed< 1 >( x, 1.0 );
    B::template seed< 2 >( x, 1.0 );
    const S s      = 1. - x;
    const S P      = 1. + a2 * x * ( 1. + a3 * x );
    const S Q      = a1 * x * P;
    const S result = pow( s, p ) / ( pow( s, p ) + Q );
    return { B::value( result ), B::template derivative< 1 >( result ), B::template derivative< 2 >( result ) };
  }

  /// Post-Widder L_k(tau) for a list of retardation times (KelvinChain::computeElasticModuli)
  template < typename B, int k, typename Phi >
  std::vector< double > postWidderModuli( const Phi& phi, const std::vector< double >& taus )
  {
    std::vector< double > E( taus.size() );
    double                fact = 1.;
    for ( int i = 2; i < k; i++ )
      fact *= i;
    const double spacing = taus[1] / taus[0];
    for ( std::size_t i = 0; i < taus.size(); i++ ) {
      const double tau = taus[i];
      double       val = -std::pow( -tau * k, k ) / fact;
      val *= B::template taylorDerivatives< k >( phi, tau * k )[k];
      E[i] = 1. / ( std::log( spacing ) * val );
    }
    return E;
  }

  /// B4::computeStress drying creep: Post-Widder of order 5 on the drying compliance, evaluated per stress update
  template < typename B >
  std::vector< double > kernelB4DryingKelvinModuli()
  {
    constexpr int         order = 5;
    const double          b = 8. * ( 1. - 0.6 ), xiZero = -0.3;
    std::vector< double > taus;
    for ( int i = 0; i < 8; i++ )
      taus.push_back( 1e-3 * std::pow( std::sqrt( 10. ), i ) );
    auto phi = [&]( typename B::template Taylor< order > xi ) {
      return sqrt( exp( tanh( sqrt( xi - xiZero ) ) * b ) - exp( tanh( sqrt( -xiZero ) ) * b ) );
    };
    return postWidderModuli< B, order >( phi, taus );
  }

  /// LinearViscoelasticOrthotropicPowerLaw: power-law compliance, Post-Widder order 7
  template < typename B >
  std::vector< double > kernelPowerLawKelvinModuli()
  {
    constexpr int         order = 7;
    std::vector< double > taus;
    for ( int i = 0; i < 12; i++ )
      taus.push_back( 1e-4 * std::pow( 10., i ) );
    auto phi = []( typename B::template Taylor< order > tau ) {
      return Marmot::ContinuumMechanics::Viscoelasticity::ComplianceFunctions::powerLaw( tau, 1e-4, 0.25 );
    };
    return postWidderModuli< B, order >( phi, taus );
  }

} // namespace ADBench
