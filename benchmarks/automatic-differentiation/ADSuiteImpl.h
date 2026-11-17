#pragma once
#include "ADKernels.h"
#include "ADSuite.h"
#include <tuple>

namespace ADBench {

  namespace detail {
    template < std::size_t... Rest >
    void append( std::vector< double >& out, const Fastor::Tensor< double, Rest... >& t )
    {
      for ( Fastor::FASTOR_INDEX i = 0; i < t.size(); ++i )
        out.push_back( t.data()[t.get_mem_index( i )] );
    }
    template < typename D >
    void append( std::vector< double >& out, const Eigen::MatrixBase< D >& m )
    {
      for ( int i = 0; i < m.size(); ++i )
        out.push_back( m.derived().data()[i] );
    }
    inline void append( std::vector< double >& out, double d )
    {
      out.push_back( d );
    }
    inline void append( std::vector< double >& out, const std::vector< double >& c )
    {
      out.insert( out.end(), c.begin(), c.end() );
    }
    template < std::size_t N >
    void append( std::vector< double >& out, const std::array< double, N >& c )
    {
      out.insert( out.end(), c.begin(), c.end() );
    }
  } // namespace detail

  inline Tensor33d testDeformationGradient()
  {
    Tensor33d F = { { 1.05, 0.02, -0.01 }, { 0.03, 0.97, 0.015 }, { -0.02, 0.01, 1.02 } };
    return F;
  }

  inline Tensor33d testRightCauchyGreen()
  {
    const Tensor33d F = testDeformationGradient();
    Tensor33d       C;
    for ( int i = 0; i < 3; i++ )
      for ( int j = 0; j < 3; j++ ) {
        C( i, j ) = 0;
        for ( int k = 0; k < 3; k++ )
          C( i, j ) += F( k, i ) * F( k, j );
      }
    return C;
  }

  inline Vector6d testStressOld()
  {
    Vector6d s;
    s << 10., 5., -3., 2., 1., 0.5;
    return s;
  }

  inline Vector6d testStrainIncrement()
  {
    Vector6d e;
    e << 1e-3, -2e-4, 3e-4, 5e-4, -1e-4, 2e-4;
    return e;
  }

  template < typename... T >
  std::vector< double > flatten( const T&... parts )
  {
    std::vector< double > out;
    ( detail::append( out, parts ), ... );
    return out;
  }

  template < typename Tuple >
  std::vector< double > flattenTuple( const Tuple& t )
  {
    return std::apply( []( const auto&... parts ) { return flatten( parts... ); }, t );
  }

  inline const char* nameOf( const std::string& id )
  {
    if ( id == "tangent" )
      return "finite-strain tangent dTau/dF (9 dirs)";
    if ( id == "hypoJ2" )
      return "hypoelastic J2 tangent, Eigen (6 dirs)";
    if ( id == "hessian" )
      return "energy Hessian d2psi/dC2";
    if ( id == "third" )
      return "energy 3rd derivative d3psi/dC3";
    if ( id == "phasefield" )
      return "phase-field degradation (dual2nd scalar)";
    if ( id == "b4" )
      return "B4 drying Kelvin moduli (Taylor 5)";
    if ( id == "powerlaw" )
      return "power-law Kelvin moduli (Taylor 7)";
    return id.c_str();
  }

  template < typename B >
  std::vector< Suite > makeSuites()
  {
    std::vector< Suite > suites;
    auto                 add = []( Suite& s, const std::string& id, std::function< std::vector< double >() > f ) {
      s.cases.push_back( { id, nameOf( id ), std::move( f ) } );
    };

    // scalar mode: one direction per pass and nesting level, the way Marmot works today
    if constexpr ( B::maxScalarOrder >= 1 ) {
      Suite scalar{ B::name, {} };
      add( scalar, "tangent", [] {
        return flattenTuple( kernelFiniteStrainTangent< B >( testDeformationGradient() ) );
      } );
      add( scalar, "hypoJ2", [] {
        return flattenTuple( kernelHypoElasticJ2< B >( testStressOld(), testStrainIncrement() ) );
      } );
      if constexpr ( B::maxScalarOrder >= 2 ) {
        add( scalar, "hessian", [] { return flattenTuple( kernelEnergyHessian< B >( testRightCauchyGreen() ) ); } );
      }
      if constexpr ( B::maxScalarOrder >= 3 ) {
        add( scalar, "third", [] {
          return flattenTuple( kernelEnergyThirdDerivative< B >( testRightCauchyGreen() ) );
        } );
      }
      if constexpr ( B::maxScalarOrder >= 2 ) {
        add( scalar, "phasefield", [] { return flatten( kernelPhaseFieldDegradation< B >( 0.35 ) ); } );
      }
      if constexpr ( B::hasTaylor ) {
        add( scalar, "b4", [] { return flatten( kernelB4DryingKelvinModuli< B >() ); } );
        add( scalar, "powerlaw", [] { return flatten( kernelPowerLawKelvinModuli< B >() ); } );
      }
      suites.push_back( std::move( scalar ) );
    }

    // vector mode: all directions in one pass
    if constexpr ( B::hasVectorMode ) {
      Suite vec{ std::string( B::name ) + "-vec", {} };
      add( vec, "tangent", [] {
        return flattenTuple( kernelFiniteStrainTangentVector< B >( testDeformationGradient() ) );
      } );
      add( vec, "hypoJ2", [] {
        return flattenTuple( kernelHypoElasticJ2Vector< B >( testStressOld(), testStrainIncrement() ) );
      } );
      if constexpr ( B::maxVectorOrder >= 2 ) {
        add( vec, "hessian", [] { return flattenTuple( kernelEnergyHessianVector< B >( testRightCauchyGreen() ) ); } );
      }
      if constexpr ( B::maxVectorOrder >= 3 ) {
        add( vec, "third", [] {
          return flattenTuple( kernelEnergyThirdDerivativeVector< B >( testRightCauchyGreen() ) );
        } );
      }
      suites.push_back( std::move( vec ) );
    }

    // reverse mode (taped): first-order Jacobians from one recording
    if constexpr ( B::hasReverseMode ) {
      Suite rev{ std::string( B::name ) + "-rev", {} };
      if constexpr ( B::reverseFastorEinsumIsCorrect ) {
        add( rev, "tangent", [] { return flattenTuple( B::finiteStrainTangent( testDeformationGradient() ) ); } );
      }
      add( rev, "hypoJ2", [] { return flattenTuple( B::hypoElasticJ2( testStressOld(), testStrainIncrement() ) ); } );
      suites.push_back( std::move( rev ) );
    }
    return suites;
  }

} // namespace ADBench
