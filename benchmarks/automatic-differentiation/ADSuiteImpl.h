#pragma once
#include "ADKernels.h"
#include "ADSuite.h"

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
    inline void append( std::vector< double >& out, double d ) { out.push_back( d ); }
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

  template < typename B >
  Suite makeSuite()
  {
    using detail::append;
    Suite s;
    s.backend = B::name;

    s.cases.push_back( { "finite-strain tangent (dual, 9 passes)",
                         "ADCompressibleNeoHooke via MaterialFiniteStrainAD::dF_dT",
                         [] {
                           std::vector< double > out;
                           auto [tau, dTau_dF] = kernelFiniteStrainTangent< B >( testDeformationGradient() );
                           append( out, tau );
                           append( out, dTau_dF );
                           return out;
                         } } );

    s.cases.push_back( { "hypoelastic J2 tangent (Eigen dual, 6 passes)",
                         "MaterialHypoElasticAD::dF_dX (ADVonMises-like)",
                         [] {
                           std::vector< double > out;
                           Vector6d              sOld, dE;
                           sOld << 10., 5., -3., 2., 1., 0.5;
                           dE << 1e-3, -2e-4, 3e-4, 5e-4, -1e-4, 2e-4;
                           auto [S, C] = kernelHypoElasticJ2< B >( sOld, dE );
                           append( out, S );
                           append( out, C );
                           return out;
                         } } );

    s.cases.push_back( { "energy Hessian (dual2nd, 81 passes)",
                         "SecondOrder::d2f_dT2 on PenceGouPotentialB",
                         [] {
                           std::vector< double > out;
                           auto [psi, dPsi, d2Psi] = kernelEnergyHessian< B >( testRightCauchyGreen() );
                           append( out, psi );
                           append( out, dPsi );
                           append( out, d2Psi );
                           return out;
                         } } );

    s.cases.push_back( { "energy 3rd derivative (dual3rd, 729 passes)",
                         "ThirdOrder::d3f_dT3 (CompressibleFiniteStrainLinearViscoelasticity)",
                         [] {
                           std::vector< double > out;
                           auto [psi, d1, d2, d3] = kernelEnergyThirdDerivative< B >( testRightCauchyGreen() );
                           append( out, psi );
                           append( out, d1 );
                           append( out, d2 );
                           append( out, d3 );
                           return out;
                         } } );

    s.cases.push_back( { "phase-field degradation (dual2nd scalar)",
                         "EnergyDegradationFunctions::SecondOrderDerived::generic",
                         [] {
                           std::vector< double > out;
                           append( out, kernelPhaseFieldDegradation< B >( 0.35 ) );
                           return out;
                         } } );

    s.cases.push_back( { "B4 drying Kelvin moduli (Taylor 5, 8 taus)",
                         "B4::computeStress -> KelvinChain::computeElasticModuli<5>",
                         [] {
                           std::vector< double > out;
                           append( out, kernelB4DryingKelvinModuli< B >() );
                           return out;
                         } } );

    s.cases.push_back( { "power-law Kelvin moduli (Taylor 7, 12 taus)",
                         "LinearViscoelasticOrthotropicPowerLaw -> computeElasticModuli<7>",
                         [] {
                           std::vector< double > out;
                           append( out, kernelPowerLawKelvinModuli< B >() );
                           return out;
                         } } );

    return s;
  }

} // namespace ADBench
