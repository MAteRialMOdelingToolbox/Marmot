#include "Marmot/MarmotBondSlipLawFactory.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/ModelCode2010BondSlip.h"
#include <cmath>
#include <memory>
#include <string>
#include <utility>
#include <vector>

using namespace Marmot;
using namespace Marmot::Testing;
using namespace Marmot::Math;
using Marmot::Materials::ModelCode2010BondSlip;

namespace {
  // MC2010 Table 6.1-1, good bond conditions, unconfined, pull-out failure, fck = 30 MPa (units: N, mm)
  const double tauMax  = 2.5 * std::sqrt( 30. );
  const double Kn      = 1000.;
  const double props[] = { tauMax, 1.0, 2.0, 10.0, 0.4, 0.4 * tauMax, 200., Kn };
  // with Coulomb friction mu = 0.4 and crushing at 3 fck
  const double mu = 0.4, pMax = 90.;
  const double frictionProps[] = { tauMax, 1.0, 2.0, 10.0, 0.4, 0.4 * tauMax, 200., Kn, mu, pMax };

  std::unique_ptr< ModelCode2010BondSlip > create( const double* p = props, int n = 8 )
  {
    return std::unique_ptr< ModelCode2010BondSlip >( dynamic_cast< ModelCode2010BondSlip* >(
      Marmot::Registry::BondSlipLawFactory::createBondSlipLaw( "MODELCODE2010BONDSLIP", p, n, 1 ) ) );
  }

  std::unique_ptr< ModelCode2010BondSlip > createWithFriction()
  {
    return create( frictionProps, 10 );
  }

  /// a bond point with its own state
  struct BondPoint {
    const Marmot::Materials::BondSlipLaw& law;
    std::vector< double >                 stateVars = std::vector< double >( law.getNumberOfRequiredStateVars(), 0.0 );
    Vector3d                              stress    = Vector3d::Zero();
    double                                energy = 0.0, dissipation = 0.0;
    Vector3d                              slip = Vector3d::Zero();
    Matrix3d                              C    = Matrix3d::Zero();

    void update( const Vector3d& s )
    {
      Marmot::Materials::BondSlipLaw::State state{ stress, energy, dissipation, stateVars.data() };
      law.computeBondStress( state, C, s, s - slip, { 0.0, 1.0 } );
      stress      = state.bondStress;
      energy      = state.elasticEnergyDensity;
      dissipation = state.dissipation;
      slip        = s;
    }

    Vector3d trial( const Vector3d& s ) const
    {
      auto                                  sv = stateVars;
      Marmot::Materials::BondSlipLaw::State state{ stress, energy, dissipation, sv.data() };
      Matrix3d                              C_;
      law.computeBondStress( state, C_, s, s - slip, { 0.0, 1.0 } );
      return state.bondStress;
    }

    void checkTangent( const Vector3d& s, const std::string& where ) const
    {
      auto                                  sv = stateVars;
      Marmot::Materials::BondSlipLaw::State state{ stress, energy, dissipation, sv.data() };
      Matrix3d                              C_;
      law.computeBondStress( state, C_, s, s - slip, { 0.0, 1.0 } );
      Matrix3d     numC;
      const double h = 1e-8;
      for ( int j = 0; j < 3; j++ ) {
        Vector3d dh   = Vector3d::Zero();
        dh( j )       = h;
        numC.col( j ) = ( trial( s + dh ) - trial( s - dh ) ) / ( 2 * h );
      }
      throwExceptionOnFailure( checkIfEqual< double >( C_, numC, 1e-5 ),
                               Marmot::Utilities::MakeString() << "tangent at " << where << ":\n"
                                                               << C_ << "\nvs numerical\n"
                                                               << numC );
    }
  };
} // namespace

/// monotonic loading in small steps reproduces the envelope, including the ascending power law, the plateau, the
/// softening and the residual branch
void testMonotonicEnvelope()
{
  auto      law = create();
  BondPoint p{ *law };
  for ( int i = 1; i <= 1500; i++ ) {
    const double s = i * 0.01;
    p.update( Vector3d( s, 0, 0 ) );

    double expected;
    if ( s <= law->s0 )
      expected = law->K0 * s;
    else if ( s <= 1.0 )
      expected = tauMax * std::pow( s, 0.4 );
    else if ( s <= 2.0 )
      expected = tauMax;
    else if ( s <= 10.0 )
      expected = tauMax - 0.6 * tauMax * ( s - 2.0 ) / 8.0;
    else
      expected = 0.4 * tauMax;

    throwExceptionOnFailure( checkIfEqual( p.stress( 0 ), expected, 1e-10 ),
                             Marmot::Utilities::MakeString()
                               << "envelope at s=" << s << ": " << p.stress( 0 ) << " != " << expected );
  }
  throwExceptionOnFailure( p.dissipation > 0, "monotonic slip must dissipate" );
}

/// unloading is elastic with K0, leaves a residual slip, and the reversed direction is bounded by the already reached
/// (damaged) bond stress
void testUnloadingAndReversal()
{
  auto      law = create();
  BondPoint p{ *law };
  for ( int i = 1; i <= 50; i++ )
    p.update( Vector3d( i * 0.1, 0, 0 ) ); // s = 5: softening branch
  const double tauAt5 = p.stress( 0 );
  const double sP     = *law->getStateView( "plastic slip", p.stateVars.data() ).stateLocation;
  throwExceptionOnFailure( checkIfEqual( sP, 5.0 - tauAt5 / law->K0, 1e-12 ), "plastic slip" );

  p.update( Vector3d( 4.95, 0, 0 ) );
  throwExceptionOnFailure( checkIfEqual( p.stress( 0 ), tauAt5 - law->K0 * 0.05, 1e-10 ) &&
                             checkIfEqual( p.C( 0, 0 ), law->K0, 1e-12 ),
                           "elastic unloading with K0" );

  // full reversal: frictional, bounded by the envelope at the largest slip
  for ( int i = 1; i <= 100; i++ )
    p.update( Vector3d( 4.95 - i * 0.05, 0, 0 ) );
  throwExceptionOnFailure( checkIfEqual( p.stress( 0 ), -tauAt5, 1e-10 ) && checkIfEqual( p.C( 0, 0 ), 0.0, 1e-14 ),
                           Marmot::Utilities::MakeString()
                             << "reversed friction " << p.stress( 0 ) << " != " << -tauAt5 );
}

/// consistent tangents on every branch, from committed states
void testTangents()
{
  auto      law = create();
  BondPoint p{ *law };
  p.checkTangent( Vector3d( 0.001, 0.01, -0.02 ), "virgin linear branch" );
  p.checkTangent( Vector3d( 0.3, 0.01, -0.02 ), "virgin ascending branch" );
  p.update( Vector3d( 0.3, 0, 0 ) );
  p.checkTangent( Vector3d( 0.5, 0, 0 ), "ascending branch" );
  p.checkTangent( Vector3d( 0.29, 0, 0 ), "unloading" );
  p.update( Vector3d( 3.0, 0, 0 ) );
  p.checkTangent( Vector3d( 3.5, 0, 0 ), "softening" );
  p.checkTangent( Vector3d( -1.0, 0, 0 ), "reversed friction" );
  p.checkTangent( Vector3d( -4.0, 0, 0 ), "reversed virgin softening" );
  p.update( Vector3d( 12.0, 0, 0 ) );
  p.checkTangent( Vector3d( 13.0, 0, 0 ), "residual" );
}

void testInvalidPropertiesThrow()
{
  const std::vector< double > good = { tauMax, 1.0, 2.0, 10.0, 0.4, 0.4 * tauMax, 200., 1000. };
  // each case breaks one requirement
  std::vector< std::pair< std::string, std::vector< double > > > cases = {
    { "too few properties", { tauMax, 1.0, 2.0 } },
  };
  auto with = [&]( int index, double value ) {
    auto p   = good;
    p[index] = value;
    return p;
  };
  cases.push_back( { "s1 <= 0", with( 1, 0.0 ) } );
  cases.push_back( { "s2 >= s3", with( 2, 20.0 ) } );
  cases.push_back( { "alpha > 1", with( 4, 1.5 ) } );
  cases.push_back( { "tauF > tauMax", with( 5, 2 * tauMax ) } );
  cases.push_back( { "K0 <= tauMax / s1", with( 6, 0.5 * tauMax ) } );
  cases.push_back( { "Kn < 0", with( 7, -1.0 ) } );
  for ( const auto& [what, badProps] : cases ) {
    bool threw = false;
    try {
      Marmot::Registry::BondSlipLawFactory::createBondSlipLaw( "MODELCODE2010BONDSLIP",
                                                               badProps.data(),
                                                               static_cast< int >( badProps.size() ),
                                                               1 );
    }
    catch ( const std::invalid_argument& ) {
      threw = true;
    }
    throwExceptionOnFailure( threw, what + " must be rejected" );
  }
}

/// explicit mu = 0 and no cap (10 properties) equals the default of 8 properties on a loading/unloading/reversal path
/// with normal slip
void testOptionalFrictionPropertiesDefaultToZero()
{
  const double p0[] = { tauMax, 1.0, 2.0, 10.0, 0.4, 0.4 * tauMax, 200., Kn, 0.0, 0.0 };
  auto         a    = create( p0, 10 );
  auto         b    = create();
  BondPoint    pa{ *a }, pb{ *b };
  for ( int i = 1; i <= 400; i++ ) {
    const double   s = i <= 200 ? i * 0.03 : 6.0 - ( i - 200 ) * 0.05;
    const Vector3d slip( s, 0.01 * std::sin( 0.1 * i ), -0.02 * std::cos( 0.07 * i ) );
    pa.update( slip );
    pb.update( slip );
    throwExceptionOnFailure( checkIfEqual< double >( pa.stress, pb.stress, 1e-12 ) &&
                               checkIfEqual< double >( pa.C, pb.C, 1e-12 ),
                             Marmot::Utilities::MakeString()
                               << "explicit mu = 0 differs from the default at step " << i );
  }
}

/// constant normal pressure raises the whole envelope by mu * p; the cap limits p (and the friction) to pMax
void testCoulombEnvelopeAndCap()
{
  auto law = createWithFriction();
  for ( const double sn : { 0.0, 0.02, -0.05, 0.5 } ) {
    BondPoint    p{ *law };
    const double pressure = std::min( Kn * std::abs( sn ), pMax );
    for ( int i = 1; i <= 1500; i++ ) {
      const double s = i * 0.01;
      p.update( Vector3d( s, sn, 0.0 ) );
      const double env = law->envelope( s ).first;
      // the elastic K0 branch: below the bound, tau = K0 s
      const double expected = std::min( 200. * s, env + mu * pressure );
      throwExceptionOnFailure( checkIfEqual( p.stress( 0 ), expected, 1e-9 ) &&
                                 checkIfEqual( std::abs( p.stress( 1 ) ), pressure, 1e-12 ),
                               Marmot::Utilities::MakeString() << "envelope at s=" << s << ", sn=" << sn << ": "
                                                               << p.stress( 0 ) << " != " << expected );
    }
  }
  // crushed concrete around the bar: unloading the normal slip from 0.5 leaves a plastic normal slip 0.5 - pMax / Kn
  BondPoint p{ *law };
  p.update( Vector3d( 0.0, 0.5, 0.0 ) );
  p.update( Vector3d( 0.0, 0.45, 0.0 ) );
  throwExceptionOnFailure( checkIfEqual( p.stress( 1 ), pMax - Kn * 0.05, 1e-10 ), "elastic normal unloading" );
  p.update( Vector3d( 0.0, 0.5 - pMax / Kn, 0.0 ) );
  throwExceptionOnFailure( checkIfEqual( p.stress( 1 ), 0.0, 1e-10 ), "plastic normal slip" );
  throwExceptionOnFailure( p.dissipation > 0, "crushing dissipates" );
}

/// consistent tangents on every branch, from committed states
void testTangentsWithFriction()
{
  auto      law = createWithFriction();
  BondPoint p{ *law };
  p.checkTangent( Vector3d( 0.001, 0.01, -0.02 ), "virgin linear branch" );
  p.checkTangent( Vector3d( 0.3, 0.01, -0.02 ), "virgin ascending branch with pressure" );
  p.update( Vector3d( 0.3, 0.01, -0.02 ) );
  p.checkTangent( Vector3d( 0.5, 0.012, -0.01 ), "ascending branch, pressure changing" );
  p.checkTangent( Vector3d( 0.29, 0.01, -0.02 ), "unloading" );
  p.update( Vector3d( 3.0, 0.03, 0.01 ) );
  p.checkTangent( Vector3d( 3.5, 0.02, 0.015 ), "softening with pressure" );
  p.checkTangent( Vector3d( -1.0, 0.03, 0.01 ), "reversed friction with pressure" );
  p.checkTangent( Vector3d( 3.5, 0.2, 0.1 ), "softening, normal cap" );
  p.update( Vector3d( 12.0, 0.15, -0.05 ) );
  p.checkTangent( Vector3d( 13.0, 0.16, -0.05 ), "residual, crushed" );
  p.checkTangent( Vector3d( 13.0, 0.10, -0.05 ), "residual, normal unloading" );
}

void testInvalidFrictionPropertiesThrow()
{
  const std::vector< double > good = { tauMax, 1.0, 2.0, 10.0, 0.4, 0.4 * tauMax, 200., Kn, 0.5, pMax };
  std::vector< std::pair< std::string, std::vector< double > > > cases = {
    { "too many properties", { tauMax, 1.0, 2.0, 10.0, 0.4, 0.4 * tauMax, 200., Kn, mu, pMax, 1.0 } },
  };
  auto with = [&]( int index, double value ) {
    auto p   = good;
    p[index] = value;
    return p;
  };
  cases.push_back( { "Kn = 0 with a normal cap", with( 7, 0.0 ) } );
  cases.push_back( { "mu < 0", with( 8, -0.1 ) } );
  for ( const auto& [what, badProps] : cases ) {
    bool threw = false;
    try {
      Marmot::Registry::BondSlipLawFactory::createBondSlipLaw( "MODELCODE2010BONDSLIP",
                                                               badProps.data(),
                                                               static_cast< int >( badProps.size() ),
                                                               1 );
    }
    catch ( const std::invalid_argument& ) {
      threw = true;
    }
    throwExceptionOnFailure( threw, what + " must be rejected" );
  }
}

/// alpha = 1 (accepted): no regularization, s0 = 0, the ascending branch is linear up to the peak
void testAlphaOne()
{
  const double p1[] = { tauMax, 1.0, 2.0, 10.0, 1.0, 0.4 * tauMax, 200., Kn };
  auto         law  = create( p1 );
  throwExceptionOnFailure( law->s0 == 0.0, "s0 = 0 for alpha = 1" );
  BondPoint p{ *law };
  p.update( Vector3d( 0.5, 0, 0 ) );
  throwExceptionOnFailure( checkIfEqual( p.stress( 0 ), 0.5 * tauMax, 1e-12 ), "linear ascending branch" );
  p.checkTangent( Vector3d( 0.6, 0, 0 ), "alpha = 1" );
}

/// Kn = 0 without a cap: no normal stress, no division by Kn, finite energy
void testZeroNormalStiffness()
{
  const double p0[] = { tauMax, 1.0, 2.0, 10.0, 0.4, 0.4 * tauMax, 200., 0.0, mu };
  auto         law  = create( p0, 9 );
  BondPoint    p{ *law };
  p.update( Vector3d( 3.0, 0.1, -0.2 ) );
  throwExceptionOnFailure( checkIfEqual( p.stress( 1 ), 0.0, 1e-14 ) && checkIfEqual( p.stress( 2 ), 0.0, 1e-14 ) &&
                             std::isfinite( p.energy ) && std::isfinite( p.dissipation ),
                           "Kn = 0" );
  throwExceptionOnFailure( checkIfEqual( p.stress( 0 ), law->envelope( 3.0 ).first, 1e-12 ),
                           "no friction without pressure" );
}

int main()
{
  executeTestsAndCollectExceptions( { testMonotonicEnvelope,
                                      testUnloadingAndReversal,
                                      testTangents,
                                      testInvalidPropertiesThrow,
                                      testOptionalFrictionPropertiesDefaultToZero,
                                      testCoulombEnvelopeAndCap,
                                      testTangentsWithFriction,
                                      testInvalidFrictionPropertiesThrow,
                                      testAlphaOne,
                                      testZeroNormalStiffness } );
  return 0;
}
