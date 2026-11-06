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
  // MC2010, good bond conditions, unconfined, pull-out failure, fck = 30 MPa (units: N, mm)
  const double tauMax  = 2.5 * std::sqrt( 30. );
  const double props[] = { tauMax, 1.0, 2.0, 10.0, 0.4, 0.4 * tauMax, 200., 1000. };

  std::unique_ptr< ModelCode2010BondSlip > create()
  {
    return std::unique_ptr< ModelCode2010BondSlip >( dynamic_cast< ModelCode2010BondSlip* >(
      Marmot::Registry::BondSlipLawFactory::createBondSlipLaw( "MODELCODE2010BONDSLIP", props, 8, 1 ) ) );
  }

  /// a bond point with its own state
  struct BondPoint {
    const ModelCode2010BondSlip& law;
    std::vector< double >        stateVars = std::vector< double >( 2, 0.0 );
    Vector3d                     stress    = Vector3d::Zero();
    double                       energy = 0.0, dissipation = 0.0;
    Vector3d                     slip = Vector3d::Zero();
    Matrix3d                     C    = Matrix3d::Zero();

    /// update to the slip s, the state is committed
    void update( const Vector3d& s )
    {
      Marmot::Materials::BondSlipLaw::State state{ stress, energy, dissipation, stateVars.data() };
      law.computeBondStress( state, C, s, s - slip, { 0.0, 1.0 } );
      stress      = state.bondStress;
      energy      = state.elasticEnergyDensity;
      dissipation = state.dissipation;
      slip        = s;
    }

    /// trial update without committing, for numerical differentiation
    Vector3d trial( const Vector3d& s ) const
    {
      auto                                  sv = stateVars;
      Marmot::Materials::BondSlipLaw::State state{ stress, energy, dissipation, sv.data() };
      Matrix3d                              C_;
      law.computeBondStress( state, C_, s, s - slip, { 0.0, 1.0 } );
      return state.bondStress;
    }

    /// checks the tangent of an update from the current state to s
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
                               Marmot::Utilities::MakeString()
                                 << "tangent at " << where << ": " << C_( 0, 0 ) << " vs " << numC( 0, 0 ) );
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
  for ( const auto& [what, props] : cases ) {
    bool threw = false;
    try {
      Marmot::Registry::BondSlipLawFactory::createBondSlipLaw( "MODELCODE2010BONDSLIP",
                                                               props.data(),
                                                               static_cast< int >( props.size() ),
                                                               1 );
    }
    catch ( const std::invalid_argument& ) {
      threw = true;
    }
    throwExceptionOnFailure( threw, what + " must be rejected" );
  }
}

int main()
{
  executeTestsAndCollectExceptions(
    { testMonotonicEnvelope, testUnloadingAndReversal, testTangents, testInvalidPropertiesThrow } );
  return 0;
}
