#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotKelvinChain.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotViscoelasticity.h"
#include "autodiff/forward/real.hpp"

using namespace Marmot::Testing;
using namespace Marmot::Materials::KelvinChain;
using namespace Marmot::ContinuumMechanics::Viscoelasticity;
using namespace Marmot;

void evaluateKCandUpdateStateVarsTestFunction()
{
  double factor = 1.35;

  // properties kelvin chain
  int    n       = 2;
  double min     = 10.;
  double spacing = 5.;

  // arbitrary moduli
  Properties elasticModuli( n );
  for ( int i = 0; i < n; i++ )
    elasticModuli( i ) = 3 * std::pow( 10, i );

  // time increment
  int dT = 10;

  // arbitrary initial state vars
  StateVarMatrix stateVars( 6, n );
  stateVars << 0.01, 0.1, 0.02, 0.2, 0.03, 0.3, 0.04, 0.4, 0.05, 0.5, 0.06, 0.6;

  // initialize compliance and strain
  double   uniaxialCompliance = 0;
  Vector6d dStrain            = { 0., 0., 0., 0., 0., 0. };

  // test retardation times
  Properties retardationTimes = generateRetardationTimes( n, min, spacing );
  for ( int i = 0; i < n; i++ )
    throwExceptionOnFailure( checkIfEqual( retardationTimes( i ), min * std::pow( spacing, i ) ),
                             MakeString() << __PRETTY_FUNCTION__ << " error in generation of retardation times" );

  double   corrCompliance = 0.16976016796969498;
  Vector6d corrDStrain    = { 0.033004975878657986,
                              0.06600995175731597,
                              0.09901492763597397,
                              0.13201990351463194,
                              0.1650248793932899,
                              0.1980298552719479 };

  // test evaluation of kelvin chain
  evaluateKelvinChain( dT, elasticModuli, retardationTimes, stateVars, uniaxialCompliance, dStrain, factor );
  throwExceptionOnFailure( checkIfEqual( uniaxialCompliance, corrCompliance ),
                           "error in uniaxial compliance of kelvin chain" );
  throwExceptionOnFailure( checkIfEqual< double >( dStrain, corrDStrain ),
                           MakeString() << __PRETTY_FUNCTION__ << " error in dStrain of kelvin chain" );

  // test state var update
  Vector6d       dStress              = { 0.1, 0.2, 0.3, 0.06, 0.04, 0.02 }; // arbitrary stress increment
  double         nu                   = 0.2;
  Matrix6d       unitComplianceMatrix = ContinuumMechanics::Elasticity::Isotropic::complianceTensor( 1.0, nu );
  StateVarMatrix corrStateVars( 6, n );
  corrStateVars << 0.0036787944117144234, 0.0818730753077982, 0.03264241117657116, 0.16737153555403675,
    0.061606027941427874, 0.2528699958002753, 0.04505696447062846, 0.3318427631573212, 0.03862182994108596,
    0.4122656844897432, 0.03218669541154347, 0.49268860582216517;

  updateStateVarMatrix( dT, elasticModuli, retardationTimes, stateVars, dStress, unitComplianceMatrix );
  throwExceptionOnFailure( checkIfEqual< double >( stateVars, corrStateVars ),
                           MakeString() << __PRETTY_FUNCTION__ << " error in state var update" );
}

void computeLambdaAndBetaTestFunction()
{
  double lambda, beta;
  double dT = 30;

  // case dT_tau >= 30.0
  double tau = 1 / ( 30 / dT );
  computeLambdaAndBeta( dT, tau, lambda, beta );

  throwExceptionOnFailure( checkIfEqual( beta, 0 ),
                           MakeString() << __PRETTY_FUNCTION__ << " error in beta with dT_tau >= 30.0" );
  throwExceptionOnFailure( checkIfEqual( lambda, tau / dT ),
                           MakeString() << __PRETTY_FUNCTION__ << " error in lambda with dT_tau >= 30.0" );

  // case dT_tau < 1e-6
  tau = 1 / ( 1e-7 / dT );
  computeLambdaAndBeta( dT, tau, lambda, beta );

  double corrLam = 1 - 0.5 * dT / tau + 1. / 6 * dT / tau * dT / tau;

  throwExceptionOnFailure( checkIfEqual( beta, 1 ),
                           MakeString() << __PRETTY_FUNCTION__ << " error in beta with dT_tau < 1e-6" );
  throwExceptionOnFailure( checkIfEqual( lambda, corrLam ),
                           MakeString() << __PRETTY_FUNCTION__ << " error in lambda with dT_tau < 1e-6" );

  // case else
  tau = 1 / ( 10 / dT );
  computeLambdaAndBeta( dT, tau, lambda, beta );

  double corrBeta = std::exp( -dT / tau );
  corrLam         = ( 1 - beta ) * ( tau / dT );

  throwExceptionOnFailure( checkIfEqual( beta, corrBeta ),
                           MakeString() << __PRETTY_FUNCTION__ << " error in beta with 1e-6 <= dT_tau < 30.0" );
  throwExceptionOnFailure( checkIfEqual( lambda, corrLam ),
                           MakeString() << __PRETTY_FUNCTION__ << " error in lambda with 1e-6 <= dT_tau < 30.0" );
}

void computeElasticModuliTestFunction()
{
  // approximation order of Post-Widder-Formula
  static constexpr int maxDerivativeOrder = 3;
  // Parameters for log-power-law from ModelCode
  double m = 1;
  double n = 1;
  // compliance function
  auto phi = [&]( autodiff::Real< maxDerivativeOrder, double > t ) {
    return ComplianceFunctions::logPowerLaw( t, m, n );
  };
  // Properties retardationTimes
  int        nRetardationTimes = 6;
  int        minTau            = 7;
  int        nTau              = 7;
  Properties retardationTimes  = generateRetardationTimes( nRetardationTimes, minTau, nTau );
  // using Gauss integration
  bool       gaussQuadrature    = true;
  Properties elasticModuliGauss = computeElasticModuli< maxDerivativeOrder >( phi, retardationTimes, gaussQuadrature );
  // values from python script
  Properties corrModuliGauss( nRetardationTimes );
  corrModuliGauss << 0.601432770397282, 0.526130944473972, 0.515640199270429, 0.514147063269368, 0.513933871562560,
    0.513903417920470;

  throwExceptionOnFailure( checkIfEqual< double >( elasticModuliGauss, corrModuliGauss ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " error in computation of elastic moduli using Gauss integration" );

  // using midpoint rule
  Properties elasticModuli = computeElasticModuli< maxDerivativeOrder >( phi, retardationTimes );
  // values from python script
  Properties corrModuli( nRetardationTimes );
  corrModuli << 0.590863788959411, 0.524457570465417, 0.515398044758814, 0.514112407213863, 0.513928919423037,
    0.513902710445964;

  throwExceptionOnFailure( checkIfEqual< double >( elasticModuli, corrModuli ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " error in computation of elastic moduli using midpoint rule" );
}

void evaluatePostWidderFormulaTestFunction()
{
  // approximation order of Post-Widder-Formula
  static constexpr int maxDerivativeOrder = 3;
  // Parameters for log-power-law from ModelCode
  double m = 1;
  double n = 1;
  // compliance function
  auto phi = [&]( autodiff::Real< maxDerivativeOrder, double > t ) {
    return ComplianceFunctions::logPowerLaw( t, m, n );
  };
  // retardation time
  double tau = 7;
  // evaluation of Post-Widder-Formula
  double pw = evaluatePostWidderFormula< maxDerivativeOrder >( phi, tau );
  // value from python script
  double correctValue = 0.8697407963936890;

  throwExceptionOnFailure( checkIfEqual( pw, correctValue ),
                           MakeString() << __PRETTY_FUNCTION__ << " error in Post-Widder-Formula" );
}

void approximateZerothComplianceTestFunction()
{ // approximation order of Post-Widder-Formula
  static constexpr int maxDerivativeOrder = 3;
  // Parameters for log-power-law from ModelCode
  double m = 1;
  double n = 1;
  // compliance function
  auto phi = [&]( autodiff::Real< maxDerivativeOrder, double > t ) {
    return ComplianceFunctions::logPowerLaw( t, m, n );
  };
  // min retardation time
  double tauMin = 7;
  // zeroth compliance
  double E0inv = approximateZerothCompliance< maxDerivativeOrder >( phi, tauMin );
  // value from python script ( with symbolic integration -> tol = 1e-6 )
  double corrZerothCompliance = 0.7866890;

  throwExceptionOnFailure( checkIfEqual( E0inv, corrZerothCompliance, 1e-6 ),
                           MakeString() << __PRETTY_FUNCTION__ << " error in zeroth compliance" );
}

int main()
{

  auto tests = std::vector< std::function< void() > >{ evaluateKCandUpdateStateVarsTestFunction,
                                                       computeLambdaAndBetaTestFunction,
                                                       computeElasticModuliTestFunction,
                                                       approximateZerothComplianceTestFunction,
                                                       evaluatePostWidderFormulaTestFunction };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
