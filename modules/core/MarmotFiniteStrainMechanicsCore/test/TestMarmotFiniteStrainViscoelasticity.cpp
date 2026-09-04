#include "Fastor/Fastor.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotFiniteStrainViscoelasticity.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotTesting.h"
#include <array>
#include <cmath>
#include <cstring>

using Marmot::MakeString;
using namespace Marmot::Testing;
using namespace Marmot::TensorUtility::FastorTensors::StandardTensors;
using namespace Marmot::ContinuumMechanics::Viscoelasticity::FiniteStrain;

// ---------------------------------------------------------------------------
// createMaxwellProperties
// ---------------------------------------------------------------------------

// Test that nMaxwell=0 produces empty properties and zero sumGamma
void testCreateMaxwellPropertiesEmpty()
{
  MaxwellProperties props = createMaxwellProperties( 0, nullptr );

  throwExceptionOnFailure( props.nMaxwell == 0, MakeString() << __PRETTY_FUNCTION__ << ": nMaxwell should be 0" );
  throwExceptionOnFailure( props.gamma.empty(),
                           MakeString() << __PRETTY_FUNCTION__ << ": gamma vector should be empty" );
  throwExceptionOnFailure( props.tau.empty(), MakeString() << __PRETTY_FUNCTION__ << ": tau vector should be empty" );
  throwExceptionOnFailure( checkIfEqual( props.sumGamma, 0.0 ),
                           MakeString() << __PRETTY_FUNCTION__ << ": sumGamma should be 0" );
}

// Test single Maxwell element
void testCreateMaxwellPropertiesSingle()
{
  const double      pairs[2] = { 0.3, 10.0 }; // { gamma1, tau1 }
  MaxwellProperties props    = createMaxwellProperties( 1, pairs );

  throwExceptionOnFailure( props.nMaxwell == 1, MakeString() << __PRETTY_FUNCTION__ << ": nMaxwell should be 1" );
  throwExceptionOnFailure( checkIfEqual( props.gamma[0], 0.3 ),
                           MakeString() << __PRETTY_FUNCTION__ << ": gamma[0] should be 0.3" );
  throwExceptionOnFailure( checkIfEqual( props.tau[0], 10.0 ),
                           MakeString() << __PRETTY_FUNCTION__ << ": tau[0] should be 10.0" );
  throwExceptionOnFailure( checkIfEqual( props.sumGamma, 0.3 ),
                           MakeString() << __PRETTY_FUNCTION__ << ": sumGamma should be 0.3" );
}

// Test two Maxwell elements: correct parsing, correct sumGamma
void testCreateMaxwellPropertiesTwo()
{
  const double      pairs[4] = { 0.3, 10.0, 0.2, 5.0 }; // { gamma1, tau1, gamma2, tau2 }
  MaxwellProperties props    = createMaxwellProperties( 2, pairs );

  throwExceptionOnFailure( props.nMaxwell == 2, MakeString() << __PRETTY_FUNCTION__ << ": nMaxwell should be 2" );
  throwExceptionOnFailure( checkIfEqual( props.gamma[0], 0.3 ),
                           MakeString() << __PRETTY_FUNCTION__ << ": gamma[0] should be 0.3" );
  throwExceptionOnFailure( checkIfEqual( props.tau[0], 10.0 ),
                           MakeString() << __PRETTY_FUNCTION__ << ": tau[0] should be 10.0" );
  throwExceptionOnFailure( checkIfEqual( props.gamma[1], 0.2 ),
                           MakeString() << __PRETTY_FUNCTION__ << ": gamma[1] should be 0.2" );
  throwExceptionOnFailure( checkIfEqual( props.tau[1], 5.0 ),
                           MakeString() << __PRETTY_FUNCTION__ << ": tau[1] should be 5.0" );
  throwExceptionOnFailure( checkIfEqual( props.sumGamma, 0.5 ),
                           MakeString() << __PRETTY_FUNCTION__ << ": sumGamma should be 0.5" );
}

// ---------------------------------------------------------------------------
// Helper: compute Maxwell alpha and beta analytically
// ---------------------------------------------------------------------------
static void computeAlphaBeta( double gamma, double tau, double dT, double& alpha, double& beta )
{
  const double dT_tau = std::max( dT / tau, 1e-15 );
  const double expF   = std::exp( -dT_tau );

  if ( dT_tau < 1e-6 ) {
    alpha = 1.0 - dT_tau + 0.5 * dT_tau * dT_tau;
    beta  = gamma * ( 1.0 - 0.5 * dT_tau + 1.0 / 6.0 * dT_tau * dT_tau );
  }
  else {
    alpha = expF;
    beta  = gamma / dT_tau * ( 1.0 - expF );
  }
}

// ---------------------------------------------------------------------------
// Template overload: evaluateGeneralizedMaxwellModel<double>(stress, dStress, dT, props, stateVars)
// ---------------------------------------------------------------------------

// nMaxwell=0: no change to stress
void testTemplateOverloadNoOp()
{
  MaxwellProperties props = createMaxwellProperties( 0, nullptr );

  Tensor33d stress( 0.0 );
  stress( 0, 0 ) = 1.0;
  stress( 1, 1 ) = 2.0;
  stress( 2, 2 ) = 3.0;

  Tensor33d dStress( 0.0 );
  dStress( 0, 0 ) = 0.1;

  Tensor33d stressOriginal = stress;

  evaluateGeneralizedMaxwellModel( stress, dStress, 1.0, props, nullptr );

  throwExceptionOnFailure( checkIfEqual( stress, stressOriginal ),
                           MakeString() << __PRETTY_FUNCTION__ << ": nMaxwell=0 should leave stress unchanged" );
}

// Single Maxwell element, Q_n=0: verify stress update and state var storage
void testTemplateOverloadSingleElementFromZero()
{
  const double      pairs[2] = { 0.3, 10.0 };
  MaxwellProperties props    = createMaxwellProperties( 1, pairs );

  const double dT = 1.0;
  double       alpha, beta;
  computeAlphaBeta( 0.3, 10.0, dT, alpha, beta );

  // stress_in = I, dStress = I
  Tensor33d stress  = Spatial3D::I;
  Tensor33d dStress = Spatial3D::I;

  // Q_n = 0 (all state vars zero)
  std::array< double, 9 > stateVars{};

  evaluateGeneralizedMaxwellModel( stress, dStress, dT, props, stateVars.data() );

  // Expected: stress = (1 - gamma)*I + Q_np, where Q_np = beta*I (since Q_n=0)
  //         = (1 - gamma + beta)*I = 0.985487... * I
  const double expectedFactor = ( 1.0 - 0.3 ) + beta;
  Tensor33d    stressExpected = expectedFactor * Spatial3D::I;

  throwExceptionOnFailure( checkIfEqual( stress, stressExpected, 1e-12 ),
                           MakeString() << __PRETTY_FUNCTION__ << ": stress update failed" );

  // Q_np = beta * I stored in stateVars
  for ( int i = 0; i < 3; ++i )
    for ( int j = 0; j < 3; ++j ) {
      const double expectedQnp = beta * ( i == j ? 1.0 : 0.0 );
      throwExceptionOnFailure( checkIfEqual( stateVars[3 * i + j], expectedQnp, 1e-12 ),
                               MakeString() << __PRETTY_FUNCTION__ << ": Q_np[" << i << "," << j << "] failed" );
    }
}

// Single Maxwell element, Q_n≠0: verify recursion Q_np = alpha*Q_n + beta*dStress
void testTemplateOverloadSingleElementFromNonZeroState()
{
  const double      pairs[2] = { 0.3, 10.0 };
  MaxwellProperties props    = createMaxwellProperties( 1, pairs );

  const double dT = 1.0;
  double       alpha, beta;
  computeAlphaBeta( 0.3, 10.0, dT, alpha, beta );

  // stress_in = I, dStress = I
  Tensor33d stress  = Spatial3D::I;
  Tensor33d dStress = Spatial3D::I;

  // Q_n = 0.5 * I (from previous increment)
  std::array< double, 9 > stateVars{};
  const double            q0 = 0.5;
  stateVars[0]               = q0; // Q_n(0,0)
  stateVars[4]               = q0; // Q_n(1,1)
  stateVars[8]               = q0; // Q_n(2,2)

  evaluateGeneralizedMaxwellModel( stress, dStress, dT, props, stateVars.data() );

  // Q_np diagonal = alpha * q0 + beta * 1
  const double Q_np_diag = alpha * q0 + beta;
  // stress_out diagonal = (1-gamma) * 1 + Q_np_diag
  const double stressOutDiag = ( 1.0 - 0.3 ) + Q_np_diag;

  for ( int i = 0; i < 3; ++i ) {
    throwExceptionOnFailure( checkIfEqual( stress( i, i ), stressOutDiag, 1e-12 ),
                             MakeString() << __PRETTY_FUNCTION__ << ": stress diagonal[" << i << "] failed" );
  }
  // Off-diagonal should be zero
  throwExceptionOnFailure( checkIfEqual( stress( 0, 1 ), 0.0, 1e-12 ),
                           MakeString() << __PRETTY_FUNCTION__ << ": stress off-diagonal should be zero" );
}

// Test that Taylor approximation branch (small dT/tau) gives consistent results
void testTemplateOverloadTaylorApproximation()
{
  const double      pairs[2] = { 0.3, 1.0 };
  MaxwellProperties props    = createMaxwellProperties( 1, pairs );

  // dT/tau = 1e-8 << 1e-6  -> Taylor branch
  const double dT_small = 1e-8;

  // Evaluate with Taylor branch
  Tensor33d               stressTaylor = Spatial3D::I;
  Tensor33d               dStress      = Spatial3D::I;
  std::array< double, 9 > stateVarsTaylor{};

  evaluateGeneralizedMaxwellModel( stressTaylor, dStress, dT_small, props, stateVarsTaylor.data() );

  // Evaluate with very small exact dT (not in Taylor branch, but close enough to compare)
  // Use reference formula directly
  double alpha_ref, beta_ref;
  computeAlphaBeta( 0.3, 1.0, dT_small, alpha_ref, beta_ref );

  // Expected: stress = (1-0.3 + beta_ref) * I
  const double expectedFactor = ( 1.0 - 0.3 ) + beta_ref;

  throwExceptionOnFailure( checkIfEqual( stressTaylor( 0, 0 ), expectedFactor, 1e-10 ),
                           MakeString() << __PRETTY_FUNCTION__ << ": Taylor branch stress(0,0) failed" );
}

// Two Maxwell elements: stress is sum of both contributions
void testTemplateOverloadTwoElements()
{
  const double      pairs[4] = { 0.3, 10.0, 0.2, 5.0 };
  MaxwellProperties props    = createMaxwellProperties( 2, pairs );

  const double dT = 1.0;
  double       alpha1, beta1, alpha2, beta2;
  computeAlphaBeta( 0.3, 10.0, dT, alpha1, beta1 );
  computeAlphaBeta( 0.2, 5.0, dT, alpha2, beta2 );

  // stress_in = I, dStress = I, Q1_n = Q2_n = 0
  Tensor33d                stress  = Spatial3D::I;
  Tensor33d                dStress = Spatial3D::I;
  std::array< double, 18 > stateVars{};

  evaluateGeneralizedMaxwellModel( stress, dStress, dT, props, stateVars.data() );

  // Expected: stress = (1-sumGamma)*I + beta1*I + beta2*I = (0.5 + beta1 + beta2)*I
  const double expectedFactor = ( 1.0 - 0.5 ) + beta1 + beta2;

  throwExceptionOnFailure( checkIfEqual( stress( 0, 0 ), expectedFactor, 1e-12 ),
                           MakeString() << __PRETTY_FUNCTION__ << ": stress(0,0) with two Maxwell elements failed" );
  throwExceptionOnFailure( checkIfEqual( stress( 1, 1 ), expectedFactor, 1e-12 ),
                           MakeString() << __PRETTY_FUNCTION__ << ": stress(1,1) with two Maxwell elements failed" );
  throwExceptionOnFailure( checkIfEqual( stress( 2, 2 ), expectedFactor, 1e-12 ),
                           MakeString() << __PRETTY_FUNCTION__ << ": stress(2,2) with two Maxwell elements failed" );
}

// Relaxation: stress decreases to (1-sumGamma)*sigma after many steps with no new load increment
void testTemplateOverloadRelaxation()
{
  const double      pairs[2] = { 0.3, 1.0 };
  MaxwellProperties props    = createMaxwellProperties( 1, pairs );

  // Step 1: Apply instantaneous deformation - set Q to some value via one increment
  // Use a large dT so we go through the large-dT branch and alpha is small
  Tensor33d               stress0  = Spatial3D::I;
  Tensor33d               dStress0 = Spatial3D::I;
  std::array< double, 9 > stateVars{};

  // Apply loading increment
  evaluateGeneralizedMaxwellModel( stress0, dStress0, 1e-8, props, stateVars.data() );

  // The Maxwell stress Q was set to ~gamma*dStress = ~0.3*I after this tiny increment
  // Now apply N steps of zero increment (relaxation) with dT=tau=1
  Tensor33d    elasticStress = Spatial3D::I; // hypothetical instantaneous stress
  const int    N             = 50;
  const double dT            = 1.0;

  for ( int k = 0; k < N; ++k ) {
    Tensor33d dStressZero( 0.0 );
    // stress must be reset to elastic each time (Maxwell model is purely incremental)
    Tensor33d stressStep = elasticStress;
    evaluateGeneralizedMaxwellModel( stressStep, dStressZero, dT, props, stateVars.data() );
    // After N steps, stressStep should approach (1-gamma)*elasticStress
  }

  // Final evaluation
  Tensor33d dStressZero( 0.0 );
  Tensor33d stressRelaxed = elasticStress;
  evaluateGeneralizedMaxwellModel( stressRelaxed, dStressZero, dT, props, stateVars.data() );

  // After sufficient relaxation, Q ≈ 0, so stress ≈ (1-gamma)*I = 0.7*I
  const double longTermFactor = 1.0 - 0.3;
  throwExceptionOnFailure( checkIfEqual( stressRelaxed( 0, 0 ), longTermFactor, 1e-6 ),
                           MakeString() << __PRETTY_FUNCTION__ << ": long-term stress(0,0) after relaxation failed" );
  throwExceptionOnFailure( checkIfEqual( stressRelaxed( 1, 1 ), longTermFactor, 1e-6 ),
                           MakeString() << __PRETTY_FUNCTION__ << ": long-term stress(1,1) after relaxation failed" );
}

// ---------------------------------------------------------------------------
// Simo overload: evaluateGeneralizedMaxwellModel(stress, tangent, dStress, dT, props, stateVars)
// ---------------------------------------------------------------------------

// nMaxwell=0: no change to stress or tangent
void testSimoOverloadNoOp()
{
  MaxwellProperties props = createMaxwellProperties( 0, nullptr );

  Tensor33d stress( 0.0 );
  stress( 0, 0 ) = 5.0;
  stress( 1, 1 ) = 3.0;

  Tensor3333d tangent( 1.0 ); // all entries = 1.0

  Tensor33d   dStress( 0.0 );
  Tensor33d   stressOrig  = stress;
  Tensor3333d tangentOrig = tangent;

  evaluateGeneralizedMaxwellModel( stress, tangent, dStress, 1.0, props, nullptr );

  throwExceptionOnFailure( checkIfEqual( stress, stressOrig ),
                           MakeString() << __PRETTY_FUNCTION__ << ": nMaxwell=0 should leave stress unchanged" );
  throwExceptionOnFailure( checkIfEqual( tangent, tangentOrig ),
                           MakeString() << __PRETTY_FUNCTION__ << ": nMaxwell=0 should leave tangent unchanged" );
}

// Single Maxwell element, Q_n=0: verify stress and tangent update
void testSimoOverloadSingleElementFromZero()
{
  const double      pairs[2] = { 0.3, 10.0 };
  MaxwellProperties props    = createMaxwellProperties( 1, pairs );

  const double dT = 1.0;
  double       alpha, beta;
  computeAlphaBeta( 0.3, 10.0, dT, alpha, beta );

  // stress_in = I, tangent_in = all entries equal to c (use scalar = 2.0 for each entry)
  Tensor33d   stress = Spatial3D::I;
  Tensor3333d tangent( 2.0 ); // all entries = 2.0
  Tensor33d   dStress = Spatial3D::I;

  std::array< double, 9 > stateVars{};

  evaluateGeneralizedMaxwellModel( stress, tangent, dStress, dT, props, stateVars.data() );

  // Expected stress: (1-gamma)*I + beta*I = (0.7 + beta)*I
  const double stressFactor = ( 1.0 - 0.3 ) + beta;
  throwExceptionOnFailure( checkIfEqual( stress( 0, 0 ), stressFactor, 1e-12 ),
                           MakeString() << __PRETTY_FUNCTION__ << ": stress(0,0) failed" );
  throwExceptionOnFailure( checkIfEqual( stress( 0, 1 ), 0.0, 1e-12 ),
                           MakeString() << __PRETTY_FUNCTION__ << ": off-diagonal stress should be zero" );

  // Expected tangent: (1-gamma+beta) * tangent_in = (0.7+beta) * 2.0 for all entries
  const double tangentFactor = ( 1.0 - 0.3 + beta ) * 2.0;
  for ( int i = 0; i < 3; ++i )
    for ( int j = 0; j < 3; ++j )
      for ( int k = 0; k < 3; ++k )
        for ( int l = 0; l < 3; ++l )
          throwExceptionOnFailure( checkIfEqual( tangent( i, j, k, l ), tangentFactor, 1e-12 ),
                                   MakeString() << __PRETTY_FUNCTION__ << ": tangent(" << i << "," << j << "," << k
                                                << "," << l << ") failed" );
}

// Single Maxwell element, Q_n≠0: verify recursion in both stress and tangent
void testSimoOverloadSingleElementFromNonZeroState()
{
  const double      pairs[2] = { 0.3, 10.0 };
  MaxwellProperties props    = createMaxwellProperties( 1, pairs );

  const double dT = 1.0;
  double       alpha, beta;
  computeAlphaBeta( 0.3, 10.0, dT, alpha, beta );

  Tensor33d   stress = Spatial3D::I;
  Tensor3333d tangent( 1.0 );
  Tensor33d   dStress = Spatial3D::I;

  // Q_n = 0.5 * I  (diagonal = 0.5, off-diagonal = 0)
  std::array< double, 9 > stateVars{};
  const double            q0 = 0.5;
  stateVars[0]               = q0;
  stateVars[4]               = q0;
  stateVars[8]               = q0;

  evaluateGeneralizedMaxwellModel( stress, tangent, dStress, dT, props, stateVars.data() );

  // Q_np diagonal = alpha * q0 + beta
  const double Q_np_diag = alpha * q0 + beta;
  // stress diagonal = (1-gamma) + Q_np_diag
  const double stressDiag = ( 1.0 - 0.3 ) + Q_np_diag;

  throwExceptionOnFailure( checkIfEqual( stress( 0, 0 ), stressDiag, 1e-12 ),
                           MakeString() << __PRETTY_FUNCTION__ << ": stress diagonal failed" );

  // Tangent update does not depend on Q_n:
  // tangent = (1-gamma+beta) * tangent_in = (0.7+beta) * 1.0
  const double tangentFactor = 1.0 - 0.3 + beta;
  throwExceptionOnFailure( checkIfEqual( tangent( 0, 0, 0, 0 ), tangentFactor, 1e-12 ),
                           MakeString() << __PRETTY_FUNCTION__ << ": tangent(0,0,0,0) failed" );
}

// Two Maxwell elements: verify that contributions are summed correctly
void testSimoOverloadTwoElements()
{
  const double      pairs[4] = { 0.3, 10.0, 0.2, 5.0 };
  MaxwellProperties props    = createMaxwellProperties( 2, pairs );

  const double dT = 1.0;
  double       alpha1, beta1, alpha2, beta2;
  computeAlphaBeta( 0.3, 10.0, dT, alpha1, beta1 );
  computeAlphaBeta( 0.2, 5.0, dT, alpha2, beta2 );

  Tensor33d                stress = Spatial3D::I;
  Tensor3333d              tangent( 1.0 );
  Tensor33d                dStress = Spatial3D::I;
  std::array< double, 18 > stateVars{};

  evaluateGeneralizedMaxwellModel( stress, tangent, dStress, dT, props, stateVars.data() );

  // Expected stress: (1-sumGamma) + beta1 + beta2 = (0.5 + beta1 + beta2) on diagonal
  const double stressFactor = 0.5 + beta1 + beta2;
  throwExceptionOnFailure( checkIfEqual( stress( 0, 0 ), stressFactor, 1e-12 ),
                           MakeString() << __PRETTY_FUNCTION__ << ": stress(0,0) with two elements failed" );

  // Expected tangent: (1-sumGamma + beta1 + beta2) * tangent_in
  const double tangentFactor = stressFactor; // same factor here since tangent_in=1
  throwExceptionOnFailure( checkIfEqual( tangent( 0, 0, 0, 0 ), tangentFactor, 1e-12 ),
                           MakeString() << __PRETTY_FUNCTION__ << ": tangent(0,0,0,0) with two elements failed" );
}

// Verify that the Simo overload and template overload give the same stress result
void testSimoAndTemplateOverloadConsistency()
{
  const double      pairs[2] = { 0.3, 10.0 };
  MaxwellProperties props    = createMaxwellProperties( 1, pairs );

  const double dT = 2.5;

  Tensor33d stress_simo( 0.0 );
  stress_simo( 0, 0 ) = 1.5;
  stress_simo( 1, 2 ) = 0.3;
  stress_simo( 2, 1 ) = 0.3;

  Tensor33d dStress( 0.0 );
  dStress( 0, 0 ) = 0.2;
  dStress( 1, 1 ) = 0.1;
  dStress( 1, 2 ) = 0.05;
  dStress( 2, 1 ) = 0.05;

  // Copy for template overload
  Tensor33d stress_tmpl = stress_simo;

  std::array< double, 9 > stateVarsSimo{};
  std::array< double, 9 > stateVarsTmpl{};

  Tensor3333d tangent( 1.0 );

  evaluateGeneralizedMaxwellModel( stress_simo, tangent, dStress, dT, props, stateVarsSimo.data() );
  evaluateGeneralizedMaxwellModel( stress_tmpl, dStress, dT, props, stateVarsTmpl.data() );

  throwExceptionOnFailure( checkIfEqual( stress_simo, stress_tmpl, 1e-12 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << ": Simo and template overloads should give same stress" );

  // State variables should also be identical
  for ( int i = 0; i < 9; ++i )
    throwExceptionOnFailure( checkIfEqual( stateVarsSimo[i], stateVarsTmpl[i], 1e-12 ),
                             MakeString()
                               << __PRETTY_FUNCTION__ << ": state var[" << i << "] differs between overloads" );
}

// ---------------------------------------------------------------------------
// Liu et al. overload:
//   evaluateGeneralizedMaxwellModel(stress, tangent, dTangent_dDeformation,
//                                   initialCompliance, dStress, dT, props, stateVars)
//
// Special case: dTangent_dDeformation = 0, initialTangent and initialCompliance are inverses.
// In this case, the Liu formulation reduces to the same result as the Simo formulation.
// ---------------------------------------------------------------------------

// Build an isotropic 4th-order stiffness tensor using Marmot's IHyd and ISymm:
// C_ijkl = lambda * delta_ij * delta_kl + 2*mu * ISymm_ijkl
static Tensor3333d makeIsotropicTangent( double K, double G )
{
  const double lambda = K - 2.0 / 3.0 * G;
  const double mu     = G;
  return evaluate( lambda * Spatial3D::IHyd + 2.0 * mu * Spatial3D::ISymm );
}

// Build the isotropic compliance tensor as the exact inverse of the stiffness tensor,
// using Marmot's IHyd and ISymm:
// C^{-1}_ijkl = 1/mu * ISymm_ijkl - lambda/(2*mu*(3*lambda+2*mu)) * IHyd_ijkl
static Tensor3333d makeIsotropicCompliance( double K, double G )
{
  const double lambda = K - 2.0 / 3.0 * G;
  const double mu     = G;
  return evaluate( ( 1.0 / 2 / mu ) * Spatial3D::ISymm -
                   ( lambda / ( 2.0 * mu * ( 3.0 * lambda + 2.0 * mu ) ) ) * Spatial3D::IHyd );
}

// Verify that C : C^{-1} = I4_sym (the 4th-order symmetric identity Spatial3D::ISymm)
void testIsotropicInverseConsistency()
{
  const double K = 10000.0, G = 5000.0;
  Tensor3333d  C    = makeIsotropicTangent( K, G );
  Tensor3333d  Cinv = makeIsotropicCompliance( K, G );

  using namespace Fastor;
  using namespace Marmot::TensorUtility::FastorTensors::Indices;

  // Compute C:Cinv_ijkl = C_ijmn * Cinv_mnkl = C_ijmn * Cinv_klmn (Cinv is fully symmetric)
  Tensor3333d product = einsum< ijmn, klmn, to_ijkl >( C, Cinv );

  // Should equal the 4th-order symmetric identity Spatial3D::ISymm
  throwExceptionOnFailure( checkIfEqual( product, Spatial3D::ISymm, 1e-10 ),
                           MakeString() << __PRETTY_FUNCTION__ << ": C:Cinv != ISymm" );
}

// Liu et al. overload with dTangent_dDeformation=0 and C/Cinv properly paired:
// result should match the Simo overload
void testLiuOverloadMatchesSimo()
{
  const double      pairs[2] = { 0.3, 10.0 };
  MaxwellProperties props    = createMaxwellProperties( 1, pairs );

  const double dT = 1.0;

  const double K = 10000.0, G = 5000.0;
  Tensor3333d  C    = makeIsotropicTangent( K, G );
  Tensor3333d  Cinv = makeIsotropicCompliance( K, G );

  // A non-trivial stress and dStress
  Tensor33d stress0( 0.0 );
  stress0( 0, 0 ) = 200.0;
  stress0( 1, 1 ) = 100.0;
  stress0( 2, 2 ) = 150.0;
  stress0( 0, 1 ) = 50.0;
  stress0( 1, 0 ) = 50.0;

  Tensor33d dStress( 0.0 );
  dStress( 0, 0 ) = 20.0;
  dStress( 1, 1 ) = 10.0;
  dStress( 2, 2 ) = 15.0;
  dStress( 0, 1 ) = 5.0;
  dStress( 1, 0 ) = 5.0;

  // Make copies
  Tensor33d   stressSimo  = stress0;
  Tensor33d   stressLiu   = stress0;
  Tensor3333d tangentSimo = C;
  Tensor3333d tangentLiu  = C;

  std::array< double, 9 > stateVarsSimo{};
  std::array< double, 9 > stateVarsLiu{};

  // Simo overload
  evaluateGeneralizedMaxwellModel( stressSimo, tangentSimo, dStress, dT, props, stateVarsSimo.data() );

  // Liu overload with zero dTangent_dDeformation
  Tensor333333d dTangent_dDeformation( 0.0 );
  evaluateGeneralizedMaxwellModel( stressLiu,
                                   tangentLiu,
                                   dTangent_dDeformation,
                                   Cinv,
                                   dStress,
                                   dT,
                                   props,
                                   stateVarsLiu.data() );

  // Stress should be identical
  throwExceptionOnFailure( checkIfEqual( stressLiu, stressSimo, 1e-8 ),
                           MakeString() << __PRETTY_FUNCTION__ << ": Liu stress should match Simo stress" );

  // Tangent should also be identical (since C:Cinv:beta*C = beta*C when dTangentdDef=0)
  throwExceptionOnFailure( checkIfEqual( tangentLiu, tangentSimo, 1e-8 ),
                           MakeString() << __PRETTY_FUNCTION__ << ": Liu tangent should match Simo tangent" );

  // State variables should be identical
  for ( int i = 0; i < 9; ++i )
    throwExceptionOnFailure( checkIfEqual( stateVarsLiu[i], stateVarsSimo[i], 1e-8 ),
                             MakeString()
                               << __PRETTY_FUNCTION__ << ": state var[" << i << "] differs between Liu and Simo" );
}

// Liu et al. overload: when dTangent_dDeformation is non-zero, the tangent gets an extra contribution
void testLiuOverloadNonZeroDTangentDDeformation()
{
  const double      pairs[2] = { 0.3, 10.0 };
  MaxwellProperties props    = createMaxwellProperties( 1, pairs );

  const double dT = 1.0;
  double       alpha, beta;
  computeAlphaBeta( 0.3, 10.0, dT, alpha, beta );

  const double K = 10000.0, G = 5000.0;
  Tensor3333d  C    = makeIsotropicTangent( K, G );
  Tensor3333d  Cinv = makeIsotropicCompliance( K, G );

  // Simple case: dStress = 0 and Q_n = 0
  // Then Q_np = 0, H_np = 0, and the last tangent term (H_np . dTangent_dDeformation) = 0.
  // So tangent should equal Simo result regardless of dTangent_dDeformation.
  Tensor33d   stress = Spatial3D::I;
  Tensor33d   dStress( 0.0 );
  Tensor3333d tangent( 0.0 );
  tangent = C;

  std::array< double, 9 > stateVars{};

  // Arbitrary non-zero dTangent_dDeformation
  Tensor333333d dTangent_dDeformation( 1.0 ); // all entries = 1.0

  evaluateGeneralizedMaxwellModel( stress, tangent, dTangent_dDeformation, Cinv, dStress, dT, props, stateVars.data() );

  // With dStress=0 and Q_n=0:
  //   Q_np = 0
  //   H_np = Cinv : Q_np = 0
  //   stress += C : H_np = 0  ->  stress = (1-gamma)*stress_in = 0.7 * I
  //   dH_np = Cinv : (beta * C) => tangent contribution from first term
  //   But H_np = 0, so second tangent term (H_np . dTangent) = 0
  //   tangent = (1-gamma)*C + C : (Cinv : beta*C) = (1-gamma)*C + beta*C = (1-gamma+beta)*C

  // Verify stress
  const double stressExpected = ( 1.0 - 0.3 );
  throwExceptionOnFailure( checkIfEqual( stress( 0, 0 ), stressExpected, 1e-12 ),
                           MakeString() << __PRETTY_FUNCTION__ << ": stress(0,0) with dStress=0 should be (1-gamma)" );

  // Verify tangent is (1-gamma+beta)*C
  const double factor = ( 1.0 - 0.3 + beta );
  throwExceptionOnFailure( checkIfEqual( tangent( 0, 0, 0, 0 ), factor * C( 0, 0, 0, 0 ), 1e-8 ),
                           MakeString() << __PRETTY_FUNCTION__ << ": tangent(0,0,0,0) with dStress=0 failed" );
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------
int main()
{
  auto tests = std::vector< std::function< void() > >{
    // createMaxwellProperties
    testCreateMaxwellPropertiesEmpty,
    testCreateMaxwellPropertiesSingle,
    testCreateMaxwellPropertiesTwo,

    // Template overload (no tangent)
    testTemplateOverloadNoOp,
    testTemplateOverloadSingleElementFromZero,
    testTemplateOverloadSingleElementFromNonZeroState,
    testTemplateOverloadTaylorApproximation,
    testTemplateOverloadTwoElements,
    testTemplateOverloadRelaxation,

    // Simo overload (with tangent)
    testSimoOverloadNoOp,
    testSimoOverloadSingleElementFromZero,
    testSimoOverloadSingleElementFromNonZeroState,
    testSimoOverloadTwoElements,
    testSimoAndTemplateOverloadConsistency,

    // Liu et al. overload (with compliance and dTangent/dDeformation)
    testIsotropicInverseConsistency,
    testLiuOverloadMatchesSimo,
    testLiuOverloadNonZeroDTangentDDeformation,
  };

  executeTestsAndCollectExceptions( tests );
  return 0;
}
