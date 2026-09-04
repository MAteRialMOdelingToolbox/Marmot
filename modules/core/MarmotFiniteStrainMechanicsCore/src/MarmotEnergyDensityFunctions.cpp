#include "Marmot/MarmotEnergyDensityFunctions.h"
#include "Marmot/MarmotAutomaticDifferentiationForFastor.h"
namespace Marmot::ContinuumMechanics::EnergyDensityFunctions::ThirdOrderDerived {
  using namespace Fastor;
  using namespace Marmot::TensorUtility::FastorTensors::StandardTensors;
  std::tuple< double,
              TensorUtility::FastorTensors::StandardTensors::Tensor33d,
              TensorUtility::FastorTensors::StandardTensors::Tensor3333d,
              TensorUtility::FastorTensors::StandardTensors::Tensor333333d >
  standardNeoHooke( const TensorUtility::FastorTensors::StandardTensors::Tensor33d& C,
                    const double&                                                   K,
                    const double&                                                   G )
  {

    const double lambda = K - 2.0 / 3.0 * G;

    /*
     * we use the potential in terms of C
     * Psi = G/2 ( tr(C) - 3 - 2 ln(J) ) + lambda/2 ( 0.5 ( J^2 - 1 ) - ln(J) )
     * where J = sqrt( det(C) ) and ln(J) = 0.5 ln( det(C) )
     *
     * we can the rewrite as
     * Psi = G/2 ( tr(C) - 3 - ln(det(C)) ) + lambda/2 ( 0.5 ( det(C) - 1 ) - 0.5 ln(det(C)) )
     *
     */

    const double trC    = trace( C );
    const double detC   = determinant( C );
    const double lnDetC = log( detC );
    // energy density
    const double psi = G / 2 * ( trC - 3.0 - lnDetC ) + lambda / 4 * ( detC - 1.0 - lnDetC );

    // first derivative quantities
    const Tensor33d& dTrC_dC    = Spatial3D::I;
    const Tensor33d  invCt      = transpose( inverse( C ) );
    const Tensor33d  dDetC_dC   = detC * invCt;
    const Tensor33d& dLnDetC_dC = invCt;

    // first derivative with respect to C
    const Tensor33d dPsi_dC = G / 2 * ( dTrC_dC - dLnDetC_dC ) + lambda / 4 * ( dDetC_dC - dLnDetC_dC );

    // second derivative quantities
    using namespace TensorUtility::FastorTensors::Indices;
    using lj                        = Index< l_, j_ >;
    using li                        = Index< l_, i_ >;
    const Tensor3333d  dInvCt_dC    = -0.5 * ( einsum< ik, lj, to_ijkl >( invCt, invCt ) +
                                           einsum< jk, li, to_ijkl >( invCt, invCt ) );
    const Tensor3333d  d2DetC_dC2   = detC * dInvCt_dC + einsum< ij, kl, to_ijkl >( invCt, dDetC_dC );
    const Tensor3333d& d2LnDetC_dC2 = dInvCt_dC;

    // second derivative with respect to C
    const Tensor3333d d2Psi_dC2 = lambda / 4.0 * d2DetC_dC2 - ( G / 2.0 + lambda / 4.0 ) * d2LnDetC_dC2;

    // // third derivative quantities
    // using nk = Index<n_, k_>;
    // using nj = Index<n_, j_>;
    // using ni = Index<n_, i_>;
    // using im = Index<i_, m_>;
    // using jm = Index<j_, m_>;
    // using nl = Index<n_, l_>;
    // using il = Index<i_, l_>;
    // using klmn = Index<k_, l_, m_, n_>;
    // using to_ijklmn = OIndex<i_, j_, k_, l_, m_, n_>;
    // const Tensor333333d d2InvCt_dC2_analytical = 0.25 * (
    //          // derivative of C^-1_ik C^-1_lj w.r.t. C_mn
    //          einsum< im, nk, jl, to_ijklmn >( invCt, invCt, invCt )
    //         + einsum< km, ni, jl, to_ijklmn >( invCt, invCt, invCt )
    //         + einsum< jm, nl, ik, to_ijklmn >( invCt, invCt, invCt )
    //         + einsum< lm, nj, ik, to_ijklmn >( invCt, invCt, invCt )
    //         // derivative of C^-1_jk C^-1_li w.r.t. C_mn
    //         // replace indices accordingly
    //         // iklj -> jkli
    //         + einsum< jm, nk, il, to_ijklmn >( invCt, invCt, invCt )
    //         + einsum< km, nj, il, to_ijklmn >( invCt, invCt, invCt )
    //         + einsum< im, nl, jk, to_ijklmn >( invCt, invCt, invCt )
    //         + einsum< lm, ni, jk, to_ijklmn >( invCt, invCt, invCt )
    //         );
    // const Tensor333333d& d3LnDetC_dC3_analytical = d2InvCt_dC2_analytical;
    // const Tensor333333d d3DetC_dC3_analytical = detC * d2InvCt_dC2_analytical
    //         + einsum< mn, ijkl, to_ijklmn >( dDetC_dC, dInvCt_dC )
    //         + einsum< ij, klmn, to_ijklmn >( invCt,d2DetC_dC2 );
    //         + einsum< kl, ijmn, to_ijklmn >( dDetC_dC, dInvCt_dC );

    // workaround with autodiff for third derivative
    auto detC_func = [=]( const Tensor33t< autodiff::dual3rd >& C_var ) {
      autodiff::dual3rd res = determinant( C_var );
      return res;
    };

    auto lnDetC_func = [=]( const Tensor33t< autodiff::dual3rd >& C_var ) { return log( determinant( C_var ) ); };

    Tensor333333d d3DetC_dC3 = std::get< 3 >(
      Marmot::NumericalAlgorithms::AutomaticDifferentiation::ThirdOrder::d3f_dT3< 3 >( detC_func, C ) );

    Tensor333333d d3LnDetC_dC3 = std::get< 3 >(
      Marmot::NumericalAlgorithms::AutomaticDifferentiation::ThirdOrder::d3f_dT3< 3 >( lnDetC_func, C ) );

    // compare with analytical
    // std::cout << "Max abs diff d3DetC_dC3: " <<  d3DetC_dC3 - d3DetC_dC3_analytical << std::endl;
    // std::cout << "Max abs diff d3LnDetC_dC3: " << d3LnDetC_dC3 - d3LnDetC_dC3_analytical << std::endl;

    // third derivative with respect to C
    const Tensor333333d d3Psi_dC3 = lambda / 4.0 * d3DetC_dC3 - ( G / 2.0 + lambda / 4.0 ) * d3LnDetC_dC3;

    return { psi, dPsi_dC, d2Psi_dC2, d3Psi_dC3 };
  }
} // namespace Marmot::ContinuumMechanics::EnergyDensityFunctions::ThirdOrderDerived
