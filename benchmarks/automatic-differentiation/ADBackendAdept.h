#pragma once
// Backend adapter: Adept-2 (https://github.com/rjhogan/Adept-2, Apache-2.0), operator-overloading AD with a tape.
// One evaluation is recorded on a (thread-local) adept::Stack; the full Jacobian then follows from forward or reverse
// sweeps over the recording (Adept picks the cheaper direction). First order only; no higher derivatives and no
// Taylor mode, so it implements the two first-order kernels itself instead of the number-type vocabulary.
#include "ADKernels.h"
#include <Eigen/Core>
#include <adept.h>

namespace Eigen {
  template <>
  struct NumTraits< adept::adouble > : NumTraits< double > {
    using Real       = adept::adouble;
    using NonInteger = adept::adouble;
    using Literal    = double;
    using Nested     = adept::adouble;
    enum {
      IsComplex             = 0,
      IsInteger             = 0,
      IsSigned              = 1,
      RequireInitialization = 1,
      ReadCost              = 1,
      AddCost               = 3,
      MulCost               = 3
    };
  };
  template < typename BinOp >
  struct ScalarBinaryOpTraits< adept::adouble, double, BinOp > {
    using ReturnType = adept::adouble;
  };
  template < typename BinOp >
  struct ScalarBinaryOpTraits< double, adept::adouble, BinOp > {
    using ReturnType = adept::adouble;
  };
  namespace internal {
    // adouble's converting constructor accepts any type, so Eigen's scalar-promotion fallback would treat a whole
    // matrix as a scalar and make `A * B` ambiguous (same fix as for Boost's fvar). Only arithmetic types promote.
    template < typename T, bool >
    struct PromoteToAdoubleIf {};
    template < typename T >
    struct PromoteToAdoubleIf< T, true > {
      using type = adept::adouble;
    };
    template < typename T >
    struct promote_scalar_arg< adept::adouble, T, false > : PromoteToAdoubleIf< T, std::is_arithmetic_v< T > > {};
  } // namespace internal
} // namespace Eigen

namespace ADBench {

  struct AdeptBackend {
    static constexpr const char* name           = "adept";
    static constexpr bool        hasReverseMode = true;
    static constexpr bool        hasTaylor      = false;
    static constexpr int         maxScalarOrder = 0;
    static constexpr bool        hasVectorMode  = false;
    static constexpr int         maxVectorOrder = 0;
    // Recording through Fastor's einsum yields correct values but wrong Jacobians (even a two-operand
    // einsum< iI, IJ, to_ij >), so the finite-strain kernel, which contracts with einsum, is not run. The Eigen-based
    // hypoelastic kernel is exact.
    static constexpr bool reverseFastorEinsumIsCorrect = false;

    template < typename S >
    static S constant( double v )
    {
      return S( v );
    }

    template < typename S >
    static double value( const S& x )
    {
      return adept::value( x );
    }

    // an adept::Stack activates itself on construction; one per thread, as a multi-threaded FE code would need
    static adept::Stack& stack()
    {
      static thread_local adept::Stack s;
      return s;
    }

    static std::pair< Tensor33d, Tensor3333d > finiteStrainTangent( const Tensor33d& F )
    {
      using S              = adept::adouble;
      adept::Stack&     st = stack();
      Tensor< S, 3, 3 > F_;
      for ( std::size_t a = 0; a < 9; a++ )
        F_.data()[a] = F.data()[a];
      st.new_recording();
      Tensor< S, 3, 3 > tau = kirchhoffStressPenceGou< AdeptBackend, S >( F_ );
      st.independent( F_.data(), 9 );
      st.dependent( tau.data(), 9 );
      double jac[81]; // column major: jac[ j + 9 a ] = d tau_j / d F_a
      st.jacobian( jac );

      Tensor33d   tau_;
      Tensor3333d dTau_dF;
      for ( std::size_t j = 0; j < 9; j++ ) {
        tau_.data()[j] = adept::value( tau.data()[j] );
        for ( std::size_t a = 0; a < 9; a++ )
          dTau_dF.data()[j * 9 + a] = jac[j + 9 * a];
      }
      return { tau_, dTau_dF };
    }

    static std::pair< Vector6d, Matrix6d > hypoElasticJ2( const Vector6d& stressOld, const Vector6d& dStrain )
    {
      using S                     = adept::adouble;
      adept::Stack&            st = stack();
      Eigen::Matrix< S, 6, 1 > dE;
      for ( int j = 0; j < 6; j++ )
        dE( j ) = dStrain( j );
      st.new_recording();
      Eigen::Matrix< S, 6, 1 > stress = hypoElasticJ2Stress< AdeptBackend, S >( stressOld, dE );
      st.independent( dE.data(), 6 );
      st.dependent( stress.data(), 6 );
      Matrix6d C; // column major, as Adept's default layout
      st.jacobian( C.data() );

      Vector6d S_;
      for ( int i = 0; i < 6; i++ )
        S_( i ) = adept::value( stress( i ) );
      return { S_, C };
    }
  };

} // namespace ADBench
