#pragma once
// Backend adapter: Marmot's own vector mode, Marmot::Math::Differentiation::Automatic::VectorMode::Dual (static,
// expression templates with eagerly evaluated values and partial derivatives). First derivatives only: it neither
// nests nor has a scalar (one direction per level) or Taylor mode.
#include "ADBackendNestedFad.h"
#include "Marmot/MarmotAutomaticDifferentiationVectorMode.h"
#include <cstddef>

namespace ADBench {

  namespace MarmotDetail {
    struct Access {
      template < typename S >
      static auto& val( S& x )
      {
        return x.val();
      }
      template < typename S >
      static auto val( const S& x )
      {
        return x.val();
      }
      template < typename S >
      static auto& dx( S& x, std::size_t i )
      {
        return x.dx( i );
      }
      template < typename S >
      static auto dx( const S& x, std::size_t i )
      {
        return x.dx( i );
      }
      template < typename S >
      static S constant( double v )
      {
        return S( v );
      }
    };
  } // namespace MarmotDetail

  struct MarmotBackend : NestedFadOps< MarmotDetail::Access > {
    static constexpr const char* name           = "marmot";
    static constexpr bool        hasReverseMode = false;
    static constexpr bool        hasTaylor      = false;
    static constexpr bool        hasVectorMode  = true;
    static constexpr int         maxScalarOrder = 0;
    static constexpr int         maxVectorOrder = 1;

    template < std::size_t N, std::size_t depth >
      requires( depth == 1 )
    using VecDual = Marmot::Math::Differentiation::Automatic::VectorMode::Dual< N >;

    template < typename S >
    static S constant( double v )
    {
      return S( v );
    }
  };

} // namespace ADBench
