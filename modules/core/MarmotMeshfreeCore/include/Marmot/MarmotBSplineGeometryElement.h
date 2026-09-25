#pragma once
#include "Eigen/Core"
#include "Marmot/MarmotBSpline.h"
#include "Marmot/MarmotFiniteElement.h"
#include <iostream>
#include <map>

namespace Marmot::FiniteElement {

  template < int nDim, int order >
  class MarmotBSplineGeometryElement {

    static int constexpr pow( int base, int exponent ) { return exponent == 0 ? 1 : base * pow( base, exponent - 1 ); }

  public:
    constexpr static int nKnotsPerDir = 2 * order + 2;
    using KnotVector                  = Eigen::Matrix< double, nKnotsPerDir, 1 >;
    using KnotVectors                 = Eigen::Matrix< double, nKnotsPerDir, nDim >;

    constexpr static int nNodes = pow( order + 1, nDim );
    using NSized                = Eigen::Matrix< double, 1, nNodes >;
    using CoordinateVector      = Eigen::Matrix< double, nDim * nNodes, 1 >;
    using JacobianSized         = Eigen::Matrix< double, nDim, nDim >;

    using dNdXSized = Eigen::Matrix< double, nDim, nNodes >;
    using XiSized   = Eigen::Matrix< double, nDim, 1 >;

    const Eigen::Map< const CoordinateVector > _mapCoordinates;
    const Marmot::FiniteElement::ElementShapes _shape;

    KnotVectors _knotVectors;
    KnotVectors _knotVectorsNormalized;

    MarmotBSplineGeometryElement( const double* nodeCoordinates_,
                                  int           sizeNodeCoordinates,
                                  const double* knotVectors_,
                                  int           sizeKnotVectors )
      : _mapCoordinates( Eigen::Map< const CoordinateVector >( nodeCoordinates_ ) ),
        _shape( Marmot::FiniteElement::getElementShapeByMetric( nDim, nNodes ) )
    {
      _knotVectors = Eigen::Map< const KnotVectors >( knotVectors_ );

      for ( int d = 0; d < nDim; d++ ) {

        /* knotVectors_ */
        _knotVectorsNormalized = _knotVectors;
        /* _knotVectorsNormalized.col( */
        /*   d ) = ( _knotVectors.array().col( d ) - _knotVectors.col( d )( order ) ) * */
        /*         ( 2 / ( _knotVectors.col( d )( nKnotsPerDir - order - 1 ) - _knotVectors.col( d )( order ) ) ); */

        /* _knotVectorsNormalized.array().col( d ) -= 1; */
      }
    };

    std::string getElementShape() const
    {
      using namespace Marmot::FiniteElement;
      static std::map< ElementShapes, std::string > shapes = {
        { Bar2, "bar2" },
        { Quad4, "quad4" },
        { Quad8, "quad8" },
        { Quad9, "quad9" },
        { Hexa8, "hexa8" },
        { Hexa20, "hexa20" },
        { Hexa27, "hexa27" },
      };

      return shapes[this->_shape];
    }

    NSized N( const XiSized& xi ) const
    {
      NSized        N_;
      constexpr int nN = order + 1;

      if constexpr ( nDim == 1 )
        for ( int p = 0; p < nN; p++ )
          N_( p ) = B< order >( xi( 0 ), this->_knotVectorsNormalized.col( 0 ).data(), p );

      else if constexpr ( nDim == 2 )
        for ( int q = 0; q < nN; q++ )
          for ( int p = 0; p < nN; p++ )
            N_( p + q * ( nN ) ) = B< order >( xi( 0 ), this->_knotVectorsNormalized.col( 0 ).data(), p ) *
                                   B< order >( xi( 1 ), this->_knotVectorsNormalized.col( 1 ).data(), q );

      else if constexpr ( nDim == 3 )
        for ( int r = 0; r < nN; r++ )
          for ( int q = 0; q < nN; q++ )
            for ( int p = 0; p < nN; p++ )
              N_( p + q * nN + r * nN * nN ) = B< order >( xi( 0 ), this->_knotVectorsNormalized.col( 0 ).data(), p ) *
                                               B< order >( xi( 1 ), this->_knotVectorsNormalized.col( 1 ).data(), q ) *
                                               B< order >( xi( 2 ), this->_knotVectorsNormalized.col( 2 ).data(), r );

      return N_;
    }

    dNdXSized dNdXi( const XiSized& xi ) const
    {
      dNdXSized     dN_dXi_;
      constexpr int nN = order + 1;

      if constexpr ( nDim == 1 )
        for ( int p = 0; p < nN; p++ )
          dN_dXi_( p ) = dB_dU< order >( xi( 0 ), this->_knotVectorsNormalized[0].data(), p );

      else if constexpr ( nDim == 2 )
        for ( int q = 0; q < nN; q++ )
          for ( int p = 0; p < nN; p++ ) {
            dN_dXi_( 0, p + q * nN ) = dB_dU< order >( xi( 0 ), this->_knotVectorsNormalized.col( 0 ).data(), p ) *
                                       B< order >( xi( 1 ), this->_knotVectorsNormalized.col( 1 ).data(), q );
            dN_dXi_( 1, p + q * nN ) = B< order >( xi( 0 ), this->_knotVectorsNormalized.col( 0 ).data(), p ) *
                                       dB_dU< order >( xi( 1 ), this->_knotVectorsNormalized.col( 1 ).data(), q );
          }

      else if constexpr ( nDim == 3 )
        for ( int r = 0; r < nN; r++ )
          for ( int q = 0; q < nN; q++ )
            for ( int p = 0; p < nN; p++ ) {
              dN_dXi_( 0,
                       p + q * nN +
                         r * nN * nN ) = dB_dU< order >( xi( 0 ), this->_knotVectorsNormalized.col( 0 ).data(), p ) *
                                         B< order >( xi( 1 ), this->_knotVectorsNormalized.col( 1 ).data(), q ) *
                                         B< order >( xi( 2 ), this->_knotVectorsNormalized.col( 2 ).data(), r );
              dN_dXi_( 1,
                       p + q * nN +
                         r * nN * nN ) = B< order >( xi( 0 ), this->_knotVectorsNormalized.col( 0 ).data(), p ) *
                                         dB_dU< order >( xi( 1 ), this->_knotVectorsNormalized.col( 1 ).data(), q ) *
                                         B< order >( xi( 2 ), this->_knotVectorsNormalized.col( 2 ).data(), r );
              dN_dXi_( 2,
                       p + q * nN +
                         r * nN * nN ) = B< order >( xi( 0 ), this->_knotVectorsNormalized.col( 0 ).data(), p ) *
                                         B< order >( xi( 1 ), this->_knotVectorsNormalized.col( 1 ).data(), q ) *
                                         dB_dU< order >( xi( 2 ), this->_knotVectorsNormalized.col( 2 ).data(), r );
            }

      return dN_dXi_;
    }

    JacobianSized Jacobian( const dNdXSized& dNdXi ) const
    {
      return Marmot::FiniteElement::Jacobian< nDim, nNodes >( dNdXi, _mapCoordinates );
    }

    dNdXSized dNdX( const dNdXSized& dNdXi, const JacobianSized& JacobianInverse ) const
    {
      return ( dNdXi.transpose() * JacobianInverse ).transpose();
    }
  };

} // namespace Marmot::FiniteElement
