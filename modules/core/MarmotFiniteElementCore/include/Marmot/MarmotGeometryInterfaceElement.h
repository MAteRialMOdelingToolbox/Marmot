#pragma once

#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotTypedefs.h"

#include <Eigen/Dense>
#include <cmath>
#include <map>
#include <stdexcept>
#include <string>

template < int nDim, int nNodes >
class MarmotGeometryInterfaceElement {
public:
  static_assert( nDim == 2 || nDim == 3,
                 "MarmotGeometryInterfaceElement supports interfaces embedded in 2D or 3D only." );

  static_assert( nNodes % 2 == 0, "Interface element must have bottom and top sides, so nNodes must be even." );

  static constexpr int nXi             = nDim - 1;
  static constexpr int nInterfaceNodes = nNodes / 2;

  static constexpr int nDofPerNodeU = nDim;
  static constexpr int nDofPerSide  = nInterfaceNodes * nDim;
  static constexpr int nDofElement  = nNodes * nDim;

  using XiSized    = Eigen::Matrix< double, nXi, 1 >;
  using NSized     = Eigen::Matrix< double, 1, nInterfaceNodes >;
  using dNdXiSized = Eigen::Matrix< double, nXi, nInterfaceNodes >;

  using CoordinateVector     = Eigen::Matrix< double, nDim * nNodes, 1 >;
  using SideCoordinateVector = Eigen::Matrix< double, nDim * nInterfaceNodes, 1 >;

  using SurfaceJacobianSized = Eigen::Matrix< double, nDim, nXi >;
  using MetricSized          = Eigen::Matrix< double, nXi, nXi >;

  using GradSized = Eigen::Matrix< double, nDim, nInterfaceNodes >;

  using VectorDim = Eigen::Matrix< double, nDim, 1 >;
  using TensorDim = Eigen::Matrix< double, nDim, nDim >;

  using NMatrixSized     = Eigen::Matrix< double, nDim, nDofPerSide >;
  using NJumpMatrixSized = Eigen::Matrix< double, nDim, nDofElement >;

  using BSurfaceSized    = Eigen::Matrix< double, nDim * nDim, nDofPerSide >;
  using BAvgSurfaceSized = Eigen::Matrix< double, nDim * nDim, nDofElement >;

  Eigen::Map< const CoordinateVector > coordinates;

  const Marmot::FiniteElement::ElementShapes shape;

  MarmotGeometryInterfaceElement() : coordinates( nullptr ), shape( getInterfaceElementShape() ) {}

  static constexpr Marmot::FiniteElement::ElementShapes getInterfaceElementShape()
  {
    using namespace Marmot::FiniteElement;

    if constexpr ( nDim == 2 && nInterfaceNodes == 2 ) {
      return Bar2;
    }
    else if constexpr ( nDim == 3 && nInterfaceNodes == 4 ) {
      return Quad4;
    }
    else {
      static_assert( nDim == -1, "Unsupported MarmotGeometryInterfaceElement interpolation in the current build." );
    }
  }

  std::string getElementShape() const
  {
    using namespace Marmot::FiniteElement;

    static std::map< ElementShapes, std::string > shapes = {
      { Bar2, "iline2" },
      { Quad4, "iquad4" },
    };

    return shapes[this->shape];
  }

  void assignNodeCoordinates( const double* coords )
  {
    new ( &coordinates ) Eigen::Map< const CoordinateVector >( coords );
  }

  NSized     N( const XiSized& xi ) const;
  dNdXiSized dNdXi( const XiSized& xi ) const;

  SideCoordinateVector getSideCoordinates( const int side = 0 ) const
  {
    if ( side != 0 && side != 1 )
      throw std::invalid_argument( "Interface side must be 0 or 1." );

    SideCoordinateVector xSide;

    const int nodeOffset = side * nInterfaceNodes;

    for ( int A = 0; A < nInterfaceNodes; ++A ) {
      for ( int i = 0; i < nDim; ++i ) {
        xSide( A * nDim + i ) = coordinates( ( nodeOffset + A ) * nDim + i );
      }
    }

    return xSide;
  }

  SurfaceJacobianSized surfaceJacobian( const dNdXiSized& dN, const int side = 0 ) const
  {
    const SideCoordinateVector xSide = getSideCoordinates( side );

    SurfaceJacobianSized J;
    J.setZero();

    for ( int A = 0; A < nInterfaceNodes; ++A ) {
      for ( int i = 0; i < nDim; ++i ) {
        const double XAi = xSide( A * nDim + i );

        for ( int alpha = 0; alpha < nXi; ++alpha ) {
          J( i, alpha ) += XAi * dN( alpha, A );
        }
      }
    }

    return J;
  }

  MetricSized metric( const SurfaceJacobianSized& J ) const { return J.transpose() * J; }

  double sqrtDetMetric( const SurfaceJacobianSized& J ) const
  {
    const MetricSized G    = metric( J );
    const double      detG = G.determinant();

    if ( detG <= 0.0 )
      throw std::runtime_error( "Degenerate interface element: det(G) <= 0." );

    return std::sqrt( detG );
  }

  GradSized surfaceGradient( const dNdXiSized& dN, const SurfaceJacobianSized& J ) const
  {
    const MetricSized GInv = metric( J ).inverse();

    /*
     * Surface gradient of scalar shape functions:
     *
     *   grad_s N_A = J * G^{-1} * dN_A/dxi
     *
     * Shape:
     *   J       : nDim x nXi
     *   GInv    : nXi  x nXi
     *   dN      : nXi  x nInterfaceNodes
     *   gradN   : nDim x nInterfaceNodes
     */
    return J * GInv * dN;
  }

  VectorDim normal( const SurfaceJacobianSized& J ) const
  {
    VectorDim n;
    n.setZero();

    if constexpr ( nDim == 2 ) {
      const VectorDim t = J.col( 0 );

      n( 0 ) = -t( 1 );
      n( 1 ) = t( 0 );
    }
    else if constexpr ( nDim == 3 ) {
      n = J.col( 0 ).cross( J.col( 1 ) );
    }

    const double normN = n.norm();

    if ( normN < 1e-16 )
      throw std::runtime_error( "Degenerate interface element: normal vector is zero." );

    return n / normN;
  }

  TensorDim normalProjector( const VectorDim& n ) const { return n * n.transpose(); }

  TensorDim tangentProjector( const VectorDim& n ) const { return TensorDim::Identity() - normalProjector( n ); }

  NMatrixSized NMatrix( const NSized& N_ ) const
  {
    NMatrixSized Nmat;
    Nmat.setZero();

    for ( int A = 0; A < nInterfaceNodes; ++A ) {
      for ( int i = 0; i < nDim; ++i ) {
        Nmat( i, A * nDim + i ) = N_( 0, A );
      }
    }

    return Nmat;
  }

  NJumpMatrixSized NJumpMatrix( const NSized& N_ ) const
  {
    NJumpMatrixSized NJump;
    NJump.setZero();

    for ( int A = 0; A < nInterfaceNodes; ++A ) {
      for ( int i = 0; i < nDim; ++i ) {
        const int colBottom = A * nDim + i;
        const int colTop    = ( nInterfaceNodes + A ) * nDim + i;

        NJump( i, colBottom ) = -N_( 0, A );
        NJump( i, colTop )    = N_( 0, A );
      }
    }

    return NJump;
  }

  BSurfaceSized BSurfaceMatrix( const GradSized& gradN, const TensorDim& T ) const
  {
    /*
     * Python-compatible surface displacement-gradient operator.
     *
     * This intentionally does NOT build the classical symmetric/Cauchy strain
     * matrix. The Python EdelweissFE interface element uses the full projected
     * displacement gradient:
     *
     *   surface_grad[a, i, k, q]
     *
     * where:
     *
     *   a = interface node index
     *   i = displacement component
     *   k = gradient direction
     *   q = quadrature point
     *
     * For dim == 2, Python does:
     *
     *   surface_grad[a, i, k] = gradN[j, A] * T[j, k]
     *
     * For dim == 3, Python does:
     *
     *   surface_grad[a, i, k] = T[i, m] * gradN[j, A] * T[j, k]
     *
     * Then calculate_B_surface_grad writes it into the B matrix as:
     *
     *   row = i * nDim + k
     *   col = A * nDim + i
     *
     * Therefore each displacement component i only writes into its own
     * component column, not into all vector columns.
     */
    BSurfaceSized B;
    B.setZero();

    for ( int A = 0; A < nInterfaceNodes; ++A ) {
      for ( int i = 0; i < nDim; ++i ) {
        for ( int k = 0; k < nDim; ++k ) {

          double value = 0.0;

          if constexpr ( nDim == 2 ) {
            for ( int j = 0; j < nDim; ++j ) {
              value += gradN( j, A ) * T( j, k );
            }
          }
          else if constexpr ( nDim == 3 ) {
            for ( int m = 0; m < nDim; ++m ) {
              for ( int j = 0; j < nDim; ++j ) {
                value += T( i, m ) * gradN( j, A ) * T( j, k );
              }
            }
          }

          const int row = i * nDim + k;
          const int col = A * nDim + i;

          B( row, col ) = value;
        }
      }
    }

    return B;
  }

  BAvgSurfaceSized BAverageSurfaceMatrix( const BSurfaceSized& Bside ) const
  {
    BAvgSurfaceSized BAvg;
    BAvg.setZero();

    for ( int row = 0; row < nDim * nDim; ++row ) {
      for ( int col = 0; col < nDofPerSide; ++col ) {
        BAvg( row, col )               = 0.5 * Bside( row, col );
        BAvg( row, nDofPerSide + col ) = 0.5 * Bside( row, col );
      }
    }

    return BAvg;
  }

  BSurfaceSized BSurfaceMatrixFullyProjected( const GradSized& gradN, const TensorDim& T ) const
  {
    /*
     * Backward-compatible alias.
     *
     * The Python implementation already uses the fully projected
     * displacement-gradient operator. Keeping this function avoids changing
     * callers that still pass fullyProjectedB=true.
     */
    return BSurfaceMatrix( gradN, T );
  }

  struct QuadratureGeometry {
    NSized     N;
    dNdXiSized dNdXi;

    SurfaceJacobianSized J;
    MetricSized          G;
    double               sqrtDetG;

    GradSized gradN;

    VectorDim n;
    TensorDim normalProjection;
    TensorDim tangentProjection;

    NMatrixSized     NmatSide;
    NJumpMatrixSized NmatJump;

    BSurfaceSized    BmatSide;
    BAvgSurfaceSized BmatAverage;
  };

  QuadratureGeometry evaluateAt( const XiSized& xi,
                                 const int      sideForGeometry = 0,
                                 const bool     fullyProjectedB = false ) const
  {
    QuadratureGeometry q;

    q.N     = N( xi );
    q.dNdXi = dNdXi( xi );

    q.J        = surfaceJacobian( q.dNdXi, sideForGeometry );
    q.G        = metric( q.J );
    q.sqrtDetG = sqrtDetMetric( q.J );

    q.gradN = surfaceGradient( q.dNdXi, q.J );

    q.n = normal( q.J );

    q.normalProjection  = normalProjector( q.n );
    q.tangentProjection = tangentProjector( q.n );

    q.NmatSide = NMatrix( q.N );
    q.NmatJump = NJumpMatrix( q.N );

    /*
     * Both branches now intentionally produce the Python-compatible projected
     * displacement-gradient operator.
     */
    if ( fullyProjectedB )
      q.BmatSide = BSurfaceMatrixFullyProjected( q.gradN, q.tangentProjection );
    else
      q.BmatSide = BSurfaceMatrix( q.gradN, q.tangentProjection );

    q.BmatAverage = BAverageSurfaceMatrix( q.BmatSide );

    return q;
  }
};