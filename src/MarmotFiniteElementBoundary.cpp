#include "Marmot/MarmotFiniteElement.h"
#include <iostream>

using namespace Eigen;
namespace Marmot {
  namespace FiniteElement {

    BoundaryElement::BoundaryElement( ElementShapes   parentShape,
                                      int             parentFaceNumber,
                                      int             nDim,
                                      const VectorXd& parentCoordinates )
      : nDim( nDim )
    {
      // get boundary shape dependent on parental shape
      switch ( parentShape ) {
        /*
         * Extend for future elements here ... (scroll down) ..
         *
         * */

      case Quad4: {
        boundaryShape                = Bar2;
        nNodes                       = Spatial1D::Bar2::nNodes;
        boundaryIndicesInParentNodes = Spatial2D::Quad4::getBoundaryElementIndices( parentFaceNumber );
        break;
      }

      case Quad8: {
        boundaryShape                = Bar3;
        nNodes                       = Spatial1D::Bar3::nNodes;
        boundaryIndicesInParentNodes = Spatial2D::Quad8::getBoundaryElementIndices( parentFaceNumber );
        break;
      }
      case Hexa8: {
        boundaryShape                = Quad4;
        nNodes                       = Spatial2D::Quad4::nNodes;
        boundaryIndicesInParentNodes = Spatial3D::Hexa8::getBoundaryElementIndices( parentFaceNumber );
        break;
      }
      case Hexa20: {
        boundaryShape                = Quad8;
        nNodes                       = Spatial2D::Quad8::nNodes;
        boundaryIndicesInParentNodes = Spatial3D::Hexa20::getBoundaryElementIndices( parentFaceNumber );
        break;
      }

      default: throw std::invalid_argument( "Boundary Element currently not implemented" );
      }

      // get the 'condensed' boundary element coordinates
      nParentCoordinates                 = parentCoordinates.size();
      boundaryIndicesInParentCoordinates = expandNodeIndicesToCoordinateIndices( boundaryIndicesInParentNodes, nDim );
      coordinates                        = condenseParentToBoundaryVector( parentCoordinates );

      // get the proper gausspoints for the boundary element

      // compute for each quadraturePoint the
      // - dNdXi
      // - Jacobian pp(x,xi)
      // - areaVector
      for ( const auto& quadraturePointInfo : FiniteElement::Quadrature::
              getGaussPointInfo( boundaryShape, FiniteElement::Quadrature::IntegrationTypes::FullIntegration ) ) {

        BoundaryElementQuadraturePoint qp;
        qp.xi     = quadraturePointInfo.xi;
        qp.weight = quadraturePointInfo.weight;

        // N
        // dNdXi
        switch ( boundaryShape ) {
        /*
         * ... and extend for future elements here!
         *
         * */
        case Bar2: {
          qp.N     = Spatial1D::Bar2::N( qp.xi( 0 ) );
          qp.dNdXi = Spatial1D::Bar2::dNdXi( qp.xi( 0 ) );
          break;
        }
        case Bar3: {
          qp.N     = Spatial1D::Bar3::N( qp.xi( 0 ) );
          qp.dNdXi = Spatial1D::Bar3::dNdXi( qp.xi( 0 ) );
          break;
        }
        case Quad4: {
          qp.N     = Spatial2D::Quad4::N( qp.xi );
          qp.dNdXi = Spatial2D::Quad4::dNdXi( qp.xi );
          break;
        }
        case Quad8: {
          qp.N     = Spatial2D::Quad8::N( qp.xi );
          qp.dNdXi = Spatial2D::Quad8::dNdXi( Vector2d( qp.xi ) );
          break;
        }

        default: break; // exception handling already in first switch
        }

        // Jacobian
        qp.J = Jacobian( qp.dNdXi, coordinates );

        // areaVector and integration area
        if ( nDim == 2 ) {
          // 90deg rotation
          Vector2d n;
          n << qp.J( 1 ), -qp.J( 0 );
          qp.areaVector = n;
        }
        else {
          // cross product
          typedef Eigen::Ref< const Vector3d > vector3_cr;
          Vector3d                             n = vector3_cr( qp.J.col( 0 ) ).cross( vector3_cr( qp.J.col( 1 ) ) );
          qp.areaVector                          = n;
        }

        quadraturePoints.push_back( std::move( qp ) );
      }
    }

    VectorXd BoundaryElement::computeNormalLoadVector()
    {
      /* compute the load vector for a constant distributed load (e.g. pressure)
       * Attention: result =  boundary-element-sized!
       *  -> use expandBoundaryToParentVector to obtain the parent-element-sized load vector
       * */

      VectorXd Pk = VectorXd::Zero( coordinates.size() );

      for ( const auto& qp : quadraturePoints )
        Pk += qp.weight * qp.areaVector.transpose() * NB( qp.N, nDim );

      return Pk;
    }

    MatrixXd BoundaryElement::computeNormalLoadVectorStiffness()
    {
      MatrixXd K = MatrixXd::Zero( coordinates.size(), coordinates.size() );

      if ( nDim == 2 ) {
        // Neuner, November 2018
        Matrix2d R;
        // clang-format off
                R << 0, 1, 
                    -1, 0;
        // clang-format on

        for ( const auto& qp : quadraturePoints )
          for ( int I = 0; I < nNodes; I++ )
            for ( int J = 0; J < nNodes; J++ )
              K.block< 2, 2 >( I * 2, J * 2 ) += qp.N( I ) * qp.dNdXi( J ) * R * qp.weight;
      }

      else if ( nDim == 3 ) {
        // Belytschko et al. 2014, pp.364
        Matrix3d HXi0, HXi1;
        for ( const auto& qp : quadraturePoints ) {
          const MatrixXd& J = qp.J;

          // clang-format off
                    HXi0 << 0,       J(2,0), -J(1,0), 
                           -J(2,0),  0,       J(0,0), 
                            J(1,0), -J(0,0),  0;

                    HXi1 << 0,       J(2,1), -J(1,1), 
                           -J(2,1),  0,       J(0,1), 
                            J(1,1), -J(0,1),  0;
          // clang-format on

          for ( int I = 0; I < nNodes; I++ )
            for ( int J = 0; J < nNodes; J++ )
              K.block< 3, 3 >( I * 3, J * 3 ) += qp.N( I ) * ( qp.dNdXi( 0, J ) * HXi1 - qp.dNdXi( 1, J ) * HXi0 ) *
                                                 qp.weight;
        }
      }

      return K;
    }

    VectorXd BoundaryElement::condenseParentToBoundaryVector( const VectorXd& parentVector )
    {
      /*  condense any parent vector to the corresponding boundary child vector (e.g.
       * coordinates) dependent on the underlying indices mapping
       * */
      VectorXd boundaryVector( nNodes * nDim );

      for ( int i = 0; i < boundaryIndicesInParentCoordinates.size(); i++ )
        boundaryVector( i ) = parentVector( boundaryIndicesInParentCoordinates( i ) );

      return boundaryVector;
    }

    void BoundaryElement::assembleIntoParentVector( const Eigen::VectorXd&        boundaryVector,
                                                    Eigen::Ref< Eigen::VectorXd > parentVector )
    {
      /* assemble any boundary vector (e.g. pressure load) into the corresponding parental vector
       * */

      for ( int i = 0; i < boundaryVector.size(); i++ )
        parentVector( boundaryIndicesInParentCoordinates( i ) ) += boundaryVector( i );
    }

    void BoundaryElement::assembleIntoParentStiffness( const Eigen::MatrixXd&        KBoundary,
                                                       Eigen::Ref< Eigen::MatrixXd > KParent )
    {
      // mind the negative sign for a proper assembly of the residual(!) stiffness !
      for ( int i = 0; i < KBoundary.cols(); i++ )
        for ( int j = 0; j < KBoundary.cols(); j++ )
          KParent( boundaryIndicesInParentCoordinates( i ), boundaryIndicesInParentCoordinates( j ) ) -= KBoundary( i,
                                                                                                                    j );
    }

  } // namespace FiniteElement
} // namespace Marmot
