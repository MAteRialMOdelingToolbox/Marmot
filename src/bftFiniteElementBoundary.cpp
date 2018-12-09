#include "bftFiniteElement.h"
#include <iostream>

using namespace Eigen;
namespace bft {
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
                boundaryShape                = Truss2;
                nNodes                       = Spatial1D::Truss2::nNodes;
                boundaryIndicesInParentNodes = Spatial2D::Quad4::getBoundaryElementIndices( parentFaceNumber );
                break;
            }

            case Quad8: {
                boundaryShape                = Truss3;
                nNodes                       = Spatial1D::Truss3::nNodes;
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
            boundaryIndicesInParentCoordinates = expandNodeIndicesToCoordinateIndices( boundaryIndicesInParentNodes,
                                                                                       nDim );
            coordinates                        = condenseParentToBoundaryVector( parentCoordinates );

            // get the proper gausspoints for the boundary element

            // compute for each gaussPt the
            // - dNdXi
            // - Jacobian pp(x,xi)
            // - normalVector
            // - integrationArea
            for ( const auto& gaussPtInfo :
                  NumIntegration::getGaussPointInfo( boundaryShape,
                                                     NumIntegration::IntegrationTypes::FullIntegration ) ) {

                BoundaryElementGaussPt gpt;
                gpt.xi     = gaussPtInfo.xi;
                gpt.weight = gaussPtInfo.weight;

                // N
                // dNdXi
                switch ( boundaryShape ) {
                /*
                 * ... and extend for future elements here!
                 *
                 * */
                case Truss2: {
                    gpt.N     = Spatial1D::Truss2::N( gpt.xi( 0 ) );
                    gpt.dNdXi = Spatial1D::Truss2::dNdXi( gpt.xi( 0 ) );
                    break;
                }
                case Truss3: {
                    gpt.N     = Spatial1D::Truss3::N( gpt.xi( 0 ) );
                    gpt.dNdXi = Spatial1D::Truss3::dNdXi( gpt.xi( 0 ) );
                    break;
                }
                case Quad4: {
                    gpt.N     = Spatial2D::Quad4::N( gpt.xi );
                    gpt.dNdXi = Spatial2D::Quad4::dNdXi( gpt.xi );
                    break;
                }
                case Quad8: {
                    gpt.N     = Spatial2D::Quad8::N( gpt.xi );
                    gpt.dNdXi = Spatial2D::Quad8::dNdXi( Vector2d( gpt.xi ) );
                    break;
                }

                default: break; // exception handling already in first switch
                }

                // Jacobian
                gpt.J = Jacobian( gpt.dNdXi, coordinates );

                // normalVector and integration area
                if ( nDim == 2 ) {
                    // 90deg rotation
                    Vector2d n;
                    n << gpt.J( 1 ), -gpt.J( 0 );
                    gpt.integrationArea = n.norm();
                    gpt.normalVector    = n / gpt.integrationArea;
                }
                else {
                    // cross product
                    typedef Eigen::Ref<const Vector3d> vector3_cr;
                    Vector3d n          = vector3_cr( gpt.J.col( 0 ) ).cross( vector3_cr( gpt.J.col( 1 ) ) );
                    gpt.integrationArea = n.norm();
                    gpt.normalVector    = n / gpt.integrationArea;
                }

                gaussPts.push_back( std::move( gpt ) );
            }
        }

        VectorXd BoundaryElement::computeNormalLoadVector()
        {
            /* compute the load vector for a constant distributed load (e.g. pressure)
             * Attention: result =  boundary-element-sized!
             *  -> use expandBoundaryToParentVector to obtain the parent-element-sized load vector
             * */

            VectorXd Pk = VectorXd::Zero( coordinates.size() );

            for (const auto& gPt : gaussPts)
                Pk += gPt.integrationArea * gPt.weight * gPt.normalVector.transpose() * NB( gPt.N, nDim );

            return Pk;
        }

        MatrixXd BoundaryElement::computeNormalLoadVectorStiffness()
        {
            MatrixXd K = MatrixXd::Zero( coordinates.size(), coordinates.size() );

            if ( nDim == 2 ) {
                Matrix2D R;
                R << 0, 1, 
                    -1, 0;

                for ( const auto& gPt : gaussPts )
                    for ( int I = 0; I < nNodes; I++ )
                        for ( int J = 0; J < nNodes; J++ )
                            K.block( I*2, J*2, I*2 + 2, J*2 + 2 ) -= gPt.N( I ) * gPt.dNdXi( J ) * R * gPt.weight;
            }

            else if ( nDim == 3 ) {
                Matrix3d Hxi0, Hxi1;
                for ( const auto& gPt : gaussPts ) {
                    HXi0 << 0, J( 2, 0 ), -J( 1, 0 ), 
                            -J( 2, 0 ), 0, J( 0, 0 ), 
                            J( 1, 0 ), -J( 0, 0 ), 0;

                    HXi1 << 0, J( 2, 1 ), -J( 1, 1 ), 
                         -J( 2, 1 ), 0, J( 0, 1 ), 
                         J( 1, 1 ), -J( 0, 1 ), 0;

                    for ( int I = 0; I < nNodes; I++ )
                        for ( int J = 0; J < nNodes; J++ )
                            K.block( I*3, J*3, I*3 + 3, J*3 + 3 ) -= gPt.N( I ) *
                                                                      ( gPt.dNdXi( 0, J ) * HXi0 -
                                                                      gPt.dNdXi( 1, J ) * HXi1 ) *
                                                                      gPt.weight;
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

        VectorXd BoundaryElement::expandBoundaryToParentVector( const VectorXd& boundaryVector )
        {
            /* expand any boundary vector (e.g. pressure load) to the corresponding parental vector
             * */

            VectorXd parentVector( nParentCoordinates );
            parentVector.setZero();

            for ( int i = 0; i < boundaryVector.size(); i++ )
                parentVector( boundaryIndicesInParentCoordinates( i ) ) = boundaryVector( i );

            return parentVector;
        }

    } // namespace FiniteElement
} // namespace bft
