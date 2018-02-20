#include "bftFiniteElement.h"
#include <iostream>

namespace bft{
    namespace FiniteElement{

        BoundaryElement::BoundaryElement (
                ElementShapes parentShape,
                int parentFaceNumber, 
                int nDim, 
                const Ref<const VectorXd>& parentCoordinates
                ):
            nDim(nDim)
        {
            //get boundary shape dependent on parental shape 
            switch (parentShape)
            {
                /* 
                 * Extend for future elements here ... (scroll down) .. 
                 *
                 * */

                case Quad4: {   boundaryShape = Truss2; 
                                nNodes = bft::FiniteElement::Spatial2D::Truss2::nNodes;
                                boundaryIndicesInParentNodes = bft::FiniteElement::Spatial2D::Quad4::getBoundaryElementIndices( parentFaceNumber );
                                break;}

                case Quad8: {   boundaryShape = Truss3; 
                                nNodes = bft::FiniteElement::Spatial2D::Truss3::nNodes;
                                boundaryIndicesInParentNodes = bft::FiniteElement::Spatial2D::Quad8::getBoundaryElementIndices( parentFaceNumber );
                                break;}
                case Hexa8: {   boundaryShape =  Quad4; 
                                nNodes = bft::FiniteElement::Spatial2D::Quad4::nNodes;
                                boundaryIndicesInParentNodes = bft::FiniteElement::Spatial3D::Hexa8::getBoundaryElementIndices( parentFaceNumber );
                                break;}
                case Hexa20: {  boundaryShape =  Quad8; 
                                nNodes = bft::FiniteElement::Spatial2D::Quad8::nNodes;
                                boundaryIndicesInParentNodes = bft::FiniteElement::Spatial3D::Hexa20::getBoundaryElementIndices( parentFaceNumber );
                                break;}

                default: throw std::invalid_argument("Boundary Element currently not implemented");

            }

            // get the 'condensed' boundary element coordinates
            nParentCoordinates =                    parentCoordinates.size();
            boundaryIndicesInParentCoordinates =    expandNodeIndicesToCoordinateIndices ( boundaryIndicesInParentNodes, nDim );
            coordinates =                           condenseParentToBoundaryVector ( parentCoordinates );

            // get the proper gausspoints for the boundary element
            MatrixXd gaussPointList =       NumIntegration::getGaussPointList( boundaryShape , NumIntegration::IntegrationTypes::FullIntegration);
            VectorXd gaussPointWeights =    NumIntegration::getGaussWeights( boundaryShape , NumIntegration::IntegrationTypes::FullIntegration);

            // compute for each gaussPt the
            // - dNdXi
            // - Jacobian pp(x,xi)
            // - normalVector
            // - integrationArea 
            for(int i = 0; i < gaussPointList.rows(); i++){

                const VectorXd& xi = gaussPointList.row(i);

                BoundaryElementGaussPt gpt ;
                gpt.xi = xi;
                gpt.weight = gaussPointWeights(i);
    
                // N
                // dNdXi
                switch ( boundaryShape ) 
                {
                    /* 
                     * ... and extend for future elements here!
                     *
                     * */
                    case Truss2: {  gpt.N = Spatial2D::Truss2::N ( xi(0 ) ); 
                                     gpt.dNdXi = Spatial2D::Truss2::dNdXi( xi(0) ); 
                                     break;}
                    case Truss3: {  gpt.N = Spatial2D::Truss3::N ( xi(0) ); 
                                     gpt.dNdXi = Spatial2D::Truss3::dNdXi( xi(0) ); 
                                     break;}
                    case Quad4: {  
                                    gpt.N =     Spatial2D::Quad4::N ( Ref<const Vector2d> ( xi ) ); 
                                    gpt.dNdXi = Spatial2D::Quad4::dNdXi( Ref<const Vector2d> ( xi ) ); 
                                    break;}
                    case Quad8: {  
                                    gpt.N =     Spatial2D::Quad8::N ( Ref<const Vector2d> ( xi ) ); 
                                    gpt.dNdXi = Spatial2D::Quad8::dNdXi( Ref<const Vector2d> ( xi ) ); 
                                    break;}

                    default: break; // exception handling already in first switch
                }
            
                // Jacobian
                gpt.J = Jacobian(gpt.dNdXi, coordinates); 

                // normalVector and integration area
                if(nDim == 2)
                {
                    // 90deg rotation 
                    Vector2d n; n <<  gpt.J(1), -gpt.J(0); 
                    gpt.integrationArea = n.norm(); 
                    gpt.normalVector = n / gpt.integrationArea;

                }
                else
                {
                    // cross product
                    typedef Eigen::Ref<const Vector3d> vector3_cr;
                    Vector3d n  = vector3_cr ( gpt.J.col(0)) .cross ( vector3_cr( gpt.J.col(1)) );
                    gpt.integrationArea = n.norm();
                    gpt.normalVector = n / gpt.integrationArea;
                }

                gaussPts.push_back ( std::move (gpt) );
            }

        }

        VectorXd BoundaryElement::computeNormalLoadVector()
        {
            /* compute the load vector for a constant distributed load (e.g. pressure) 
             * Attention: result =  boundary-element-sized! 
             *  -> use expandBoundaryToParentVector to obtain the parent-element-sized load vector 
             * */

            VectorXd Pk = VectorXd::Zero( coordinates.size() );

            for(size_t i = 0; i < gaussPts.size(); i++){
                const BoundaryElementGaussPt& gpt = gaussPts[i];
                Pk += gpt.integrationArea * gpt.weight * gpt.normalVector.transpose() * NB (gpt.N, nDim);
            }

            return Pk;
        }

        VectorXd BoundaryElement::condenseParentToBoundaryVector(const Ref<const VectorXd>& parentVector)
        {
            /*  condense any parent vector to the corresponding boundary child vector (e.g. coordinates)
             *  dependent on the underlying indices mapping
             * */
            VectorXd boundaryVector( nNodes * nDim );

            for(int i = 0; i < boundaryIndicesInParentCoordinates.size(); i++)
               boundaryVector(i) = parentVector( boundaryIndicesInParentCoordinates(i) );

            return boundaryVector;
            
        }

        VectorXd BoundaryElement::expandBoundaryToParentVector(const Ref<const VectorXd>& boundaryVector)
        {   
            /* expand any boundary vector (e.g. pressure load) to the corresponding parental vector
             * */

            VectorXd parentVector ( nParentCoordinates  );
            parentVector.setZero();

            for(int i = 0; i < boundaryVector.size(); i++)
               parentVector( boundaryIndicesInParentCoordinates(i) ) = boundaryVector ( i );

            return parentVector ;

        }

    }
}

