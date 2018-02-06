#include "bftFiniteElement.h"
#include <iostream>

namespace bft{
    namespace FiniteElement{

        BoundaryElement::BoundaryElement (
                ElementShapes parentShape,
                int parentFaceNumber, 
                //int parentFaceID,
                int nDim, 
                const Ref<const VectorXd>& parentCoordinates
                //const Ref<const MatrixXd>& gaussPointList,
                //const Ref<const VectorXd>& gaussPointWeights
                ):
            //shape(shape),
            nDim(nDim)
            //coordinates( coordinates )
        {
            // compute for each gaussPt the
            // - dNdXi
            // - Jacobian pp(x,xi)
            // - normalVector
            // - integrationArea 
            //
            //ElementShapes shape;

            switch (parentShape)
            {
                case Quad4: {   shape = Truss2; 
                                nNodes = bft::FiniteElement::Spatial2D::Truss2::nNodes;
                                boundaryIndicesInParentNodes = bft::FiniteElement::Spatial2D::Quad4::getBoundaryElementIndices( parentFaceNumber );
                                break;}

                case Quad8: {   shape = Truss3; 
                                nNodes = bft::FiniteElement::Spatial2D::Truss3::nNodes;
                                boundaryIndicesInParentNodes = bft::FiniteElement::Spatial2D::Quad8::getBoundaryElementIndices( parentFaceNumber );
                                break;}
                case Hexa8: {   shape = Hexa8; 
                                nNodes = bft::FiniteElement::Spatial3D::Hexa8::nNodes;
                                boundaryIndicesInParentNodes = bft::FiniteElement::Spatial3D::Hexa8::getBoundaryElementIndices( parentFaceNumber );
                                break;}

                default: throw std::invalid_argument("Boundary Element currently not implemented");

            }

            nParentCoordinates =                    parentCoordinates.size();
            boundaryIndicesInParentCoordinates =    expandNodeIndicesToCoordinateIndices ( boundaryIndicesInParentNodes, nDim );


            coordinates =                   condenseParentToBoundaryVector ( parentCoordinates );
            MatrixXd gaussPointList =       NumIntegration::getGaussPointList( shape , NumIntegration::IntegrationTypes::FullIntegration);
            VectorXd gaussPointWeights =    NumIntegration::getGaussWeights( shape , NumIntegration::IntegrationTypes::FullIntegration);

            for(int i = 0; i < gaussPointList.rows(); i++){

                const VectorXd& xi = gaussPointList.row(i);


                //std::cout << " gpt " << i << xi.transpose() << std::endl;

                BoundaryElementGaussPt gpt ;
                gpt.xi = xi;
                gpt.weight = gaussPointWeights(i);

                // dNdXi
                switch ( shape  ) 
                {
                    case Truss2: {  gpt.N = Spatial2D::Truss2::N ( xi(0 ) ); 
                                     gpt.dNdXi = Spatial2D::Truss2::dNdXi( xi(0) ); 
                                     break;}
                    case Truss3: {  gpt.N = Spatial2D::Truss3::N ( xi(0) ); 
                                     gpt.dNdXi = Spatial2D::Truss3::dNdXi( xi(0) ); 
                                     break;}
                    case Quad4: {   gpt.N =     Spatial2D::Quad4::N ( xi.head(2) ); 
                                    gpt.dNdXi = Spatial2D::Quad4::dNdXi( xi.head(2) ); 
                                    break;}
                    default: break; // exception handling already in first switch
                }
            
                //std::cout << "N"  <<gpt.N.transpose() << std::endl;
                //std::cout << "dNdXi" << gpt.dNdXi.transpose() << std::endl;
                //// Jacobian
                //std::cout << "coordSParent" << coordinates.transpose() << std::endl;
                //std::cout << "coordS" << parentCoordinates.transpose() << std::endl;
                //std::cout << "indicees" << boundaryIndicesInParentNodes.transpose() << std::endl;
                //std::cout << "indices2" << boundaryIndicesInParentCoordinates.transpose() << std::endl;
                //std::cout << "nDim " << nDim << std::endl;
                //std::cout << "nNodes" << nNodes<< std::endl;
                gpt.J = Jacobian(gpt.dNdXi, coordinates); 

                //std::cout << "J" << gpt.J << std::endl;

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
                    gpt.integrationArea = gpt.normalVector.norm();
                    gpt.normalVector = n / gpt.integrationArea;
                }

                //std::cout << "n" <<  gpt.normalVector.transpose() << std::endl;
                //std::cout << "intA" << gpt.integrationArea << std::endl;

                gaussPts.push_back ( std::move (gpt) );
            }

        }

        VectorXd BoundaryElement::computeNormalLoadVector()
        {
            VectorXd Pk = VectorXd::Zero( coordinates.size() );

            for(size_t i = 0; i < gaussPts.size(); i++){
                const BoundaryElementGaussPt& gpt = gaussPts[i];
                Pk += gpt.integrationArea * gpt.weight * gpt.normalVector.transpose() * NB (gpt.N, nDim);
            }

            return Pk;
        }

        VectorXd BoundaryElement::condenseParentToBoundaryVector(const Ref<const VectorXd>& parentVector)
        {
            VectorXd boundaryVector( nNodes * nDim );

            for(int i = 0; i < boundaryIndicesInParentCoordinates.size(); i++)
               boundaryVector(i) = parentVector( boundaryIndicesInParentCoordinates(i) );

            return boundaryVector;
            
        }

        VectorXd BoundaryElement::expandBoundaryToParentVector(const Ref<const VectorXd>& boundaryVector)
        {
            VectorXd parentVector ( nParentCoordinates  );
            parentVector.setZero();

            for(int i = 0; i < boundaryVector.size(); i++)
               parentVector( boundaryIndicesInParentCoordinates(i) ) = boundaryVector ( i );

            return parentVector ;

        }

    }
}

