#include "bftFiniteElement.h"
#include "bftTypedefs.h"
#include <iostream>

namespace bft{

    //****************************************************
    namespace FiniteElement
    {
        ElementShapes getElementShapeByMetric(int nDim, int nNodes){

            switch(nDim){

                case (2):{ switch(nNodes){
                             case 4:    { return ElementShapes::Quad4;}
                             case 8:    { return ElementShapes::Quad8;}
                             default:   {throw std::invalid_argument("Invalid number of nodes for nDim=2" );}
                         }}
                case (3):{ switch(nNodes){
                             case 8:    { return ElementShapes::Hexa8;}
                             case 20:   { return ElementShapes::Hexa20;}
                             default:   {throw std::invalid_argument("Invalid number of nodes for nDim=3" );}
                         }}
                default:{throw std::invalid_argument("Invalid dimension specified");}
            }
        }
        // create a large N matrix dependent on the degrees of freedom per node
        MatrixXd NB(const Ref<const VectorXd>& N, const int nDoFPerNode){

            MatrixXd N_(nDoFPerNode, N.size()*nDoFPerNode);
            N_ = MatrixXd::Zero(nDoFPerNode, N.size()*nDoFPerNode);

            for (int i=0; i<N.size(); i++){
                for (int j=0; j<nDoFPerNode; j++){
                    N_(j,nDoFPerNode*i+j) = N(i);
                }
            }

            return N_;} 

        MatrixXd Jacobian(const MatrixXd& dNdXi, const VectorXd& coordinates)
        {
            /* Notation:
             *
             * /                          \
             * | x1,xi1,  x1,xi2,  x1,xi3 |
             * |                          |
             * | x2,xi1,  x2,xi2,  x2,xi3 |
             * |                          |
             * | x3,xi1,  x3,xi2,  x3,xi3 |
             * \                          /
             *
             * */

            int nDim = dNdXi.rows();
            int nNodes = dNdXi.cols();

            MatrixXd J_ = MatrixXd::Zero(nDim, nDim);

            for(int i = 0; i < nDim; i++)		// loop over global dimensions
                for(int j=0; j < nDim; j++)		// loop over local dimensions
                    for(int k=0; k<nNodes; k++) // Loop over nodes
                        J_(i, j) += dNdXi(j, k) * coordinates(i + k*nDim);
            return J_;
        }

        namespace BoundaryElementFactory{

            using bft::FiniteElement::ElementShapes;

            VectorXd getBoundaryNodeList(bft::FiniteElement::ElementShapes shape, const int& elementFace){
                switch(shape)
                {
                    case(ElementShapes::Quad4):  { using namespace FiniteElement::Spatial2D::Quad4;
                                                     return get2DCoordinateIndicesOfBoundaryTruss(elementFace); }
                    case(ElementShapes::Quad8):  { using namespace FiniteElement::Spatial2D::Quad8;
                                                     return get2DCoordinateIndicesOfBoundaryTruss(elementFace); } 
                    default: {throw std::invalid_argument("NodeIdxList: Invalid shape combination for boundary element");}
                }
            }

            MatrixXd getNormalVector(bft::FiniteElement::ElementShapes shape, const Ref<const VectorXd>& coords,  const Ref<const VectorXd>& gp){
                switch(shape)
                {
                    case(ElementShapes::Quad4): { using namespace FiniteElement::Spatial2D::Truss2;
                                                    return NormalVector(Jacobian(dNdXi(gp(0)), coords)); }
                    case(ElementShapes::Quad8): { using namespace FiniteElement::Spatial2D::Truss3; 
                                                    return NormalVector(Jacobian(dNdXi(gp(0)), coords)); }
                    default: {throw std::invalid_argument("Jacobian: Invalid shape combination for boundary element");}
                }
            }

            MatrixXd getGaussPointList(bft::FiniteElement::ElementShapes shape)
            {
                switch(shape)
                {
                    case(ElementShapes::Quad4): {  return FiniteElement::BoundaryElementFactory::gaussPtList2; }
                    case(ElementShapes::Quad8): {  return FiniteElement::BoundaryElementFactory::gaussPtList3; } 
                    default: {throw std::invalid_argument("Gausspoints: Invalid shape/integrationType combination boundary eement");}
                }
            }

            VectorXd getGaussWeights(bft::FiniteElement::ElementShapes shape)
            {
                switch(shape)
                {
                    case(ElementShapes::Quad4): {  return FiniteElement::BoundaryElementFactory::gaussPtList2Weights; }
                    case(ElementShapes::Quad8): {  return FiniteElement::BoundaryElementFactory::gaussPtList3Weights; } 
                    default: {throw std::invalid_argument("Boundary element: invalid gauss weights");}
                }
            }

            MatrixXd getNB(bft::FiniteElement::ElementShapes shape, const Ref<const VectorXd>& gp)
            {
                using namespace bft::FiniteElement;
                switch(shape)
                {
                    case(ElementShapes::Quad4): {  return NB( Spatial2D::Truss2::N(gp(0)), 2); }
                    case(ElementShapes::Quad8): {  return NB( Spatial2D::Truss3::N(gp(0)), 2); } 
                    default: {throw std::invalid_argument("Boundary element: invalid NB shape");}
                }
            }

            double getIntVol(bft::FiniteElement::ElementShapes shape, const Ref<const VectorXd>& coords, const Ref<const VectorXd>& gp)
            {
                switch(shape)
                {
                    case(ElementShapes::Quad4): {  using namespace bft::FiniteElement::Spatial2D::Truss2; 
                                                    return Jacobian(dNdXi(gp(0)), coords).norm(); }
                    case(ElementShapes::Quad8): {  using namespace bft::FiniteElement::Spatial2D::Truss3; 
                                                    return Jacobian(dNdXi(gp(0)), coords).norm(); } 
                    default: {throw std::invalid_argument("Boundary element: invalid integration volume");}
                }
            }
        }// end of namespace BoundaryElementFactory
    }

    namespace NumIntegration
    {
        MatrixXd getGaussPointList(bft::FiniteElement::ElementShapes shape, IntegrationTypes integrationType)
        {
            using bft::FiniteElement::ElementShapes;
            switch(shape)
            {
                case(ElementShapes::Quad4): {   if(integrationType == IntegrationTypes::FullIntegration)
                                                return Spatial2D::gaussPtList2x2;
                                                else
                                                    return Spatial2D::gaussPtList1x1;}

                case(ElementShapes::Quad8): {   if(integrationType == IntegrationTypes::FullIntegration)
                                                return Spatial2D::gaussPtList3x3;
                                                else
                                                    return Spatial2D::gaussPtList2x2;}

                case(ElementShapes::Hexa8): {   if(integrationType == IntegrationTypes::FullIntegration)
                                                return Spatial3D::gaussPtList2x2x2;
                                                else
                                                    return Spatial3D::gaussPtList1x1x1;}

                case(ElementShapes::Hexa20): {   if(integrationType == IntegrationTypes::FullIntegration)
                                                 return Spatial3D::gaussPtList3x3x3;
                                                 else
                                                     return Spatial3D::gaussPtList2x2x2;}

                default: {throw std::invalid_argument("Invalid shape/integrationType combination");}
            }
        }

        VectorXd getGaussWeights (bft::FiniteElement::ElementShapes shape, IntegrationTypes integrationType){
            using bft::FiniteElement::ElementShapes;
            switch(shape)
            {
                case(ElementShapes::Quad4): {   if(integrationType == IntegrationTypes::FullIntegration)
                                                return Spatial2D::gaussPtList2x2Weights;
                                                else
                                                    return Spatial2D::gaussPtList1x1Weights;}

                case(ElementShapes::Quad8): {   if(integrationType == IntegrationTypes::FullIntegration)
                                                return Spatial2D::gaussPtList3x3Weights;
                                                else
                                                    return Spatial2D::gaussPtList2x2Weights;}

                case(ElementShapes::Hexa8): {   if(integrationType == IntegrationTypes::FullIntegration)
                                                return Spatial3D::gaussPtList2x2x2Weights;
                                                else
                                                    return Spatial3D::gaussPtList1x1x1Weights;}

                case(ElementShapes::Hexa20): {   if(integrationType == IntegrationTypes::FullIntegration)
                                                 return Spatial3D::gaussPtList3x3x3Weights;
                                                 else
                                                     return Spatial3D::gaussPtList2x2x2Weights;}

                default: {throw std::invalid_argument("Invalid shape/integrationType combination");}
            }
        }
    } // end of namespace NumIntegration
} // end of namespace bft
