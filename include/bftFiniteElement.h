#pragma once 
#include <array>
#include "bftTypedefs.h"

namespace bft{

	namespace FiniteElement
    {
        enum ElementShapes{
            Quad4,
            Quad8,
            Hexa8,
            Hexa30,
        };

        ElementShapes getElementShapeByMetric(int nDim, int nNodes);

        // 'Expanded' N , aka NBold aka multidimensional Interpolation Operator
        MatrixXd NB(const Ref<const VectorXd>& N, const int nDoFPerNode); // Dynamic version

        template <int nDim, int nNodes>
            Matrix<double, nDim, nDim*nNodes> NB(const Ref<const Matrix<double, nNodes, 1>>& N )
            {
                // Alternative Templated version of Interpolation operator NBold;
                Matrix<double, nDim, nDim*nNodes> N_ = Matrix<double, nDim, nDim*nNodes>::Zero();
                for (int i=0; i<nNodes; i++){
                    for (int j=0; j<nDim; j++){
                        N_(j,nDim*i+j) = N(i);
                    }
                }
                return N_;
            }
        

        MatrixXd Jacobian(const MatrixXd& dN_dXi, const VectorXd& coordinates); // Dynamic version

        template <int nDim, int nNodes> 
            Matrix<double, nDim, nDim> Jacobian(const Matrix<double, nDim, nNodes>&         dNdXi, 
                                                const Matrix<double, nDim * nNodes, 1>&     coordinates)
        {
            // Alternative Templated version of Jacobian for compile time known sizes
                Matrix<double, nDim, nDim> J_ = Matrix<double, nDim, nDim>::Zero();
                for(int i = 0; i < nDim; i++)		// loop over global dimensions
                    for(int j=0; j < nDim; j++)		// loop over local dimensions
                        for(int k=0; k<nNodes; k++) // Loop over nodes
                            J_(i, j) += dNdXi(j, k) * coordinates(i + k*nDim);
                return J_;
        }

        namespace Spatial2D
        {
            constexpr int nDim = 2;
            constexpr int voigtLength=3;

            template<int nNodes> 
                Matrix<double, voigtLength, nNodes*nDim> B(const Ref<const  Matrix<double, nDim, nNodes>>& dNdX) {   

                    Matrix<double, voigtLength, nNodes*nDim> B_ = Matrix<double, voigtLength, nNodes*nDim>::Zero();
                    for(int i = 0; i < nNodes; i++){
                        B_(0, 2*i) =        dNdX(0, i);
                        B_(1, 2*i +1) =     dNdX(1, i);
                        B_(2, 2*i) =        dNdX(1, i);
                        B_(2, 2*i +1) =     dNdX(0, i);}  
                    return B_;
                }
            
            namespace Quad4
            {
                constexpr int nNodes = 4;

                typedef Matrix<double, 1,       nNodes>         NSized;
                typedef Matrix<double, nDim,    nNodes>         dNdXiSized;
                typedef Matrix<double, 3, nNodes * nDim>        BSized;

                NSized N(const Ref<const Vector2d>& xi);           
                dNdXiSized dNdXi(const Ref<const Vector2d>& xi);
                NSized get2DCoordinateIndicesOfBoundaryTruss(int elementFace);
              
                // convenience functions; they are wrappers to the corresponding template functions
                Matrix2d Jacobian(const Ref<const dNdXiSized >& dNdXi, const Ref<const Matrix<double, nNodes*nDim, 1>>& coordinates);
                BSized B(const Ref<const dNdXiSized>& dNdXi);

                namespace Boundary2 
                {
                    const Matrix2d gaussPts1d_2(int elementFace);
                    NSized dNdXi(int elementFace, const Ref<const Vector2d>& xi);
                } 
            } 
            
            namespace Quad8
            {
                constexpr int nNodes = 8;

                typedef Matrix<double, 1,       nNodes>         NSized;
                typedef Matrix<double, nDim,    nNodes>         dNdXiSized;
                typedef Matrix<double, 3, nNodes * nDim>        BSized;
                
                NSized N(const Ref<const Vector2d>& xi);           
                dNdXiSized dNdXi(const Ref<const Vector2d>& xi);
                std::array<int,3> getNodesOfFace(int elementFace);   
                Vector6 get2DCoordinateIndicesOfBoundaryTruss(int elementFace);

                // convenience functions; they are wrappers to the corresponding template functions
                Matrix2d Jacobian(const Ref<const dNdXiSized >& dNdXi, const Ref<const Matrix<double, nNodes*nDim, 1>>& coordinates);
                BSized B(const Ref<const dNdXiSized>& dNdXi);
            } 
            
            namespace Truss2
            {
                constexpr int nNodes = 2;
                
                Vector2d N(double  xi);
                Vector2d dNdXi(double xi);
                
                Vector2d Jacobian(const Vector2d& dNdXi, const Vector4d& coordinates);

                Vector2d TangentialVector(const Vector2d& Jacobian);
                Vector2d NormalVector(const Vector2d& Jacobian);
            }
            
            namespace Truss3
            {
                constexpr int nNodes = 3;
                
                Vector3d N(double  xi);
                Vector3d dNdXi(double xi);

                Vector2d Jacobian(const Vector3d& dNdXi, const Vector6& coordinates);

                Vector2d TangentialVector(const Vector2d& Jacobian);
                Vector2d NormalVector(const Vector2d& Jacobian);
            }    
        }//end of namespace Spatial2D

        namespace Spatial3D
        {
            constexpr int nDim = 3;
            constexpr int voigtLength = 6;

            template<int nNodes> 
                Matrix<double, voigtLength, nNodes*nDim> B(const Ref<const  Matrix<double, nDim, nNodes>>& dNdX) {   

                    Matrix<double, voigtLength, nNodes*nDim> B_ = Matrix<double, voigtLength, nNodes*nDim>::Zero();
                    for(int i = 0; i < nNodes; i++){
                        B_(0, nDim*i) =        dNdX(0, i);
                        B_(1, nDim*i +1) =     dNdX(1, i);
                        B_(2, nDim*i +2) =     dNdX(1, i);
                        B_(3, nDim*i +0) =     dNdX(1, i);
                        B_(3, nDim*i +1) =     dNdX(0, i);
                        B_(4, nDim*i +1) =     dNdX(2, i);
                        B_(4, nDim*i +2) =     dNdX(1, i);
                        B_(5, nDim*i +2) =     dNdX(0, i);
                        B_(5, nDim*i +0) =     dNdX(2, i); }  

                    return B_;
                }
 
            namespace Hexa8
            {
                constexpr int nNodes = 8;
                typedef Matrix<double, 1,       nNodes>         NSized;
                typedef Matrix<double, nDim,    nNodes>         dNdXiSized;
                typedef Matrix<double, 6,       nNodes * nDim>  BSized;

                NSized N(const Ref<const Vector3d>& xi);           
                dNdXiSized dNdXi(const Ref<const Vector3d>& xi);
               
                // convenience functions; they are wrappers to the corresponding template functions
                Matrix3d Jacobian(const Ref<const dNdXiSized >& dNdXi, const Ref<const Matrix<double, nNodes*nDim, 1>>& coordinates);
                BSized B(const Ref<const dNdXiSized>& dNdXi);
            }
        }
    }

    namespace NumIntegration
    {
        const Vector2d gaussPts1d_2();
        const Matrix<double, 4, 2> gaussPts2d_2x2();

        constexpr double gp2 = 0.577350269189625764509;
        constexpr double gp3 = 0.774596669241483;

        enum IntegrationTypes{
            FullIntegration,
            ReducedIntegration
        };

        MatrixXd getGaussPointList(bft::FiniteElement::ElementShapes shape, IntegrationTypes integrationType);
        VectorXd getGaussWeights (bft::FiniteElement::ElementShapes shape, IntegrationTypes integrationType);

        
        namespace Spatial2D
        { 
            constexpr int nDim = 2;

            const RowVector2d gaussPtList1x1 = RowVector2d :: Zero();
            const Matrix<double,1,1> gaussPtList1x1Weights = (Matrix<double,1,1>() << 
                    4).finished();

            using Matrix42 = Matrix<double, 4, nDim>;
            const Matrix42 gaussPtList2x2 = (Matrix42() <<  
                    +gp2,  +gp2, 
                    -gp2,  +gp2,
                    -gp2,  -gp2, 
                    +gp2,  -gp2).finished();
            const Matrix<double, 4,1>   gaussPtList2x2Weights = (Matrix<double,4,1>() << 
                    1,1,1,1).finished();

            using Matrix92 = Matrix<double, 9, nDim>;
            const Matrix92 gaussPtList3x3 = (Matrix92() <<  
                    0,     0,
                    -gp3,  -gp3,
                    +gp3,  -gp3,
                    +gp3,  +gp3,
                    -gp3,  +gp3,
                    0,     -gp3,
                    gp3,   0,
                    0,     +gp3,
                    -gp3,  0).finished();
            const Matrix<double, 9,1>   gaussPtList3x3Weights = (Matrix<double,9,1>() <<   
                    64./81, 25./91, 25./81, 
                    25./81, 25./81, 40./81, 
                    40./81, 40./81, 40./81).finished();

        }

        namespace Spatial3D
        { 
            constexpr int nDim = 3;

            const RowVector3d gaussPtList1x1x1 = RowVector3d::Zero();
            const Matrix<double, 1, 1> gaussPtList1x1x1Weights = 
                (Matrix<double, 1, 1>() << 
                8.0
                ).finished();

            const Matrix<double, 8, nDim> gaussPtList2x2x2 = (
            Matrix<double, 8, nDim>() <<  
               -gp2,    -gp2,     -gp2,
               +gp2,    -gp2,     -gp2,
               +gp2,    +gp2,     -gp2,
               -gp2,    +gp2,     -gp2,
               -gp2,    -gp2,     +gp2,
               +gp2,    -gp2,     +gp2,
               +gp2,    +gp2,     +gp2,
               -gp2,    +gp2,     +gp2).finished();

            const Matrix<double, 8,1>   gaussPtList2x2x2Weights = (Matrix<double,8,1>() <<   
                    1,1,1,1,1,1,1,1).finished();
        }

    } 
} 

