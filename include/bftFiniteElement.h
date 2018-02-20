#pragma once 
#include <array>
#include <vector>
#include "bftTypedefs.h"

namespace bft{

    namespace FiniteElement
    {
        enum ElementShapes{
            Truss2,
            Truss3,
            Quad4,
            Quad8,
            Hexa8,
            Hexa20,
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

        VectorXd expandNodeIndicesToCoordinateIndices( VectorXd nodeIndices, int nDim);


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

                // convenience functions; they are wrappers to the corresponding template functions
                Matrix2d Jacobian(const Ref<const dNdXiSized >& dNdXi, const Ref<const Matrix<double, nNodes*nDim, 1>>& coordinates);
                BSized B(const Ref<const dNdXiSized>& dNdXi);

                Vector2d getBoundaryElementIndices ( int faceID);
            } 

            namespace Quad8
            {
                constexpr int nNodes = 8;

                typedef Matrix<double, 1,       nNodes>         NSized;
                typedef Matrix<double, nDim,    nNodes>         dNdXiSized;
                typedef Matrix<double, 3, nNodes * nDim>        BSized;

                NSized N(const Ref<const Vector2d>& xi);           
                dNdXiSized dNdXi(const Ref<const Vector2d>& xi);

                // convenience functions; they are wrappers to the corresponding template functions
                Matrix2d Jacobian(const Ref<const dNdXiSized >& dNdXi, const Ref<const Matrix<double, nNodes*nDim, 1>>& coordinates);
                BSized B(const Ref<const dNdXiSized>& dNdXi);

                Vector3d getBoundaryElementIndices ( int faceID);
            } 

            namespace Truss2
            {
                constexpr int nNodes = 2;
                typedef Matrix<double, 1,       nNodes>         NSized;
                typedef Matrix<double, 1,    nNodes>         dNdXiSized;

                NSized N(double  xi);
                dNdXiSized dNdXi(double xi);

            }

            namespace Truss3
            {
                constexpr int nNodes = 3;
                typedef Matrix<double, 1,       nNodes>         NSized;
                typedef Matrix<double, 1,    nNodes>         dNdXiSized;

                NSized N(double  xi);
                dNdXiSized dNdXi(double xi);

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
                        B_(2, nDim*i +2) =     dNdX(2, i);
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

                Vector4d getBoundaryElementIndices ( int faceID );
            }

            namespace Hexa20
            {
                constexpr int nNodes = 20;
                typedef Matrix<double, 1,       nNodes>         NSized;
                typedef Matrix<double, nDim,    nNodes>         dNdXiSized;
                typedef Matrix<double, 6,       nNodes * nDim>  BSized;

                NSized N(const Ref<const Vector3d>& xi);           
                dNdXiSized dNdXi(const Ref<const Vector3d>& xi);

                // convenience functions; they are wrappers to the corresponding template functions
                Matrix3d Jacobian(const Ref<const dNdXiSized >& dNdXi, const Ref<const Matrix<double, nNodes*nDim, 1>>& coordinates);
                BSized B(const Ref<const dNdXiSized>& dNdXi);

                Matrix<double, 8, 1 > getBoundaryElementIndices ( int faceID );
            }
        } // End Spatial3D

        class BoundaryElement{

            /* Boundary element, for instance for distributed surface loads
             * */

            struct BoundaryElementGaussPt{
                double weight;
                VectorXd xi ;
                VectorXd N;
                MatrixXd dNdXi;
                MatrixXd J;

                VectorXd normalVector;
                double integrationArea;
            };

            const int            nDim;

            ElementShapes           boundaryShape;
            int                     nNodes;
            int                     nParentCoordinates;

            std::vector < BoundaryElementGaussPt>  gaussPts;

            VectorXd boundaryIndicesInParentNodes;
            VectorXd boundaryIndicesInParentCoordinates;
            VectorXd coordinates;

            public:

            BoundaryElement (ElementShapes parentShape, 
                            int nDim, 
                            int parentFaceNumber,
                            const Ref<const VectorXd>& parentCoordinates
                            );

            VectorXd computeNormalLoadVector ();
            VectorXd condenseParentToBoundaryVector(const Ref<const VectorXd>& parentVector);
            VectorXd expandBoundaryToParentVector(const Ref<const VectorXd>& boundaryVector);

        };
    }

    namespace NumIntegration
    {
        constexpr double gp2 = 0.577350269189625764509;
        constexpr double gp3 = 0.774596669241483;

        enum IntegrationTypes{
            FullIntegration,
            ReducedIntegration
        };

        MatrixXd getGaussPointList(bft::FiniteElement::ElementShapes shape, IntegrationTypes integrationType);
        VectorXd getGaussWeights (bft::FiniteElement::ElementShapes shape, IntegrationTypes integrationType);
        int getNumGaussPoints(bft::FiniteElement::ElementShapes shape, IntegrationTypes integrationType);

        namespace Spatial1D
        {
            constexpr int nDim = 1;

            constexpr double gp2 = 0.577350269189625764509;
            constexpr double gp3 = 0.774596669241483;

            const Vector2d gaussPtList2 = (Vector2d() << -gp2,  +gp2).finished();
            const Vector2d gaussPtList2Weights = (Vector2d() << 1,   1).finished();

            const Vector3d gaussPtList3 = (Vector3d() << -gp3,  0,  +gp3).finished();
            const Vector3d gaussPtList3Weights = (Vector3d() << 5./9,  8./9,  5./9).finished();

        }

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
                    64./81, 25./81, 25./81, 
                    25./81, 25./81, 40./81, 
                    40./81, 40./81, 40./81).finished();

            void modifyCharElemLengthAbaqusLike(double& charElemLength, int intPoint);
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

            const Matrix<double, 3*3*3, nDim> gaussPtList3x3x3 = (Matrix<double, 3*3*3, nDim>() <<  
                    -gp3,  -gp3, -gp3,
                    0,     -gp3, -gp3,
                    +gp3,  -gp3, -gp3, -gp3,  0,    -gp3,
                    0,     0,    -gp3,
                    gp3,   0,    -gp3,
                    -gp3,  +gp3, -gp3,
                    0,     +gp3, -gp3,
                    +gp3,  +gp3, -gp3,

                    -gp3,  -gp3, 0,
                    0,     -gp3, 0,
                    +gp3,  -gp3, 0,
                    -gp3,  0,    0,
                    0,     0,    0,
                    gp3,   0,    0,
                    -gp3,  +gp3, 0,
                    0,     +gp3, 0,
                    +gp3,  +gp3, 0,

                    -gp3,  -gp3, +gp3,
                    0,     -gp3, +gp3,
                    +gp3,  -gp3, +gp3,
                    -gp3,  0,    +gp3,
                    0,     0,    +gp3,
                    gp3,   0,    +gp3,
                    -gp3,  +gp3, +gp3,
                    0,     +gp3, +gp3,
                    +gp3,  +gp3, +gp3

                        ).finished();

            const Matrix<double, 3*3*3,1>   gaussPtList3x3x3Weights = (Matrix<double,3*3*3,1>() <<   
                    0.171467764060357,  0.274348422496571,  0.171467764060357,
                    0.274348422496571,  0.438957475994513,  0.274348422496571,
                    0.171467764060357,  0.274348422496571,  0.171467764060357,
                    0.274348422496571,  0.438957475994513,  0.274348422496571,
                    0.438957475994513,  0.702331961591221,  0.438957475994513,
                    0.274348422496571,  0.438957475994513,  0.274348422496571,
                    0.171467764060357,  0.274348422496571,  0.171467764060357,
                    0.274348422496571,  0.438957475994513,  0.274348422496571,
                    0.171467764060357,  0.274348422496571,  0.171467764060357
                    ).finished();
        }
    } 
} 

