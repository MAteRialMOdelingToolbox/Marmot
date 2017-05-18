#pragma once 
#include <array>
#include "bftTypedefs.h"

namespace bft{

	//****************************************************
	namespace FiniteElement
    {
        MatrixXd NB(const Ref<const VectorXd>& N, const int nDoFPerNode);

        namespace Spatial2d
        {
            static constexpr int nDim = 2;
            
            namespace Quad4
            {
                static constexpr int nNodes = 4;
                
                Vector4d N(const Ref<const Vector2d>& xi);           
                Matrix<double, nNodes,nDim> dNdXi(const Ref <const Vector2d>& xi);
                Vector4d get2DCoordinateIndicesOfBoundaryTruss(int elementFace);
            } 
            
            namespace Quad8
            {
                static constexpr int nNodes = 8;
                
                Matrix<double, nNodes,1> N(const Ref<const Vector2d>& xi);           
                Matrix<double, nNodes,nDim> dNdXi(const Ref <const Vector2d>& xi);
                std::array<int,3> getNodesOfFace(int elementFace);   
                Vector6 get2DCoordinateIndicesOfBoundaryTruss(int elementFace);
            } 
            
            namespace Truss2
            {
                static constexpr int nNodes = 2;
                
                //typedef Eigen::Matrix<double, nDim , nDim * nNodes> Nb_SizedMatrix;
                
                Vector2d N(double  xi);
                //Nb_SizedMatrix Nb(const Vector2d& N);
                Matrix<double, nNodes, nDim> dNdXi(const double& xi);
                
                Vector2d Jacobian(const Vector2d& dNdXi, const Vector4d& coordinates);

                Vector2d TangentialVector(const Vector2d& Jacobian);
                Eigen::Vector2d NormalVector(const Vector2d& Jacobian);
            }
            
            namespace Truss3
            {
                static constexpr int nNodes = 3;
                
                //typedef Eigen::Matrix<double, nDim , nDim * nNodes> Nb_SizedMatrix;
                
                Vector3d N(double  xi);
                //Nb_SizedMatrix Nb(const Vector3d& N);
                Vector3d dNdXi(double xi);

                Vector2d Jacobian(const Vector3d& dNdXi, const Vector6& coordinates);

                Vector2d TangentialVector(const Vector2d& Jacobian);
                Vector2d NormalVector(const Vector2d& Jacobian);
            }    

            namespace Boundary2 
            {
                const Matrix2d gaussPts1d_2(int elementFace);
                Vector4d dNdXi(int elementFace, const Ref<const Vector2d>& xi);
            } 

        }//end of namespace Spatial2D

        namespace NumIntegration
        {
            const Vector2d gaussPts1d_2();
            const Matrix<double, 4, 2> gaussPts2d_2x2();
        } 
    }
} 

