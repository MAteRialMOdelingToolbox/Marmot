#pragma once 
#include "bftTypedefs.h"

namespace bft{

	//****************************************************
	namespace FiniteElement
    {
        namespace Quad4
        {
            static constexpr int nNodes = 4;
            static constexpr int nDim = 2;
            
            Matrix<double, nNodes, 1> shapeFunctions(const Ref<const Vector2d>& xi);           
            Matrix<double, nNodes, nDim> dNdXi(const Ref <const Vector2d>& xi);
        } 
        namespace Boundary2 
        {
            const Matrix2d gaussPts1d_2(int elementFace);
            Vector4d dNdXi(int elementFace, const Ref<const Vector2d>& xi);
        }
        namespace Spatial2D
        {
            static constexpr int nDim = 2;

            namespace Truss2{
                static constexpr int nNodes = 2;

                Vector2d N(double  xi);
                typedef Eigen::Matrix<double, nDim , nDim * nNodes> Nb_SizedMatrix;
                Nb_SizedMatrix Nb(const Vector2d& N);
                Vector2d dNdXi(double xi);

                Vector2d Jacobian(const Vector2d& dNdXi, const Vector4d& coordinates);

                Vector2d TangentialVector(const Vector2d& Jacobian);
                Eigen::Vector2d NormalVector(const Vector2d& Jacobian);
            }

            namespace Truss3{
                static constexpr int nNodes = 3;

                Vector3d N(double  xi);
                typedef Eigen::Matrix<double, nDim , nDim * nNodes> Nb_SizedMatrix;
                Nb_SizedMatrix Nb(const Vector3d& N);
                Vector3d dNdXi(double xi);

                Vector2d Jacobian(const Vector3d& dNdXi, const Vector6& coordinates);

                Vector2d TangentialVector(const Vector2d& Jacobian);
                Vector2d NormalVector(const Vector2d& Jacobian);

            }
            namespace Quad4{
                static constexpr int nNodes = 4;
                Vector4d get2DCoordinateIndicesOfBoundaryTruss(int elementFace);
                
            }
            namespace Quad8{
                static constexpr int nNodes = 8;

                //typedef Eigen::Matrix<double, nDim , nDim * nNodes> Nb_SizedMatrix;
                Vector6 get2DCoordinateIndicesOfBoundaryTruss(int elementFace);
            }

        }//end of namespace Spatial2D

        namespace NumIntegration
        {
            const Eigen::Matrix<double, 4, 2> gaussPts2d_2x2();
        } 
    } 
} 

