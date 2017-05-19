#pragma once 
#include <array>
#include "bftTypedefs.h"

namespace bft{

	//****************************************************
	namespace FiniteElement
    {
        MatrixXd NB(const Ref<const VectorXd>& N, const int nDoFPerNode);

        // works, but offers no real advantage
        //template <int dim, int nNodes>
            //Matrix<double, dim, dim*nNodes> Nb(const Ref<const Matrix<double, nNodes, 1>>& N )
            //{
                //Matrix<double, dim, dim*nNodes> N_ = Matrix<double, dim, dim*nNodes>::Zero();
                //for (int i=0; i<N.size(); i++){
                    //for (int j=0; j<dim; j++){
                        //N_(j,dim*i+j) = N(i);
                    //}
                //}

                //return N_;
            //}
        

        namespace Spatial2D
        {
            static constexpr int nDim = 2;
            
            namespace Quad4
            {
                static constexpr int nNodes = 4;
                
                Vector4d N(const Ref<const Vector2d>& xi);           
                Matrix<double, nNodes,nDim> dNdXi(const Ref <const Vector2d>& xi);
                Vector4d get2DCoordinateIndicesOfBoundaryTruss(int elementFace);

                namespace Boundary2 
                {
                    const Matrix2d gaussPts1d_2(int elementFace);
                    Vector4d dNdXi(int elementFace, const Ref<const Vector2d>& xi);
                } 
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
                
                Vector2d N(double  xi);
                Vector2d dNdXi(double xi);
                
                Vector2d Jacobian(const Vector2d& dNdXi, const Vector4d& coordinates);

                Vector2d TangentialVector(const Vector2d& Jacobian);
                Vector2d NormalVector(const Vector2d& Jacobian);
            }
            
            namespace Truss3
            {
                static constexpr int nNodes = 3;
                
                Vector3d N(double  xi);
                Vector3d dNdXi(double xi);

                Vector2d Jacobian(const Vector3d& dNdXi, const Vector6& coordinates);

                Vector2d TangentialVector(const Vector2d& Jacobian);
                Vector2d NormalVector(const Vector2d& Jacobian);
            }    


        }//end of namespace Spatial2D

    }
    namespace NumIntegration
    {
        const Vector2d gaussPts1d_2();
        const Matrix<double, 4, 2> gaussPts2d_2x2();
    } 
} 

