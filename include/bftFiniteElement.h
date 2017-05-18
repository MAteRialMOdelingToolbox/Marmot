#include "bftVoigt.h"
#include "bftFunctions.h"
#include "bftConstants.h"
#include <iostream>
#include <array>

namespace bft{

	//****************************************************
	namespace FiniteElement
    {
        MatrixXd createNBold(const Ref<const VectorXd>& N, const int nDoFPerNode);

        namespace Quad4
        {
            static constexpr int nNodes = 4;
            static constexpr int nDim = 2;
            
            Matrix<double, nNodes, 1> shapeFunctions(const Ref<const Vector2d>& xi);           
            Matrix<double, nNodes, nDim> dNdXi(const Ref <const Vector2d>& xi);
        } 
        namespace Quad8
        {
            static constexpr int nNodes = 8;
            static constexpr int nDim = 2;
            
            Matrix<double, nNodes, 1> shapeFunctions(const Ref<const Vector2d>& xi);           
            Matrix<double, nNodes, nDim> dNdXi(const Ref <const Vector2d>& xi);
            std::array<int,3> getNodesOfFace(int elementFace);   
        } 
        namespace Truss2
        {
            static constexpr int nNodes = 2;
            static constexpr int nDim = 1;
            
            Matrix<double, nNodes, 1> shapeFunctions(const double& xi);           
            Matrix<double, nNodes, nDim> dNdXi(const double& xi);
        }
        
        namespace Truss3
        {
            static constexpr int nNodes = 3;
            static constexpr int nDim = 1;
            
            Matrix<double, nNodes, 1> shapeFunctions(const double& xi);           
            Matrix<double, nNodes, nDim> dNdXi(const double& xi);
        }    

        namespace Boundary2 
        {
            const Matrix2d gaussPts1d_2(int elementFace);
            Vector4d dNdXi(int elementFace, const Ref<const Vector2d>& xi);
        } 
    }

    namespace NumIntegration
        {
            const Vector2d gaussPts1d_2();
            const Matrix<double, 4, 2> gaussPts2d_2x2();
        } 
} 

