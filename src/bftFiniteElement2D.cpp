#include "bftFiniteElement.h"
#include "bftTypedefs.h"
#include <iostream>

namespace bft{

    //****************************************************
    namespace FiniteElement
    {
        namespace Spatial2D
        {
            namespace Quad4
            {

                NSized N(const Ref<const Vector2d>& xi){

                    /* 
                     *   Shape functions
                     *    (1) _______(0)
                     *       |       |  
                     *       |       |  
                     *       |_______|  
                     *    (2)        (3) 
                     */

                    NSized N_;
                    N_ <<   1./4 * (1.+xi(0))*(1.+xi(1)),
                       1./4 * (1.-xi(0))*(1.+xi(1)),
                       1./4 * (1.-xi(0))*(1.-xi(1)),
                       1./4 * (1.+xi(0))*(1.-xi(1));
                    return N_;
                }

                dNdXiSized dNdXi(const Ref<const Vector2d>& xi){

                    dNdXiSized result;

                    result <<   /* 1                2                 3                 4*/
                        /* ,xi1 */  +1./4*(1+xi(1)),    -1./4*(1+xi(1)), -1./4*(1-xi(1)), +1./4*(1-xi(1)), 
                        /*, xi2 */  +1./4*(1+xi(0)),    +1./4*(1-xi(0)), -1./4*(1-xi(0)), -1./4*(1+xi(0));

                    return result;}

                Vector2d getBoundaryElementIndices ( int faceID)
                {
                    /*    
                     *               face 1
                     *           (1) _______(0)
                     *              |       |  
                     *       face 2 |       | face 4  
                     *              |_______|  
                     *           (2)        (3) 
                     *               face 3    
                     * */

                    switch(faceID){
                        case 1: {return (Vector2d() << 0,1).finished();}
                        case 2: {return (Vector2d() << 1,2).finished();}
                        case 3: {return (Vector2d() << 2,3).finished();}
                        case 4: {return (Vector2d() << 3,0).finished();}
                        default: {throw std::invalid_argument("Quad4: invalid face ID specifed");}
                    }
                }

                Matrix2d Jacobian(const Ref<const dNdXiSized >& dNdXi, const Ref<const Vector8d>& coordinates)
                {
                    // convenience wrapper to the templated version of the Jacobian
                    return FiniteElement::Jacobian<nDim, nNodes> ( dNdXi, coordinates );
                }

                BSized B(const Ref<const  dNdXiSized>& dNdX)
                {
                    // convenience wrapper to the templated version of the Spatial2D B Operator
                    return FiniteElement::Spatial2D::B<nNodes> ( dNdX ) ;
                }

            } // end of namespace Quad4
            namespace Quad8
            {
                NSized N(const Ref<const Vector2d>& xi){
                    /* Shape functions
                     *
                     *         (6)
                     *   (3) _______(2)
                     *      |       |  
                     *  (7) |       | (5)  
                     *      |_______|  
                     *   (0)        (1) 
                     *         (4)
                     *
                     * */

                    const double xi0 = xi(0);
                    const double xi1 = xi(1);
                    Matrix<double, nNodes, 1> N_;
                    N_ <<                (1.-xi0)*(1.-xi1) * (-1-xi0-xi1)/4,
                       (1.+xi0)*(1.-xi1) * (-1+xi0-xi1)/4,
                       (1.+xi0)*(1.+xi1) * (-1+xi0+xi1)/4,
                       (1.-xi0)*(1.+xi1) * (-1-xi0+xi1)/4,   
                       (1-xi0*xi0)       * (1-xi1      )/2,
                       (1+xi0      )     * (1-xi1*xi1)/2,
                       (1-xi0*xi0)       * (1+xi1     )/2,
                       (1-xi0      )     * (1-xi1*xi1)/2;

                    return N_; 
                }


                dNdXiSized dNdXi(const Ref<const Vector2d>& xi){

                    const double xi0 = xi(0);
                    const double xi1 = xi(1);
                    Matrix<double, nDim, nNodes> result;
                    result <<       
                        /*                  1                           2                           3                           4
                         *                  5                           6                           7                           8 */
                        /* ,Xi1 */     (1-xi1)*(2*xi0+xi1)/4,      (1-xi1)*(2*xi0-xi1)/4,      (1+xi1)*(2*xi0+xi1)/4,      (1+xi1)*(2*xi0-xi1)/4,
                        /* ,Xi1 */     -xi0*(1-xi1),               (1-xi1*xi1)/2,              -xi0*(1+xi1),               -(1-xi1*xi1)/2,
                        /*                  1                           2                           3                           4
                         *                  5                           6                           7                           8*/
                        /* ,Xi2 */     (1-xi0)*(xi0+2*xi1)/4,      -(1+xi0)*(xi0-2*xi1)/4,     (1+xi0)*(xi0+2*xi1)/4,      (1-xi0)*(-xi0+2*xi1)/4,
                        /* ,Xi2 */     -(1-xi0*xi0)/2,             -xi1*(1+xi0),               +(1-xi0*xi0)/2,             -xi1*(1-xi0);


                    return result;
                }

                Matrix2d Jacobian(const Ref<const dNdXiSized >& dNdXi, const Ref<const Matrix<double, nNodes*nDim, 1>>& coordinates)
                {
                    // convenience wrapper to the templated version of the Jacobian
                    return FiniteElement::Jacobian<nDim, nNodes> ( dNdXi, coordinates );
                }

                BSized B(const Ref<const  dNdXiSized>& dNdX)
                {
                    // convenience wrapper to the templated version of the Spatial2D B Operator
                    return FiniteElement::Spatial2D::B<nNodes> ( dNdX ) ;
                }

                Vector3d getBoundaryElementIndices ( int faceID)
                {
                    switch(faceID){
                        case 1: {return (Vector3d() << 0,1,4).finished();}
                        case 2: {return (Vector3d() << 1,2,5).finished();}
                        case 3: {return (Vector3d() << 2,3,6).finished();}
                        case 4: {return (Vector3d() << 3,0,7).finished();}
                        default: {throw std::invalid_argument("Quad8: invalid face ID specifed");}
                    }
                }

            } // end of namespace Quad8

            void modifyCharElemLengthAbaqusLike(double& charElemLength, int intPoint){
                switch(intPoint){
                    case 0:         // central node
                        charElemLength *= 2./3.; break;
                    case 1:         // corner nodes 
                    case 2: 
                    case 3: 
                    case 4:
                        charElemLength *= 5./12; break;
                    case 5:         // middle nodes
                    case 6:
                    case 7:
                    case 8:
                        charElemLength *= std::sqrt(5./18.); break; }
            }
        } // end of namespace Spatial2D 
    } // end of namespace FiniteElement
} // end of namespace bft
