#include "bftEnhancedAssumedStrain.h"


namespace bft{
    namespace EAS{

        MatrixXd F ( const Ref< const MatrixXd >& J )
        {
            // nDim == 2
            if ( J.cols() == 2 ){

                Matrix3d  F;
                // consistent with SimoRifai
                F <<   J(0,0)*J(0,0),    J(1,0)*J(0,1),        2*J(0,0)*J(0,1),
                        J(0,1)*J(1,0),    J(1,1)*J(1,1),        2*J(1,0)*J(1,1),
                        J(0,0)*J(1,0),    J(0,1)*J(1,1),        J(0,0)*J(1,1)+J(0,1)*J(1,0);
                return F;
            }
            else{
                throw std::invalid_argument ("Invalid Dimension for bft::EnhancedAssumedStrain!" );}
        }

        MatrixXd EASInterpolation ( EASType type, const Ref< const Vector2d >& xi )
        {
            // Implementation for 2D
            
            switch(type){
                case DeBorstEAS2: { 

                       Matrix<double, 3, 2> E_; 

                        E_ <<   xi[1],      0,
                                0,          xi[0],
                                0,          0;

                        return E_;
                    }
                case SimoRifaiEAS5: {

                       Matrix<double, 3, 5> E_; 

                        E_ <<   xi[0],      0,      0,      0,      xi[0]*xi[1],
                                0,          xi[1],  0,      0,      -xi[0]*xi[1],
                                0,          0,      xi[0],  xi[1],  xi[0]*xi[0]-xi[1]*xi[1]; 

                        return E_;
                    }
                case SimoRifaiEAS4: {

                       Matrix<double, 3, 4> E_; 

                        E_ <<   xi[0],      0,      0,      0,      
                                0,          xi[1],  0,      0,     
                                0,          0,      xi[0],  xi[1];
                        return E_;
                    }

                default:    throw std::invalid_argument("Invalid EAS Type Requested"); }

        }
    }
}
