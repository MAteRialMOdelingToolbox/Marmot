#pragma once
#include "bftTypedefs.h"

#define VOIGTSIZE(x) (x + x*x)/2

namespace bft{
    namespace EAS{
        enum EASType{
            DeBorstEAS2,
            SimoRifaiEAS5,
            SimoRifaiEAS4,
        };

        template <int nDim, int nEnhancedStrainParameters>
        Matrix<double, VOIGTSIZE(nDim), nEnhancedStrainParameters>  EASInterpolation( EASType type, const Ref<const Matrix<double, nDim, 1>>& xi)
        {
            Matrix<double, VOIGTSIZE(nDim), nEnhancedStrainParameters> E_;
            switch(type){
                case DeBorstEAS2:
                    { 
                        if(nDim != 2 || nEnhancedStrainParameters != 2)
                            throw std::invalid_argument("Invalid EAS Type Requested - DeBorstEAS2");

                        E_ <<   xi[1],      0,
                                0,          xi[0],
                                0,          0;

                        break;
                    }
                case SimoRifaiEAS5:
                    {
                        if(nDim != 2 || nEnhancedStrainParameters != 5)
                            throw std::invalid_argument("Invalid EAS Type Requested - SimoRifaiEAS5 ");

                        E_ <<   xi[0],      0,      0,      0,      xi[0]*xi[1],
                                0,          xi[1],  0,      0,      -xi[0]*xi[1],
                                0,          0,      xi[0],  xi[1],  xi[0]*xi[0]-xi[1]*xi[1]; 
                        break;
                    }
                case SimoRifaiEAS4:
                    {
                        if(nDim != 2 || nEnhancedStrainParameters != 4)
                            throw std::invalid_argument("Invalid EAS Type Requested - SimoRifaiEAS4 ");

                        E_ <<   xi[0],      0,      0,      0,      
                                0,          xi[1],  0,      0,     
                                0,          0,      xi[0],  xi[1];
                        break;
                    }

                default:    throw std::invalid_argument("Invalid EAS Type Requested"); }

            return E_;
        }

    }
}
