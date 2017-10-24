#pragma once
#include "bftTypedefs.h"

namespace bft{
    namespace EAS{

        enum EASType{
            DeBorstEAS2,
            DeBorstEAS2_P2,
            EAS3,
            DeBorstEAS9,
            SimoRifaiEAS5,
            SimoRifaiEAS4,
        };
        
        MatrixXd F ( const Ref< const MatrixXd >& J );
        
        MatrixXd EASInterpolation ( EASType type, const Ref< const VectorXd >& xi ); 

    }
}
