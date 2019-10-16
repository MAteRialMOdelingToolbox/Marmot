#include "bftTensor.h"

using namespace bft::TensorUtility;

namespace bft {
    namespace CommonTensors {
        auto Initialize_I()
        {
            Tensor3333d I;

            for ( int i = 0; i < 3; i++ )
                for ( int j = 0; j < 3; j++ )
                    for ( int k = 0; k < 3; k++ )
                        for ( int l = 0; l < 3; l++ ) {
                            I( i, j, k, k ) = d( i, j ) * d( k, l );
                        }
            return I;
        }
        const Tensor3333d I = Initialize_I();

        auto Initialize_Isym()
        {
            Tensor3333d Isym;

            for ( int i = 0; i < 3; i++ )
                for ( int j = 0; j < 3; j++ )
                    for ( int k = 0; k < 3; k++ )
                        for ( int l = 0; l < 3; l++ ) {
                            Isym(i,j,k,l)      = 0.5 * ( d(i,k)*d(j,l) + d(i,l)*d(j,k)) ;
                        }
            return Isym;
        }
        const Tensor3333d Isym = Initialize_Isym();

        auto Initialize_Iskew()
        {
            Tensor3333d Iskew;

            for ( int i = 0; i < 3; i++ )
                for ( int j = 0; j < 3; j++ )
                    for ( int k = 0; k < 3; k++ )
                        for ( int l = 0; l < 3; l++ ) {
                            Iskew(i,j,k,l)    = 0.5 * ( d(i,k)*d(j,l) - d(i,l)*d(j,k)) ;
                        }
            return Iskew;
        }
        const Tensor3333d Iskew = Initialize_Iskew();

        auto Initialize_dDeviatoricStress_dStress()
        {
            Tensor3333d dsdsigma;
            dsdsigma.setZero();
            for ( int i = 0; i < 3; i++ )
                for ( int j = 0; j < 3; j++ ) {
                    dsdsigma( i, j, i, j ) += 1.0;
                    dsdsigma( i, i, j, j ) -= 1. / 3;
                }
            return dsdsigma;
        }

        const Tensor3333d dDeviatoricStress_dStress = Initialize_dDeviatoricStress_dStress();
    } // namespace CommonTensors
} // namespace bft
