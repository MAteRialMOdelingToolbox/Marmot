#include "MarmotTensor.h"

using namespace Marmot::TensorUtility;

namespace Marmot {
    namespace CommonTensors {
        auto Initialize_I2xI2()
        {
            Tensor3333d I2xI2;

            for ( int i = 0; i < 3; i++ )
                for ( int j = 0; j < 3; j++ )
                    for ( int k = 0; k < 3; k++ )
                        for ( int l = 0; l < 3; l++ ) {
                            I2xI2( i, j, k, l ) = d( i, j ) * d( k, l );
                        }
            return I2xI2;
        }
        const Tensor3333d I2xI2 = Initialize_I2xI2();

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
                for ( int j = 0; j < 3; j++ )
                    for ( int k = 0; k < 3; k++ )
                        for ( int l = 0; l < 3; l++ ) {
                            dsdsigma( i, j, k, l ) = d( i, k ) * d( j, l ) - 1./3 * d(i,j) * d(k,l);
                        }
            return dsdsigma;
        }

        const Tensor3333d dDeviatoricStress_dStress = Initialize_dDeviatoricStress_dStress();

        Tensor333d Initialize_LeviCivita3D()
        {
            Tensor333d e;
            e.setConstant( 0.0 );

            e( 0, 1, 2 ) = 1.0;
            e( 1, 2, 0 ) = 1.0;
            e( 2, 0, 1 ) = 1.0;

            e( 2, 1, 0 ) = -1.0;
            e( 0, 2, 1 ) = -1.0;
            e( 1, 0, 2 ) = -1.0;

            return e;
        }

        Tensor122d Initialize_LeviCivita2D()
        {
            Tensor122d e;
            e.setConstant( 0.0 );

            e( 0, 0, 1 ) = 1.0;
            e( 0, 1, 0 ) = -1.0;

            return e;
        }

        const Tensor333d LeviCivita3D = Initialize_LeviCivita3D();
        const Tensor122d LeviCivita2D = Initialize_LeviCivita2D();

    } // namespace CommonTensors
} // namespace Marmot
