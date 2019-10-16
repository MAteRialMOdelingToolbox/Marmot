#include "bftCosserat.h"
#include "bftTensor.h"

using namespace Eigen;

namespace bft {
    namespace Cosserat {

        bool checkElasticParametersForThermodynamicValidity ( double G, double lameParameter, double gamma, double alpha ){
            return ( 2 * G + 3 * lameParameter < 0 ) &&  ( 2 * gamma + 3 * alpha < 0 );
        }

        double I1 (const Eigen::Matrix3d& stress)
        {
            return stress.trace();
        }

        Eigen::Matrix3d deviatoricStress( const Eigen::Matrix3d& stress )
        {
            return stress - Matrix3d::Identity() * stress.trace()/3;
        }

        double J2( const Eigen::Matrix3d& stress,
                   const Eigen::Matrix3d& coupleStress,
                   double                 a1,
                   double                 a2,
                   double                 a3,
                   double                 a4,
                   double                 l )
        {
            const auto  s = deviatoricStress( stress );
            const auto& m = coupleStress;

            return ( a1 * ( s * s.transpose() ).trace() + a2 * ( s * s ).trace() ) +
                   ( a3 * ( m * m.transpose() ).trace() + a4 * ( m * m ).trace() ) / ( l * l );
        }

        Eigen::Matrix3d dJ2_dStress( const Eigen::Matrix3d& stress, double a1, double a2 )
        {
            using namespace bft::TensorUtility;
            const Matrix3d s      = deviatoricStress( stress );
            const Matrix3d dJ2_ds = s * ( 2 * a1 + a2 ) + s.transpose() * a2;

            const Vector9d dJ2_dStress =  flatten( dJ2_ds ).transpose() * as<9, 9>( bft::CommonTensors::dDeviatoricStress_dStress );

            return as<3,3> ( dJ2_dStress );
        }

        Eigen::Matrix3d dJ2_dCoupleStress( const Eigen::Matrix3d& coupleStress, double a3, double a4, double l )
        {
            const auto& m = coupleStress;
            return 1. / ( l * l ) * ( m * ( 2 * a3 + a4 ) + a4 * m.transpose() ); }

        Tensor3333d Cel( double lameParameter, double G, double Gc )
        {
            using namespace bft::CommonTensors;
            Tensor3333d Cel = I * I.constant( lameParameter ) + Isym * Isym.constant( 2 * G ) +
                              Iskew * Iskew.constant( 2 * Gc );

            return Cel;
        }
        Tensor3333d CelM( double alpha, double beta, double gamma )
        {
            using namespace bft::CommonTensors;
            Tensor3333d Celm = I * I.constant( alpha ) + Isym * Isym.constant( 2 * gamma ) +
                               Iskew * Iskew.constant( 2 * beta );
            return Celm;
        }

    } // namespace Cosserat

} // namespace bft
