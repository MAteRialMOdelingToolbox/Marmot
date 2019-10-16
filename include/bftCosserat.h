#pragma once
#include "bftTensor.h"
#include "bftTypedefs.h"

namespace bft {
    namespace Cosserat {
        bool checkElasticParametersForThermodynamicValidity ( double G, double lameParameter, double gamma, double alpha );

        Eigen::Matrix3d deviatoricStress( const Eigen::Matrix3d& stress );

        double I1 (const Eigen::Matrix3d& stress);
        double J2( const Eigen::Matrix3d& stress, const Eigen::Matrix3d& coupleStress, 
                   double a1, double a2, double a3, double a4, double l);

        Eigen::Matrix3d dJ2_dStress( const Eigen::Matrix3d& stress, double a1, double a2); 
        Eigen::Matrix3d dJ2_dCoupleStress( const Eigen::Matrix3d& coupleStress, double a3, double a4, double l);

        Tensor3333d Cel( double lameParameter, double G, double Gc );
        Tensor3333d CelM( double alpha, double beta, double gamma );

    } // namespace Cosserat

} // namespace bft
