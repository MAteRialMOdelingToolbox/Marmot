#include "Marmot/MarmotElasticity.h"
#include <iostream>

using namespace Eigen;

namespace Marmot {
    namespace ContinuumMechanics {
        namespace Elasticity::Isotropic {

            Matrix6d stiffnessTensor( double E, double nu )
            {
                Matrix6d C;
                // clang-format off
            	C <<  (1-nu), nu, nu, 0, 0, 0,
            	    nu, (1-nu), nu, 0, 0, 0,
            	    nu, nu, (1-nu), 0, 0, 0,
            	    0, 0, 0, (1-2*nu)/2, 0, 0,
            	    0, 0, 0, 0, (1-2*nu)/2, 0,
            	    0, 0, 0, 0, 0, (1-2*nu)/2;
                // clang-format on
                C *= E / ( ( 1 + nu ) * ( 1 - 2 * nu ) );
                return C;
            }

            Matrix6d complianceTensor( double E, double nu )
            {
                Matrix6d     CInv;
                const double G = E / ( 2 * ( 1 + nu ) );
                // clang-format off
            	CInv <<   1./E,   -nu/E,  -nu/E,  0,      0,      0,
            	       -nu/E,  1./E,   -nu/E,  0,      0,      0,
            	       -nu/E,  -nu/E,  1./E,   0,      0,      0,
            	       0,      0,      0,      1./G,   0,      0,
            	       0,      0,      0,      0,      1./G,   0,
            	       0,      0,      0,      0,      0,      1./G;
                // clang-format on
                return CInv;
            }
        } // namespace Elasticity::Isotropic

        namespace Elasticity::Anisotropic {

            Matrix6d transverseIsotropicComplianceTensor( double E1, double E2, double nu1, double nu2, double G2 )
            {
                Matrix6d     CInv;
                const double G1 = E1 / ( 2 * ( 1 + nu1 ) );
                // clang-format off
            	CInv <<   1./E1,   -nu1/E1,  -nu2/E2,  0,      0,      0,
                   -nu1/E1,  1./E1,   -nu2/E2,  0,      0,      0,
                   -nu1/E1,  -nu1/E1,  1./E2,   0,      0,      0,
                   0,      0,      0,      1./G1,   0,      0,
                   0,      0,      0,      0,      1./G2,   0,
                   0,      0,      0,      0,      0,      1./G2;
                // clang-format on
                return CInv;
            }

            Matrix6d orthotropicComplianceTensor( double E1,
                                                  double E2,
                                                  double E3,
                                                  double nu12,
                                                  double nu23,
                                                  double nu13,
                                                  double G12,
                                                  double G23,
                                                  double G31 )
            {
                Matrix6d CInv;
                // clang-format off
           	 CInv <<  1./E1, -nu12/E2, -nu13/E3,  0,      0,      0,
           	              -nu12/E2,  1./E2,   -nu23/E3,  0,      0,      0,
           	              -nu13/E3,  -nu23/E3,  1./E3,   0,      0,      0,
           	                   0,      0,      0,      1./G12,   0,      0,
           		       	   0,      0,      0,      0,      1./G23,   0,
           		       	   0,      0,      0,      0,      0,      1./G31;
                // clang-format on
                return CInv;
            }
        } // namespace Elasticity::Anisotropic
    }     // namespace ContinuumMechanics
} // namespace Marmot
