#include "Marmot/MarmotElasticity.h"

using namespace Eigen;

namespace Marmot {
  namespace ContinuumMechanics {
    namespace Elasticity::Isotropic {

      Matrix6d stiffnessTensor( const double E, const double nu )
      {
        Matrix6d C;
        // clang-format off
            	C <<  (1-nu),     nu,     nu,          0,          0,          0,
            	          nu, (1-nu),     nu,          0,          0,          0,
            	          nu,     nu, (1-nu),          0,          0,          0,
            	          0,       0,      0, (1-2*nu)/2,          0,          0,
            	          0,       0,      0,          0, (1-2*nu)/2,          0,
            	          0,       0,      0,          0,          0, (1-2*nu)/2;
        // clang-format on
        C *= E / ( ( 1 + nu ) * ( 1 - 2 * nu ) );
        return C;
      }

      Matrix6d stiffnessTensorKG( const double K, const double G )
      {
        const double E_  = E( K, G );
        const double nu_ = nu( K, G );
        return stiffnessTensor( E_, nu_ );
      }

      Matrix6d complianceTensor( const double E, const double nu )
      {
        Matrix6d     CInv;
        const double G = E / ( 2 * ( 1 + nu ) );
        // clang-format off
            	  CInv <<   1./E,  -nu/E, -nu/E,     0,    0,    0,
            	           -nu/E,   1./E, -nu/E,     0,    0,    0,
            	           -nu/E,  -nu/E,  1./E,     0,    0,    0,
            	               0,      0,     0,  1./G,    0,    0,
            	               0,      0,     0,     0, 1./G,    0,
            	               0,      0,     0,     0,    0, 1./G;
        // clang-format on
        return CInv;
      }
    } // namespace Elasticity::Isotropic

    namespace Elasticity::TransverseIsotropic {

      Matrix6d complianceTensor( const double E1,
                                 const double E2,
                                 const double nu12,
                                 const double nu23,
                                 const double G12 )
      {
        Matrix6d     CInv;
        const double G23 = E2 / ( 2 * ( 1 + nu23 ) );
        // clang-format off
            	  CInv <<   1./E1, -nu12/E2, -nu12/E2,      0,      0,      0,
                       -nu12/E2,    1./E2, -nu23/E2,      0,      0,      0,
                       -nu12/E2, -nu23/E2,    1./E2,      0,      0,      0,
                              0,        0,        0, 1./G12,      0,      0,
                              0,        0,        0,      0, 1./G12,      0,
                              0,        0,        0,      0,      0, 1./G23;
        // clang-format on
        return CInv;
      }

      Matrix6d stiffnessTensor( const double E1,
                                const double E2,
                                const double nu12,
                                const double nu23,
                                const double G12 )
      {

        Matrix3d CInvTopLeft;

        // clang-format off
        CInvTopLeft <<   1./E1, -nu12/E2, -nu12/E2,
                      -nu12/E2,    1./E2, -nu23/E2,
                      -nu12/E2, -nu23/E2,    1./E2;
        // clang-format on

        Matrix6d     C          = Matrix6d::Zero();
        const double G23        = E2 / ( 2 * ( 1 + nu23 ) );
        C.block< 3, 3 >( 0, 0 ) = CInvTopLeft.inverse();
        C( 3, 3 )               = G12;
        C( 4, 4 )               = G12;
        C( 5, 5 )               = G23;

        return C;
      }

    } // namespace Elasticity::TransverseIsotropic

    namespace Elasticity::Orthotropic {

      Matrix6d complianceTensor( const double E1,
                                 const double E2,
                                 const double E3,
                                 const double nu12,
                                 const double nu23,
                                 const double nu13,
                                 const double G12,
                                 const double G23,
                                 const double G31 )
      {
        Matrix6d CInv;
        // clang-format off
           	    CInv <<   1./E1, -nu12/E2, -nu13/E3,      0,      0,      0,
           	           -nu12/E2,    1./E2, -nu23/E3,      0,      0,      0,
           	           -nu13/E3, -nu23/E3,    1./E3,      0,      0,      0,
           	                  0,        0,        0, 1./G12,      0,      0,
           	              	  0,        0,        0,      0, 1./G31,      0,
           	              	  0,        0,        0,      0,      0, 1./G23;
        // clang-format on
        return CInv;
      }

      Matrix6d stiffnessTensor( const double E1,
                                const double E2,
                                const double E3,
                                const double nu12,
                                const double nu23,
                                const double nu13,
                                const double G12,
                                const double G23,
                                const double G31 )
      {

        Matrix3d CInvTopLeft;

        // clang-format off
        CInvTopLeft <<    1./E1, -nu12/E2, -nu13/E3,
           	           -nu12/E2,    1./E2, -nu23/E3,
           	           -nu13/E3, -nu23/E3,    1./E3;
        // clang-format on

        Matrix6d C              = Matrix6d::Zero();
        C.block< 3, 3 >( 0, 0 ) = CInvTopLeft.inverse();
        C( 3, 3 )               = G12;
        C( 4, 4 )               = G31;
        C( 5, 5 )               = G23;

        return C;
      }
    } // namespace Elasticity::Orthotropic
  }   // namespace ContinuumMechanics
} // namespace Marmot
