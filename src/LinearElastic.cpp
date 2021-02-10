#include "Marmot/LinearElastic.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotUtility.h"
#include "Marmot/MarmotVoigt.h"

namespace Marmot::Materials {

    using namespace Marmot;
    using namespace Eigen;
    using namespace ContinuumMechanics::Elasticity;

    LinearElastic::LinearElastic( const double* materialProperties, int nMaterialProperties, int materialNumber )
        : MarmotMaterialHypoElastic::MarmotMaterialHypoElastic( materialProperties,
                                                                nMaterialProperties,
                                                                materialNumber ),
          // clang-format off
          anisotropicType( static_cast< Type >( nMaterialProperties ) ),
          E1(   materialProperties[0] ),
          E2(   anisotropicType == Type::Isotropic ? E1 : materialProperties[1] ),
          E3(   anisotropicType == Type::Orthotropic ? materialProperties[2] : E2 ),
          nu12( anisotropicType == Type::Isotropic ?           materialProperties[1] :
	            ( anisotropicType == Type::TransverseIsotropic ? materialProperties[2] : 
                                                               materialProperties[3] ) ),
          nu23( anisotropicType == Type::Isotropic ?           nu12 : 
              ( anisotropicType == Type::TransverseIsotropic ? materialProperties[3] : 
                                                               materialProperties[4] ) ),
          nu13( anisotropicType == Type::Orthotropic ? materialProperties[5] : nu23 ),
          G12(  anisotropicType == Type::Isotropic ?           E1 / ( 2 * ( 1 + nu12 ) ):
             (  anisotropicType == Type::TransverseIsotropic ? materialProperties[4] : 
                                                               materialProperties[6] ) ),
          G23(  anisotropicType == Type::Isotropic ?           G12 :
             (  anisotropicType == Type::TransverseIsotropic ? E2 / ( 2 * ( 1 + nu23 ) ) : 
                                                               materialProperties[7] ) ),
          G13(  anisotropicType == Type::Orthotropic ? materialProperties[8] : G12 )
    // clang-format on
    {
        assert( nMaterialProperties == 2 || nMaterialProperties == 8 || nMaterialProperties == 12 );
        switch ( anisotropicType ) {
        case Type::Isotropic: C = Isotropic::stiffnessTensor( E1, nu12 ); break;

        case Type::TransverseIsotropic:
        case Type::Orthotropic:
            Matrix6d localStiffnessTensor;
            Vector3d normalVector;

            int i = nMaterialProperties - 3;
            normalVector << materialProperties[i++], materialProperties[i++], materialProperties[i++];

            switch ( anisotropicType ) {
            case Type::TransverseIsotropic:
                localStiffnessTensor = TransverseIsotropic::stiffnessTensor( E1, E2, nu12, nu23, G12 );
                break;
            case Type::Orthotropic:
                localStiffnessTensor = Orthotropic::stiffnessTensor( E1, E2, E3, nu12, nu23, nu13, G12, G23, G13 );
                                 
                break;
            };

            Matrix3d localCoordinateSystem = Marmot::Math::orthonormalCoordinateSystem( normalVector );

            // strain and stress transformation matrices
            using namespace ContinuumMechanics::VoigtNotation::Transformations;
            Matrix6d transformationStrainInv = transformationMatrixStrainVoigt( localCoordinateSystem ).inverse();
            Matrix6d transformationStress    = transformationMatrixStressVoigt( localCoordinateSystem );

            // transformation Cel into global coordinate system
            C = transformationStrainInv * localStiffnessTensor * transformationStress;
            break;
        };
    }

    void LinearElastic::computeStress( double*       stress,
                                       double*       dStressDDStrain,
                                       const double* dStrain,
                                       const double* timeOld,
                                       const double  dT,
                                       double&       pNewDT )
    {
        mVector6d             S( stress );
        Map< const Vector6d > dE( dStrain );
        mMatrix6d             mC( dStressDDStrain );
        mC = C;
        
        // Zero strain  increment check
        if ( ( dE.array() == 0 ).all() )
            return;

        S += mC * dE;
    }
} // namespace Marmot::Materials
