#include "Marmot/LinearElastic.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"

namespace Marmot::Materials {

  using namespace Marmot;
  using namespace Eigen;
  using namespace ContinuumMechanics::Elasticity;

  LinearElastic::LinearElastic( const double* materialProperties, int nMaterialProperties, int materialNumber )
    : MarmotMaterialHypoElastic::MarmotMaterialHypoElastic( materialProperties, nMaterialProperties, materialNumber ),
      // clang-format off
          anisotropicType( nMaterialProperties >= 15 ? Type::Orthotropic : (nMaterialProperties >= 11 ? Type::TransverseIsotropic : Type::Isotropic) ),
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
    // set global stiffness tensor
    if ( anisotropicType == Type::Isotropic ) {
      globalStiffnessTensor = Isotropic::stiffnessTensor( E1, nu12 );
    }
    else {
      // set coordinate system for transversly isotropic and orthotropic materials
      const int i          = nMaterialProperties % 2 == 0 ? nMaterialProperties - 7 : nMaterialProperties - 6;
      Vector3d  direction1 = { materialProperties[i], materialProperties[i + 1], materialProperties[i + 2] };
      Vector3d  direction2 = { materialProperties[i + 3], materialProperties[i + 4], materialProperties[i + 5] };

      Matrix3d localCoordinateSystem = Marmot::Math::orthonormalCoordinateSystem( direction1, direction2 );
      using namespace ContinuumMechanics::Transformations;
      Matrix6d transformStrainToLocalSystem  = transformationMatrixStrainVoigt( localCoordinateSystem );
      Matrix6d transformStressToGlobalSystem = transformationMatrixStressVoigt( localCoordinateSystem )
                                                 .colPivHouseholderQr()
                                                 .solve( Matrix6d::Identity() );

      Matrix6d localStiffnessTensor;

      switch ( anisotropicType ) {
      case Type::Isotropic: globalStiffnessTensor = Isotropic::stiffnessTensor( E1, nu12 ); break;
      case Type::TransverseIsotropic:
        localStiffnessTensor  = TransverseIsotropic::stiffnessTensor( E1, E2, nu12, nu23, G12 );
        globalStiffnessTensor = transformStressToGlobalSystem * localStiffnessTensor * transformStrainToLocalSystem;
        break;
      case Type::Orthotropic:
        localStiffnessTensor  = Orthotropic::stiffnessTensor( E1, E2, E3, nu12, nu23, nu13, G12, G23, G13 );
        globalStiffnessTensor = transformStressToGlobalSystem * localStiffnessTensor * transformStrainToLocalSystem;
        break;
      };
    }
  }

  void LinearElastic::computeStress( state3D&        state,
                                     Matrix6d&       dStressDDStrain,
                                     const Vector6d& dStrain,
                                     const timeInfo& timeInfo ) const
  {

    // map stress, strain increment and stiffness tensor
    dStressDDStrain = globalStiffnessTensor;

    // Zero strain increment check
    if ( ( dStrain.array() == 0 ).all() )
      return;

    // Compute stress increment
    state.stress += dStressDDStrain * dStrain;
  }

  double LinearElastic::getDensity( const double* stateVars ) const
  {
    int base_props = static_cast< int >( anisotropicType );
    if ( nMaterialProperties >= base_props + 1 )
      return materialProperties[base_props];
    else {
      throw std::runtime_error(
        std::string( MakeString() << __PRETTY_FUNCTION__ << ": Density not specified for this material." ) );
    }
  }

} // namespace Marmot::Materials
