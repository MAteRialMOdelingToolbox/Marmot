#include "Marmot/GradientEnhancedFiniteStrainParticleSQCNI.h"
#include "Marmot/MarmotMaterialPoint.h"
#include "Marmot/MarmotParticleLibrary.h"

namespace Marmot::Meshfree {

  using namespace MarmotLibrary;

  const static bool GradientEnhancedFiniteStrainParticleSQCNI_PlaneStrain_Quad_isRegistered = MarmotLibrary::
    MarmotParticleFactory::registerParticle( "GradientEnhancedFiniteStrainSQCNI/PlaneStrain/Quad",
                                             []( int                cellID,
                                                 const double*      nodeCoordinates,
                                                 int                sizeNodeCoordinates,
                                                 double             volume,
                                                 const std::string& materialName,
                                                 const double*      materialProperties,
                                                 int                sizeMaterialProperties,
                                                 const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
                                               -> Marmot::Meshfree::MarmotParticle* {
                                               return new GradientEnhancedFiniteStrainParticleSQCNI<
                                                 2,
                                                 4 >( cellID,
                                                      nodeCoordinates,
                                                      sizeNodeCoordinates,
                                                      volume,
                                                      materialName,
                                                      materialProperties,
                                                      sizeMaterialProperties,
                                                      approximation,
                                                      GradientEnhancedFiniteStrainParticleSQCNI< 2, 4 >::
                                                        SmoothingDomainUpdateType::DeformationGradient );
                                             } );

  const static bool
    GradientEnhancedFiniteStrainParticleSNNI_PlaneStrain_Quad_isRegistered = MarmotLibrary::MarmotParticleFactory::
      registerParticle( "GradientEnhancedFiniteStrainSNNI/PlaneStrain/Quad",
                        []( int                                                  cellID,
                            const double*                                        nodeCoordinates,
                            int                                                  sizeNodeCoordinates,
                            double                                               volume,
                            const std::string&                                   materialName,
                            const double*                                        materialProperties,
                            int                                                  sizeMaterialProperties,
                            const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
                          -> Marmot::Meshfree::MarmotParticle* {
                          return new GradientEnhancedFiniteStrainParticleSQCNI<
                            2,
                            4 >( cellID,
                                 nodeCoordinates,
                                 sizeNodeCoordinates,
                                 volume,
                                 materialName,
                                 materialProperties,
                                 sizeMaterialProperties,
                                 approximation,
                                 GradientEnhancedFiniteStrainParticleSQCNI< 2, 4 >::SmoothingDomainUpdateType::None );
                        } );

  const static bool GradientEnhancedFiniteStrainParticleSQCNI_R_PlaneStrain_Quad_isRegistered = MarmotLibrary::
    MarmotParticleFactory::registerParticle( "GradientEnhancedFiniteStrainSQCNI_R/PlaneStrain/Quad",
                                             []( int                cellID,
                                                 const double*      nodeCoordinates,
                                                 int                sizeNodeCoordinates,
                                                 double             volume,
                                                 const std::string& materialName,
                                                 const double*      materialProperties,
                                                 int                sizeMaterialProperties,
                                                 const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
                                               -> Marmot::Meshfree::MarmotParticle* {
                                               return new GradientEnhancedFiniteStrainParticleSQCNI<
                                                 2,
                                                 4 >( cellID,
                                                      nodeCoordinates,
                                                      sizeNodeCoordinates,
                                                      volume,
                                                      materialName,
                                                      materialProperties,
                                                      sizeMaterialProperties,
                                                      approximation,
                                                      GradientEnhancedFiniteStrainParticleSQCNI< 2, 4 >::
                                                        SmoothingDomainUpdateType::RotationOnly );
                                             } );

  const static bool GradientEnhancedFiniteStrainParticleSQCNI_RU_PlaneStrain_Quad_isRegistered = MarmotLibrary::
    MarmotParticleFactory::registerParticle( "GradientEnhancedFiniteStrainSQCNI_RU/PlaneStrain/Quad",
                                             []( int                cellID,
                                                 const double*      nodeCoordinates,
                                                 int                sizeNodeCoordinates,
                                                 double             volume,
                                                 const std::string& materialName,
                                                 const double*      materialProperties,
                                                 int                sizeMaterialProperties,
                                                 const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
                                               -> Marmot::Meshfree::MarmotParticle* {
                                               return new GradientEnhancedFiniteStrainParticleSQCNI<
                                                 2,
                                                 4 >( cellID,
                                                      nodeCoordinates,
                                                      sizeNodeCoordinates,
                                                      volume,
                                                      materialName,
                                                      materialProperties,
                                                      sizeMaterialProperties,
                                                      approximation,
                                                      GradientEnhancedFiniteStrainParticleSQCNI< 2, 4 >::
                                                        SmoothingDomainUpdateType::RotationAndPrincipalStretch );
                                             } );


  const static bool GradientEnhancedFiniteStrainParticleSQCNI_3D_Hexa_isRegistered = MarmotLibrary::
    MarmotParticleFactory::registerParticle( "GradientEnhancedFiniteStrainSQCNI/3D/Hexa",
                                             []( int                cellID,
                                                 const double*      nodeCoordinates,
                                                 int                sizeNodeCoordinates,
                                                 double             volume,
                                                 const std::string& materialName,
                                                 const double*      materialProperties,
                                                 int                sizeMaterialProperties,
                                                 const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
                                               -> Marmot::Meshfree::MarmotParticle* {
                                               return new GradientEnhancedFiniteStrainParticleSQCNI<
                                                 3,
                                                 8 >( cellID,
                                                      nodeCoordinates,
                                                      sizeNodeCoordinates,
                                                      volume,
                                                      materialName,
                                                      materialProperties,
                                                      sizeMaterialProperties,
                                                      approximation,
                                                      GradientEnhancedFiniteStrainParticleSQCNI< 3, 8 >::
                                                        SmoothingDomainUpdateType::DeformationGradient );
                                             } );

  const static bool GradientEnhancedFiniteStrainParticleSNNI_3D_Hexa_isRegistered = MarmotLibrary::
    MarmotParticleFactory::registerParticle( "GradientEnhancedFiniteStrainSNNI/3D/Hexa",
                                             []( int                cellID,
                                                 const double*      nodeCoordinates,
                                                 int                sizeNodeCoordinates,
                                                 double             volume,
                                                 const std::string& materialName,
                                                 const double*      materialProperties,
                                                 int                sizeMaterialProperties,
                                                 const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
                                               -> Marmot::Meshfree::MarmotParticle* {
                                               return new GradientEnhancedFiniteStrainParticleSQCNI<
                                                 3,
                                                 8 >( cellID,
                                                      nodeCoordinates,
                                                      sizeNodeCoordinates,
                                                      volume,
                                                      materialName,
                                                      materialProperties,
                                                      sizeMaterialProperties,
                                                      approximation,
                                                      GradientEnhancedFiniteStrainParticleSQCNI< 3, 8 >::
                                                        SmoothingDomainUpdateType::None );
                                             } );

  const static bool GradientEnhancedFiniteStrainParticleSQCNI_R_3D_Hexa_isRegistered = MarmotLibrary::
    MarmotParticleFactory::registerParticle( "GradientEnhancedFiniteStrainSQCNI_R/3D/Hexa",
                                             []( int                cellID,
                                                 const double*      nodeCoordinates,
                                                 int                sizeNodeCoordinates,
                                                 double             volume,
                                                 const std::string& materialName,
                                                 const double*      materialProperties,
                                                 int                sizeMaterialProperties,
                                                 const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
                                               -> Marmot::Meshfree::MarmotParticle* {
                                               return new GradientEnhancedFiniteStrainParticleSQCNI<
                                                 3,
                                                 8 >( cellID,
                                                      nodeCoordinates,
                                                      sizeNodeCoordinates,
                                                      volume,
                                                      materialName,
                                                      materialProperties,
                                                      sizeMaterialProperties,
                                                      approximation,
                                                      GradientEnhancedFiniteStrainParticleSQCNI< 3, 8 >::
                                                        SmoothingDomainUpdateType::RotationOnly );
                                             } );

  const static bool GradientEnhancedFiniteStrainParticleSQCNI_RU_3D_Hexa_isRegistered = MarmotLibrary::
    MarmotParticleFactory::registerParticle( "GradientEnhancedFiniteStrainSQCNI_RU/3D/Hexa",
                                             []( int                cellID,
                                                 const double*      nodeCoordinates,
                                                 int                sizeNodeCoordinates,
                                                 double             volume,
                                                 const std::string& materialName,
                                                 const double*      materialProperties,
                                                 int                sizeMaterialProperties,
                                                 const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
                                               -> Marmot::Meshfree::MarmotParticle* {
                                               return new GradientEnhancedFiniteStrainParticleSQCNI<
                                                 3,
                                                 8 >( cellID,
                                                      nodeCoordinates,
                                                      sizeNodeCoordinates,
                                                      volume,
                                                      materialName,
                                                      materialProperties,
                                                      sizeMaterialProperties,
                                                      approximation,
                                                      GradientEnhancedFiniteStrainParticleSQCNI< 3, 8 >::
                                                        SmoothingDomainUpdateType::RotationAndPrincipalStretch );
                                             } );

} // namespace Marmot::Meshfree
