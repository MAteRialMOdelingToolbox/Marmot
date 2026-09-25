#include "Marmot/GradientEnhancedFiniteStrainParticleSQCNIxNSNI.h"
#include "Marmot/MarmotParticleLibrary.h"

namespace Marmot::Meshfree {

  using namespace MarmotLibrary;

  const static bool GradientEnhancedFiniteStrainParticleSQCNIxNSNI_PlaneStrain_Quad_isRegistered = MarmotLibrary::
    MarmotParticleFactory::registerParticle( "GradientEnhancedFiniteStrainSQCNIxNSNI/PlaneStrain/Quad",
                                             []( int                cellID,
                                                 const double*      nodeCoordinates,
                                                 int                sizeNodeCoordinates,
                                                 double             volume,
                                                 const std::string& materialName,
                                                 const double*      materialProperties,
                                                 int                sizeMaterialProperties,
                                                 const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
                                               -> Marmot::Meshfree::MarmotParticle* {
                                               return new GradientEnhancedFiniteStrainParticleSQCNIxNSNI<
                                                 2,
                                                 4 >( cellID,
                                                      nodeCoordinates,
                                                      sizeNodeCoordinates,
                                                      volume,
                                                      materialName,
                                                      materialProperties,
                                                      sizeMaterialProperties,
                                                      approximation,
                                                      GradientEnhancedFiniteStrainParticleSQCNIxNSNI< 2, 4 >::
                                                        SmoothingDomainUpdateType::DeformationGradient );
                                             } );

  const static bool
    GradientEnhancedFiniteStrainParticleSNNIxNSNI_PlaneStrain_Quad_isRegistered = MarmotLibrary::MarmotParticleFactory::
      registerParticle( "GradientEnhancedFiniteStrainSNNIxNSNI/PlaneStrain/Quad",
                        []( int                                                  cellID,
                            const double*                                        nodeCoordinates,
                            int                                                  sizeNodeCoordinates,
                            double                                               volume,
                            const std::string&                                   materialName,
                            const double*                                        materialProperties,
                            int                                                  sizeMaterialProperties,
                            const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
                          -> Marmot::Meshfree::MarmotParticle* {
                          return new GradientEnhancedFiniteStrainParticleSQCNIxNSNI<
                            2,
                            4 >( cellID,
                                 nodeCoordinates,
                                 sizeNodeCoordinates,
                                 volume,
                                 materialName,
                                 materialProperties,
                                 sizeMaterialProperties,
                                 approximation,
                                 GradientEnhancedFiniteStrainParticleSQCNIxNSNI< 2,
                                                                                 4 >::SmoothingDomainUpdateType::None );
                        } );

  const static bool GradientEnhancedFiniteStrainParticleSQCNIxNSNI_R_PlaneStrain_Quad_isRegistered = MarmotLibrary::
    MarmotParticleFactory::registerParticle( "GradientEnhancedFiniteStrainSQCNI_RxNSNI/PlaneStrain/Quad",
                                             []( int                cellID,
                                                 const double*      nodeCoordinates,
                                                 int                sizeNodeCoordinates,
                                                 double             volume,
                                                 const std::string& materialName,
                                                 const double*      materialProperties,
                                                 int                sizeMaterialProperties,
                                                 const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
                                               -> Marmot::Meshfree::MarmotParticle* {
                                               return new GradientEnhancedFiniteStrainParticleSQCNIxNSNI<
                                                 2,
                                                 4 >( cellID,
                                                      nodeCoordinates,
                                                      sizeNodeCoordinates,
                                                      volume,
                                                      materialName,
                                                      materialProperties,
                                                      sizeMaterialProperties,
                                                      approximation,
                                                      GradientEnhancedFiniteStrainParticleSQCNIxNSNI< 2, 4 >::
                                                        SmoothingDomainUpdateType::RotationOnly );
                                             } );

  const static bool GradientEnhancedFiniteStrainParticleSQCNIxNSNI_RU_PlaneStrain_Quad_isRegistered = MarmotLibrary::
    MarmotParticleFactory::registerParticle( "GradientEnhancedFiniteStrainSQCNI_RUxNSNI/PlaneStrain/Quad",
                                             []( int                cellID,
                                                 const double*      nodeCoordinates,
                                                 int                sizeNodeCoordinates,
                                                 double             volume,
                                                 const std::string& materialName,
                                                 const double*      materialProperties,
                                                 int                sizeMaterialProperties,
                                                 const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
                                               -> Marmot::Meshfree::MarmotParticle* {
                                               return new GradientEnhancedFiniteStrainParticleSQCNIxNSNI<
                                                 2,
                                                 4 >( cellID,
                                                      nodeCoordinates,
                                                      sizeNodeCoordinates,
                                                      volume,
                                                      materialName,
                                                      materialProperties,
                                                      sizeMaterialProperties,
                                                      approximation,
                                                      GradientEnhancedFiniteStrainParticleSQCNIxNSNI< 2, 4 >::
                                                        SmoothingDomainUpdateType::RotationAndPrincipalStretch );
                                             } );

  const static bool GradientEnhancedFiniteStrainParticleSQCNIxNSNI_3D_Hexa_isRegistered = MarmotLibrary::
    MarmotParticleFactory::registerParticle( "GradientEnhancedFiniteStrainSQCNIxNSNI/3D/Hexa",
                                             []( int                cellID,
                                                 const double*      nodeCoordinates,
                                                 int                sizeNodeCoordinates,
                                                 double             volume,
                                                 const std::string& materialName,
                                                 const double*      materialProperties,
                                                 int                sizeMaterialProperties,
                                                 const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
                                               -> Marmot::Meshfree::MarmotParticle* {
                                               return new GradientEnhancedFiniteStrainParticleSQCNIxNSNI<
                                                 3,
                                                 8 >( cellID,
                                                      nodeCoordinates,
                                                      sizeNodeCoordinates,
                                                      volume,
                                                      materialName,
                                                      materialProperties,
                                                      sizeMaterialProperties,
                                                      approximation,
                                                      GradientEnhancedFiniteStrainParticleSQCNIxNSNI< 3, 8 >::
                                                        SmoothingDomainUpdateType::DeformationGradient );
                                             } );

  const static bool
    GradientEnhancedFiniteStrainParticleSNNIxNSNI_3D_Hexa_isRegistered = MarmotLibrary::MarmotParticleFactory::
      registerParticle( "GradientEnhancedFiniteStrainSNNIxNSNI/3D/Hexa",
                        []( int                                                  cellID,
                            const double*                                        nodeCoordinates,
                            int                                                  sizeNodeCoordinates,
                            double                                               volume,
                            const std::string&                                   materialName,
                            const double*                                        materialProperties,
                            int                                                  sizeMaterialProperties,
                            const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
                          -> Marmot::Meshfree::MarmotParticle* {
                          return new GradientEnhancedFiniteStrainParticleSQCNIxNSNI<
                            3,
                            8 >( cellID,
                                 nodeCoordinates,
                                 sizeNodeCoordinates,
                                 volume,
                                 materialName,
                                 materialProperties,
                                 sizeMaterialProperties,
                                 approximation,
                                 GradientEnhancedFiniteStrainParticleSQCNIxNSNI< 3,
                                                                                 8 >::SmoothingDomainUpdateType::None );
                        } );

  const static bool GradientEnhancedFiniteStrainParticleSQCNI_RxNSNI_3D_Hexa_isRegistered = MarmotLibrary::
    MarmotParticleFactory::registerParticle( "GradientEnhancedFiniteStrainSQCNI_RxNSNI/3D/Hexa",
                                             []( int                cellID,
                                                 const double*      nodeCoordinates,
                                                 int                sizeNodeCoordinates,
                                                 double             volume,
                                                 const std::string& materialName,
                                                 const double*      materialProperties,
                                                 int                sizeMaterialProperties,
                                                 const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
                                               -> Marmot::Meshfree::MarmotParticle* {
                                               return new GradientEnhancedFiniteStrainParticleSQCNIxNSNI<
                                                 3,
                                                 8 >( cellID,
                                                      nodeCoordinates,
                                                      sizeNodeCoordinates,
                                                      volume,
                                                      materialName,
                                                      materialProperties,
                                                      sizeMaterialProperties,
                                                      approximation,
                                                      GradientEnhancedFiniteStrainParticleSQCNIxNSNI< 3, 8 >::
                                                        SmoothingDomainUpdateType::RotationOnly );
                                             } );

  const static bool GradientEnhancedFiniteStrainParticleSQCNI_RUxNSNI_3D_Hexa_isRegistered = MarmotLibrary::
    MarmotParticleFactory::registerParticle( "GradientEnhancedFiniteStrainSQCNI_RUxNSNI/3D/Hexa",
                                             []( int                cellID,
                                                 const double*      nodeCoordinates,
                                                 int                sizeNodeCoordinates,
                                                 double             volume,
                                                 const std::string& materialName,
                                                 const double*      materialProperties,
                                                 int                sizeMaterialProperties,
                                                 const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
                                               -> Marmot::Meshfree::MarmotParticle* {
                                               return new GradientEnhancedFiniteStrainParticleSQCNIxNSNI<
                                                 3,
                                                 8 >( cellID,
                                                      nodeCoordinates,
                                                      sizeNodeCoordinates,
                                                      volume,
                                                      materialName,
                                                      materialProperties,
                                                      sizeMaterialProperties,
                                                      approximation,
                                                      GradientEnhancedFiniteStrainParticleSQCNIxNSNI< 3, 8 >::
                                                        SmoothingDomainUpdateType::RotationAndPrincipalStretch );
                                             } );

} // namespace Marmot::Meshfree
