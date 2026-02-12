#include "Marmot/DisplacementParticleSQCNIxNSNI.h"
#include "Marmot/MarmotParticleLibrary.h"

namespace Marmot::Meshfree {

  using namespace MarmotLibrary;

  const static bool
    DisplacementParticleSQCNIxNSNI_PlaneStrain_Quad_isRegistered = MarmotLibrary::MarmotParticleFactory::
      registerParticle( "Displacement/SQCNIxNSNI/PlaneStrain/Quad",
                        []( int                                                  cellID,
                            const double*                                        nodeCoordinates,
                            int                                                  sizeNodeCoordinates,
                            double                                               volume,
                            const std::string&                                   materialName,
                            const double*                                        materialProperties,
                            int                                                  sizeMaterialProperties,
                            const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
                          -> Marmot::Meshfree::MarmotParticle* {
                          return new DisplacementParticleSQCNIxNSNI<
                            2,
                            4 >( cellID,
                                 nodeCoordinates,
                                 sizeNodeCoordinates,
                                 volume,
                                 materialName,
                                 materialProperties,
                                 sizeMaterialProperties,
                                 approximation,
                                 DisplacementParticleSQCNIxNSNI< 2,
                                                                 4 >::SmoothingDomainUpdateType::DeformationGradient );
                        } );

  const static bool DisplacementParticleSQCNIxNSNI_3D_Hexa_isRegistered = MarmotLibrary::MarmotParticleFactory::
    registerParticle( "Displacement/SQCNIxNSNI/3D/Hexa",
                      []( int                                                  cellID,
                          const double*                                        nodeCoordinates,
                          int                                                  sizeNodeCoordinates,
                          double                                               volume,
                          const std::string&                                   materialName,
                          const double*                                        materialProperties,
                          int                                                  sizeMaterialProperties,
                          const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
                        -> Marmot::Meshfree::MarmotParticle* {
                        return new DisplacementParticleSQCNIxNSNI<
                          3,
                          8 >( cellID,
                               nodeCoordinates,
                               sizeNodeCoordinates,
                               volume,
                               materialName,
                               materialProperties,
                               sizeMaterialProperties,
                               approximation,
                               DisplacementParticleSQCNIxNSNI< 3, 8 >::SmoothingDomainUpdateType::DeformationGradient );
                      } );

  const static bool DisplacementParticleSNNIxNSNI_PlaneStrain_Quad_isRegistered = MarmotLibrary::MarmotParticleFactory::
    registerParticle( "Displacement/SNNIxNSNI/PlaneStrain/Quad",
                      []( int                                                  cellID,
                          const double*                                        nodeCoordinates,
                          int                                                  sizeNodeCoordinates,
                          double                                               volume,
                          const std::string&                                   materialName,
                          const double*                                        materialProperties,
                          int                                                  sizeMaterialProperties,
                          const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
                        -> Marmot::Meshfree::MarmotParticle* {
                        return new DisplacementParticleSQCNIxNSNI< 2, 4 >( cellID,
                                                                           nodeCoordinates,
                                                                           sizeNodeCoordinates,
                                                                           volume,
                                                                           materialName,
                                                                           materialProperties,
                                                                           sizeMaterialProperties,
                                                                           approximation,
                                                                           DisplacementParticleSQCNIxNSNI< 2, 4 >::
                                                                             SmoothingDomainUpdateType::None );
                      } );

  const static bool DisplacementParticleSNNIxNSNI_3D_Hexa_isRegistered = MarmotLibrary::MarmotParticleFactory::
    registerParticle( "Displacement/SNNIxNSNI/3D/Hexa",
                      []( int                                                  cellID,
                          const double*                                        nodeCoordinates,
                          int                                                  sizeNodeCoordinates,
                          double                                               volume,
                          const std::string&                                   materialName,
                          const double*                                        materialProperties,
                          int                                                  sizeMaterialProperties,
                          const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
                        -> Marmot::Meshfree::MarmotParticle* {
                        return new DisplacementParticleSQCNIxNSNI< 3, 8 >( cellID,
                                                                           nodeCoordinates,
                                                                           sizeNodeCoordinates,
                                                                           volume,
                                                                           materialName,
                                                                           materialProperties,
                                                                           sizeMaterialProperties,
                                                                           approximation,
                                                                           DisplacementParticleSQCNIxNSNI< 3, 8 >::
                                                                             SmoothingDomainUpdateType::None );
                      } );

  const static bool
    DisplacementParticleSQCNIxNSNI_R_PlaneStrain_Quad_isRegistered = MarmotLibrary::MarmotParticleFactory::
      registerParticle( "Displacement/R-SNNIxNSNI/PlaneStrain/Quad",
                        []( int                                                  cellID,
                            const double*                                        nodeCoordinates,
                            int                                                  sizeNodeCoordinates,
                            double                                               volume,
                            const std::string&                                   materialName,
                            const double*                                        materialProperties,
                            int                                                  sizeMaterialProperties,
                            const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
                          -> Marmot::Meshfree::MarmotParticle* {
                          return new DisplacementParticleSQCNIxNSNI<
                            2,
                            4 >( cellID,
                                 nodeCoordinates,
                                 sizeNodeCoordinates,
                                 volume,
                                 materialName,
                                 materialProperties,
                                 sizeMaterialProperties,
                                 approximation,
                                 DisplacementParticleSQCNIxNSNI< 2, 4 >::SmoothingDomainUpdateType::RotationOnly );
                        } );

  const static bool DisplacementParticle_R_SNNIxNSNI_3D_Hexa_isRegistered = MarmotLibrary::MarmotParticleFactory::
    registerParticle( "Displacement/R-SNNIxNSNI/3D/Hexa",
                      []( int                                                  cellID,
                          const double*                                        nodeCoordinates,
                          int                                                  sizeNodeCoordinates,
                          double                                               volume,
                          const std::string&                                   materialName,
                          const double*                                        materialProperties,
                          int                                                  sizeMaterialProperties,
                          const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
                        -> Marmot::Meshfree::MarmotParticle* {
                        return new DisplacementParticleSQCNIxNSNI< 3, 8 >( cellID,
                                                                           nodeCoordinates,
                                                                           sizeNodeCoordinates,
                                                                           volume,
                                                                           materialName,
                                                                           materialProperties,
                                                                           sizeMaterialProperties,
                                                                           approximation,
                                                                           DisplacementParticleSQCNIxNSNI< 3, 8 >::
                                                                             SmoothingDomainUpdateType::RotationOnly );
                      } );

  const static bool DisplacementParticle_RS_SNNIxNSNIPlaneStrain_Quad_isRegistered = MarmotLibrary::
    MarmotParticleFactory::registerParticle( "Displacement/RS-SNNIxNSNI/PlaneStrain/Quad",
                                             []( int                cellID,
                                                 const double*      nodeCoordinates,
                                                 int                sizeNodeCoordinates,
                                                 double             volume,
                                                 const std::string& materialName,
                                                 const double*      materialProperties,
                                                 int                sizeMaterialProperties,
                                                 const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
                                               -> Marmot::Meshfree::MarmotParticle* {
                                               return new DisplacementParticleSQCNIxNSNI<
                                                 2,
                                                 4 >( cellID,
                                                      nodeCoordinates,
                                                      sizeNodeCoordinates,
                                                      volume,
                                                      materialName,
                                                      materialProperties,
                                                      sizeMaterialProperties,
                                                      approximation,
                                                      DisplacementParticleSQCNIxNSNI< 2, 4 >::
                                                        SmoothingDomainUpdateType::RotationAndPrincipalStretch );
                                             } );

  const static bool DisplacementParticle_RS_SNNIxNSNI_3D_Hexa_isRegistered = MarmotLibrary::MarmotParticleFactory::
    registerParticle( "Displacement/RS-SNNIxNSNI/3D/Hexa",
                      []( int                                                  cellID,
                          const double*                                        nodeCoordinates,
                          int                                                  sizeNodeCoordinates,
                          double                                               volume,
                          const std::string&                                   materialName,
                          const double*                                        materialProperties,
                          int                                                  sizeMaterialProperties,
                          const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
                        -> Marmot::Meshfree::MarmotParticle* {
                        return new DisplacementParticleSQCNIxNSNI< 3, 8 >( cellID,
                                                                           nodeCoordinates,
                                                                           sizeNodeCoordinates,
                                                                           volume,
                                                                           materialName,
                                                                           materialProperties,
                                                                           sizeMaterialProperties,
                                                                           approximation,
                                                                           DisplacementParticleSQCNIxNSNI< 3, 8 >::
                                                                             SmoothingDomainUpdateType::
                                                                               RotationAndPrincipalStretch );
                      } );

} // namespace Marmot::Meshfree
