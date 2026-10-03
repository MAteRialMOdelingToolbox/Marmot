#include "Marmot/DisplacementParticleSQCNIxNSNI.h"
#include "Marmot/MarmotParticleLibrary.h"

namespace Marmot::Meshfree {

  using namespace MarmotLibrary;

  namespace {

    template < int nDim, int nVertices >
    using Particle = DisplacementParticleSQCNIxNSNI< nDim, nVertices >;

    template < int nDim, int nVertices, typename Particle< nDim, nVertices >::SmoothingDomainUpdateType updateType >
    MarmotParticle* create( int                                cellID,
                            const double*                      nodeCoordinates,
                            int                                sizeNodeCoordinates,
                            double                             volume,
                            const std::string&                 materialName,
                            const double*                      materialProperties,
                            int                                sizeMaterialProperties,
                            const MarmotMeshfreeApproximation& approximation )
    {
      return new Particle< nDim, nVertices >( cellID,
                                              nodeCoordinates,
                                              sizeNodeCoordinates,
                                              volume,
                                              materialName,
                                              materialProperties,
                                              sizeMaterialProperties,
                                              approximation,
                                              updateType );
    }

    /// register under the name of the naming convention <Formulation><IntegrationScheme>/<Dimension>/<Shape> (as all
    /// other particles) and under the former name of the displacement NSNI particles, kept for existing input files
    bool registerWithFormerName( const std::string&                             name,
                                 const std::string&                             formerName,
                                 MarmotParticleFactory::particleFactoryFunction factoryFunction )
    {
      return MarmotParticleFactory::registerParticle( name, factoryFunction ) &&
             MarmotParticleFactory::registerParticle( formerName, factoryFunction );
    }

    using U2 = Particle< 2, 4 >::SmoothingDomainUpdateType;
    using U3 = Particle< 3, 8 >::SmoothingDomainUpdateType;

    const static bool registered = registerWithFormerName( "DisplacementSQCNIxNSNI/PlaneStrain/Quad",
                                                           "Displacement/SQCNIxNSNI/PlaneStrain/Quad",
                                                           create< 2, 4, U2::DeformationGradient > ) &&
                                   registerWithFormerName( "DisplacementSQCNIxNSNI/3D/Hexa",
                                                           "Displacement/SQCNIxNSNI/3D/Hexa",
                                                           create< 3, 8, U3::DeformationGradient > ) &&
                                   registerWithFormerName( "DisplacementSNNIxNSNI/PlaneStrain/Quad",
                                                           "Displacement/SNNIxNSNI/PlaneStrain/Quad",
                                                           create< 2, 4, U2::None > ) &&
                                   registerWithFormerName( "DisplacementSNNIxNSNI/3D/Hexa",
                                                           "Displacement/SNNIxNSNI/3D/Hexa",
                                                           create< 3, 8, U3::None > ) &&
                                   registerWithFormerName( "DisplacementSQCNI_RxNSNI/PlaneStrain/Quad",
                                                           "Displacement/R-SNNIxNSNI/PlaneStrain/Quad",
                                                           create< 2, 4, U2::RotationOnly > ) &&
                                   registerWithFormerName( "DisplacementSQCNI_RxNSNI/3D/Hexa",
                                                           "Displacement/R-SNNIxNSNI/3D/Hexa",
                                                           create< 3, 8, U3::RotationOnly > ) &&
                                   registerWithFormerName( "DisplacementSQCNI_RUxNSNI/PlaneStrain/Quad",
                                                           "Displacement/RS-SNNIxNSNI/PlaneStrain/Quad",
                                                           create< 2, 4, U2::RotationAndPrincipalStretch > ) &&
                                   registerWithFormerName( "DisplacementSQCNI_RUxNSNI/3D/Hexa",
                                                           "Displacement/RS-SNNIxNSNI/3D/Hexa",
                                                           create< 3, 8, U3::RotationAndPrincipalStretch > );

  } // namespace

} // namespace Marmot::Meshfree
