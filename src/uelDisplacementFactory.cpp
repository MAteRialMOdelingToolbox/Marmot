#include "uelDisplacementFactory.h"
#include "bftFiniteElementUelSpatialWrapper.h"
#include "uelDisplacement.h"

namespace UelDisplacementFactory {
    BftElement* generateUelT2D2( int elementID )
    {
        auto uelT2D2 = std::unique_ptr<BftElement>(
            new UelDisplacement<1, 2>( elementID,
                                       bft::NumIntegration::IntegrationTypes::FullIntegration,
                                       UelDisplacement<1, 2>::SectionType::UniaxialStress ) );
        constexpr static int indicesToBeWrapped[] = {0, 1};
        constexpr static int nIndicesToBeWrapped  = 2;
        return new BftElementSpatialWrapper( 2, 1, 2, 2, indicesToBeWrapped, nIndicesToBeWrapped, std::move( uelT2D2 ) );
    }
    BftElement* generateUelCPS4( int elementID )
    {
        return new UelDisplacement<2, 4>( elementID,
                                          bft::NumIntegration::IntegrationTypes::FullIntegration,
                                          UelDisplacement<2, 4>::SectionType::PlaneStress );
    }
    BftElement* generateUelCPS8( int elementID )
    {
        return new UelDisplacement<2, 8>( elementID,
                                          bft::NumIntegration::IntegrationTypes::FullIntegration,
                                          UelDisplacement<2, 8>::SectionType::PlaneStress );
    }
    BftElement* generateUelCPE4( int elementID )
    {
        return new UelDisplacement<2, 4>( elementID,
                                          bft::NumIntegration::IntegrationTypes::FullIntegration,
                                          UelDisplacement<2, 4>::SectionType::PlaneStrain );
    }
    BftElement* generateUelCPE8( int elementID )
    {
        return new UelDisplacement<2, 8>( elementID,
                                          bft::NumIntegration::IntegrationTypes::FullIntegration,
                                          UelDisplacement<2, 8>::SectionType::PlaneStrain );
    }
    BftElement* generateUelCPS4R( int elementID )
    {
        return new UelDisplacement<2, 4>( elementID,
                                          bft::NumIntegration::IntegrationTypes::ReducedIntegration,
                                          UelDisplacement<2, 4>::SectionType::PlaneStress );
    }
    BftElement* generateUelCPS8R( int elementID )
    {
        return new UelDisplacement<2, 8>( elementID,
                                          bft::NumIntegration::IntegrationTypes::ReducedIntegration,
                                          UelDisplacement<2, 8>::SectionType::PlaneStress );
    }
    BftElement* generateUelCPE4R( int elementID )
    {
        return new UelDisplacement<2, 4>( elementID,
                                          bft::NumIntegration::IntegrationTypes::ReducedIntegration,
                                          UelDisplacement<2, 4>::SectionType::PlaneStrain );
    }
    BftElement* generateUelCPE8R( int elementID )
    {
        return new UelDisplacement<2, 8>( elementID,
                                          bft::NumIntegration::IntegrationTypes::ReducedIntegration,
                                          UelDisplacement<2, 8>::SectionType::PlaneStrain );
    }
    BftElement* generateUelC3D8( int elementID )
    {
        return new UelDisplacement<3, 8>( elementID,
                                          bft::NumIntegration::IntegrationTypes::FullIntegration,
                                          UelDisplacement<3, 8>::SectionType::Solid );
    }
    BftElement* generateUelC3D8R( int elementID )
    {
        return new UelDisplacement<3, 8>( elementID,
                                          bft::NumIntegration::IntegrationTypes::ReducedIntegration,
                                          UelDisplacement<3, 8>::SectionType::Solid );
    }
    BftElement* generateUelC3D20( int elementID )
    {
        return new UelDisplacement<3, 20>( elementID,
                                           bft::NumIntegration::IntegrationTypes::FullIntegration,
                                           UelDisplacement<3, 20>::SectionType::Solid );
    }
    BftElement* generateUelC3D20R( int elementID )
    {
        return new UelDisplacement<3, 20>( elementID,
                                           bft::NumIntegration::IntegrationTypes::ReducedIntegration,
                                           UelDisplacement<3, 20>::SectionType::Solid );
    }
} // namespace UelDisplacementFactory
