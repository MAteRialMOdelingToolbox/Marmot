#include "uelDisplacementFactory.h"
#include "bftFiniteElementUelSpatialWrapper.h"
#include "uelDisplacement.h"

namespace UelDisplacementFactory {
    BftElement* generateUelT2D2( int noEl )
    {
        auto uelT2D2 = std::unique_ptr<BftElement>(
            new UelDisplacement<1, 2>( noEl,
                                       bft::NumIntegration::IntegrationTypes::FullIntegration,
                                       UelDisplacement<1, 2>::SectionType::UniaxialStress ) );
        constexpr int indicesToBeWrapped[] = {0, 1};
        constexpr int nIndicesToBeWrapped  = 2;
        return new BftElementSpatialWrapper( 2, 1, 2, 2, indicesToBeWrapped, nIndicesToBeWrapped, std::move( uelT2D2 ) );
    }
    BftElement* generateUelCPS4( int noEl )
    {
        return new UelDisplacement<2, 4>( noEl,
                                          bft::NumIntegration::IntegrationTypes::FullIntegration,
                                          UelDisplacement<2, 4>::SectionType::PlaneStress );
    }
    BftElement* generateUelCPS8( int noEl )
    {
        return new UelDisplacement<2, 8>( noEl,
                                          bft::NumIntegration::IntegrationTypes::FullIntegration,
                                          UelDisplacement<2, 8>::SectionType::PlaneStress );
    }
    BftElement* generateUelCPE4( int noEl )
    {
        return new UelDisplacement<2, 4>( noEl,
                                          bft::NumIntegration::IntegrationTypes::FullIntegration,
                                          UelDisplacement<2, 4>::SectionType::PlaneStrain );
    }
    BftElement* generateUelCPE8( int noEl )
    {
        return new UelDisplacement<2, 8>( noEl,
                                          bft::NumIntegration::IntegrationTypes::FullIntegration,
                                          UelDisplacement<2, 8>::SectionType::PlaneStrain );
    }
    BftElement* generateUelCPS4R( int noEl )
    {
        return new UelDisplacement<2, 4>( noEl,
                                          bft::NumIntegration::IntegrationTypes::ReducedIntegration,
                                          UelDisplacement<2, 4>::SectionType::PlaneStress );
    }
    BftElement* generateUelCPS8R( int noEl )
    {
        return new UelDisplacement<2, 8>( noEl,
                                          bft::NumIntegration::IntegrationTypes::ReducedIntegration,
                                          UelDisplacement<2, 8>::SectionType::PlaneStress );
    }
    BftElement* generateUelCPE4R( int noEl )
    {
        return new UelDisplacement<2, 4>( noEl,
                                          bft::NumIntegration::IntegrationTypes::ReducedIntegration,
                                          UelDisplacement<2, 4>::SectionType::PlaneStrain );
    }
    BftElement* generateUelCPE8R( int noEl )
    {
        return new UelDisplacement<2, 8>( noEl,
                                          bft::NumIntegration::IntegrationTypes::ReducedIntegration,
                                          UelDisplacement<2, 8>::SectionType::PlaneStrain );
    }
    BftElement* generateUelC3D8( int noEl )
    {
        return new UelDisplacement<3, 8>( noEl,
                                          bft::NumIntegration::IntegrationTypes::FullIntegration,
                                          UelDisplacement<3, 8>::SectionType::Solid );
    }
    BftElement* generateUelC3D8R( int noEl )
    {
        return new UelDisplacement<3, 8>( noEl,
                                          bft::NumIntegration::IntegrationTypes::ReducedIntegration,
                                          UelDisplacement<3, 8>::SectionType::Solid );
    }
    BftElement* generateUelC3D20( int noEl )
    {
        return new UelDisplacement<3, 20>( noEl,
                                           bft::NumIntegration::IntegrationTypes::FullIntegration,
                                           UelDisplacement<3, 20>::SectionType::Solid );
    }
    BftElement* generateUelC3D20R( int noEl )
    {
        return new UelDisplacement<3, 20>( noEl,
                                           bft::NumIntegration::IntegrationTypes::ReducedIntegration,
                                           UelDisplacement<3, 20>::SectionType::Solid );
    }
} // namespace UelDisplacementFactory
