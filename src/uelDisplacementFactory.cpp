#include "uelDisplacementFactory.h"
#include "bftFiniteElementUelSpatialWrapper.h"
#include "uelDisplacement.h"

namespace UelDisplacementFactory {
    BftUel* generateUelT2D2( int noEl )
    {
        auto uelT2D2 = std::unique_ptr<BftUel>(
            new UelDisplacement<1, 2>( noEl,
                                       bft::NumIntegration::IntegrationTypes::FullIntegration,
                                       UelDisplacement<1, 2>::SectionType::UniaxialStress ) );
        constexpr int indicesToBeWrapped[] = {0, 1};
        constexpr int nIndicesToBeWrapped  = 2;
        return new BftUelSpatialWrapper( 2, 1, 2, 2, indicesToBeWrapped, nIndicesToBeWrapped, std::move( uelT2D2 ) );
    }
    BftUel* generateUelCPS4( int noEl )
    {
        return new UelDisplacement<2, 4>( noEl,
                                          bft::NumIntegration::IntegrationTypes::FullIntegration,
                                          UelDisplacement<2, 4>::SectionType::PlaneStress );
    }
    BftUel* generateUelCPS8( int noEl )
    {
        return new UelDisplacement<2, 8>( noEl,
                                          bft::NumIntegration::IntegrationTypes::FullIntegration,
                                          UelDisplacement<2, 8>::SectionType::PlaneStress );
    }
    BftUel* generateUelCPE4( int noEl )
    {
        return new UelDisplacement<2, 4>( noEl,
                                          bft::NumIntegration::IntegrationTypes::FullIntegration,
                                          UelDisplacement<2, 4>::SectionType::PlaneStrain );
    }
    BftUel* generateUelCPE8( int noEl )
    {
        return new UelDisplacement<2, 8>( noEl,
                                          bft::NumIntegration::IntegrationTypes::FullIntegration,
                                          UelDisplacement<2, 8>::SectionType::PlaneStrain );
    }
    BftUel* generateUelCPS4R( int noEl )
    {
        return new UelDisplacement<2, 4>( noEl,
                                          bft::NumIntegration::IntegrationTypes::ReducedIntegration,
                                          UelDisplacement<2, 4>::SectionType::PlaneStress );
    }
    BftUel* generateUelCPS8R( int noEl )
    {
        return new UelDisplacement<2, 8>( noEl,
                                          bft::NumIntegration::IntegrationTypes::ReducedIntegration,
                                          UelDisplacement<2, 8>::SectionType::PlaneStress );
    }
    BftUel* generateUelCPE4R( int noEl )
    {
        return new UelDisplacement<2, 4>( noEl,
                                          bft::NumIntegration::IntegrationTypes::ReducedIntegration,
                                          UelDisplacement<2, 4>::SectionType::PlaneStrain );
    }
    BftUel* generateUelCPE8R( int noEl )
    {
        return new UelDisplacement<2, 8>( noEl,
                                          bft::NumIntegration::IntegrationTypes::ReducedIntegration,
                                          UelDisplacement<2, 8>::SectionType::PlaneStrain );
    }
    BftUel* generateUelC3D8( int noEl )
    {
        return new UelDisplacement<3, 8>( noEl,
                                          bft::NumIntegration::IntegrationTypes::FullIntegration,
                                          UelDisplacement<3, 8>::SectionType::Solid );
    }
    BftUel* generateUelC3D8R( int noEl )
    {
        return new UelDisplacement<3, 8>( noEl,
                                          bft::NumIntegration::IntegrationTypes::ReducedIntegration,
                                          UelDisplacement<3, 8>::SectionType::Solid );
    }
    BftUel* generateUelC3D20( int noEl )
    {
        return new UelDisplacement<3, 20>( noEl,
                                           bft::NumIntegration::IntegrationTypes::FullIntegration,
                                           UelDisplacement<3, 20>::SectionType::Solid );
    }
    BftUel* generateUelC3D20R( int noEl )
    {
        return new UelDisplacement<3, 20>( noEl,
                                           bft::NumIntegration::IntegrationTypes::ReducedIntegration,
                                           UelDisplacement<3, 20>::SectionType::Solid );
    }
} // namespace UelDisplacementFactory
