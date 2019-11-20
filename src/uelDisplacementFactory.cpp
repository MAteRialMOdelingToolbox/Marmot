#include "bftFiniteElementUelSpatialWrapper.h"
#include "uelDisplacement.h"
#include "userLibrary.h"

class UelDisplacementRegistrationManager {
    UelDisplacementRegistrationManager() = delete;

  private:
    static bool UelT2D2_isRegistered;
    static bool UelCPS4_isRegistered;
    static bool UelCPS8_isRegistered;
    static bool UelCPE4_isRegistered;
    static bool UelCPE8_isRegistered;
    static bool UelCPS4R_isRegistered;
    static bool UelCPS8R_isRegistered;
    static bool UelCPE4R_isRegistered;
    static bool UelCPE8R_isRegistered;
    static bool UelC3D8_isRegistered;
    static bool UelC3D8R_isRegistered;
    static bool UelC3D20_isRegistered;
    static bool UelC3D20R_isRegistered;
};

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
bool UelDisplacementRegistrationManager::UelT2D2_isRegistered = userLibrary::BftElementFactory::
    registerElement( "UELT2D2", userLibrary::ElementCode::UelT2D2, generateUelT2D2 );

bool UelDisplacementRegistrationManager::UelCPS4_isRegistered = userLibrary::BftElementFactory::
    registerElement( "UELCPS4", userLibrary::ElementCode::UelCPS4, []( int elementID ) -> BftElement* {
        return new UelDisplacement<2, 4>( elementID,
                                          bft::NumIntegration::IntegrationTypes::FullIntegration,
                                          UelDisplacement<2, 4>::SectionType::PlaneStress );
    } );

BftElement* generateUelCPS8( int elementID )
{
    return new UelDisplacement<2, 8>( elementID,
                                      bft::NumIntegration::IntegrationTypes::FullIntegration,
                                      UelDisplacement<2, 8>::SectionType::PlaneStress );
}
bool UelDisplacementRegistrationManager::UelCPS8_isRegistered = userLibrary::BftElementFactory::
    registerElement( "UELCPS8", userLibrary::ElementCode::UelCPS8, generateUelCPS8 );
BftElement* generateUelCPE4( int elementID )
{
    return new UelDisplacement<2, 4>( elementID,
                                      bft::NumIntegration::IntegrationTypes::FullIntegration,
                                      UelDisplacement<2, 4>::SectionType::PlaneStrain );
}
bool UelDisplacementRegistrationManager::UelCPE4_isRegistered = userLibrary::BftElementFactory::
    registerElement( "UELCPE4", userLibrary::ElementCode::UelCPE4, generateUelCPE4 );
BftElement* generateUelCPE8( int elementID )
{
    return new UelDisplacement<2, 8>( elementID,
                                      bft::NumIntegration::IntegrationTypes::FullIntegration,
                                      UelDisplacement<2, 8>::SectionType::PlaneStrain );
}
bool UelDisplacementRegistrationManager::UelCPE8_isRegistered = userLibrary::BftElementFactory::
    registerElement( "UELCPE8", userLibrary::ElementCode::UelCPE8, generateUelCPE8 );
BftElement* generateUelCPS8R( int elementID )
{
    return new UelDisplacement<2, 8>( elementID,
                                      bft::NumIntegration::IntegrationTypes::ReducedIntegration,
                                      UelDisplacement<2, 8>::SectionType::PlaneStress );
}
bool UelDisplacementRegistrationManager::UelCPS8R_isRegistered = userLibrary::BftElementFactory::
    registerElement( "UELCPS8R", userLibrary::ElementCode::UelCPS8R, generateUelCPS8R );
// BftElement* generateUelCPE4R( int elementID )
//{
// return new UelDisplacement<2, 4>( elementID,
// bft::NumIntegration::IntegrationTypes::ReducedIntegration,
// UelDisplacement<2, 4>::SectionType::PlaneStrain );
//}
// bool UelDisplacementRegistrationManager::UelCPE4R_isRegistered =
// userLibrary::BftElementFactory::registerElement("UELCPE4R", userLibrary::ElementCode::UelCPE4R, generateUelCPE4R);
BftElement* generateUelCPE8R( int elementID )
{
    return new UelDisplacement<2, 8>( elementID,
                                      bft::NumIntegration::IntegrationTypes::ReducedIntegration,
                                      UelDisplacement<2, 8>::SectionType::PlaneStrain );
}
bool UelDisplacementRegistrationManager::UelCPE8R_isRegistered = userLibrary::BftElementFactory::
    registerElement( "UELCPE8R", userLibrary::ElementCode::UelCPE8R, generateUelCPE8R );
BftElement* generateUelC3D8( int elementID )
{
    return new UelDisplacement<3, 8>( elementID,
                                      bft::NumIntegration::IntegrationTypes::FullIntegration,
                                      UelDisplacement<3, 8>::SectionType::Solid );
}
bool UelDisplacementRegistrationManager::UelC3D8_isRegistered = userLibrary::BftElementFactory::
    registerElement( "UELC3D8", userLibrary::ElementCode::UelC3D8, generateUelC3D8 );
BftElement* generateUelC3D8R( int elementID )
{
    return new UelDisplacement<3, 8>( elementID,
                                      bft::NumIntegration::IntegrationTypes::ReducedIntegration,
                                      UelDisplacement<3, 8>::SectionType::Solid );
}
bool UelDisplacementRegistrationManager::UelC3D8R_isRegistered = userLibrary::BftElementFactory::
    registerElement( "UELC3D8R", userLibrary::ElementCode::UelC3D8R, generateUelC3D8R );
BftElement* generateUelC3D20( int elementID )
{
    return new UelDisplacement<3, 20>( elementID,
                                       bft::NumIntegration::IntegrationTypes::FullIntegration,
                                       UelDisplacement<3, 20>::SectionType::Solid );
}
bool UelDisplacementRegistrationManager::UelC3D20_isRegistered = userLibrary::BftElementFactory::
    registerElement( "UELC3D20", userLibrary::ElementCode::UelC3D20, generateUelC3D20 );
BftElement* generateUelC3D20R( int elementID )
{
    return new UelDisplacement<3, 20>( elementID,
                                       bft::NumIntegration::IntegrationTypes::ReducedIntegration,
                                       UelDisplacement<3, 20>::SectionType::Solid );
}
bool UelDisplacementRegistrationManager::UelC3D20R_isRegistered = userLibrary::BftElementFactory::
    registerElement( "UELC3D20R", userLibrary::ElementCode::UelC3D20R, generateUelC3D20R );
