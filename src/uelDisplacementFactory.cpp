#include "bftFiniteElement.h"
#include "bftFiniteElementUelSpatialWrapper.h"
#include "uelDisplacement.h"
#include "userLibrary.h"

namespace UelDisplacementRegistration {

    template <class T, bft::NumIntegration::IntegrationTypes integrationType, typename T::SectionType sectionType>
    userLibrary::BftElementFactory::elementFactoryFunction makeFactoryFunction()
    {
        return []( int elementID ) -> BftElement* { return new T( elementID, integrationType, sectionType ); };
    }

    using namespace userLibrary;
    using namespace bft::NumIntegration;

    const static bool UelCPS4_isRegistered = BftElementFactory::
        registerElement( "UELCPS4",
                         userLibrary::ElementCode::UelCPS4,
                         makeFactoryFunction<UelDisplacement<2, 4>,
                                             FullIntegration,
                                             UelDisplacement<2, 4>::PlaneStress>() );

    const static bool UelCPS8_isRegistered = BftElementFactory::
        registerElement( "UELCPS8",
                         userLibrary::ElementCode::UelCPS8,
                         makeFactoryFunction<UelDisplacement<2, 8>,
                                             FullIntegration,
                                             UelDisplacement<2, 8>::PlaneStress>() );

    const static bool UelCPE4_isRegistered = BftElementFactory::
        registerElement( "UELCPE4",
                         userLibrary::ElementCode::UelCPE4,
                         makeFactoryFunction<UelDisplacement<2, 4>,
                                             FullIntegration,
                                             UelDisplacement<2, 4>::PlaneStrain>() );

    const static bool UelCPE8_isRegistered = BftElementFactory::
        registerElement( "UELCPE8",
                         userLibrary::ElementCode::UelCPE8,
                         makeFactoryFunction<UelDisplacement<2, 8>,
                                             FullIntegration,
                                             UelDisplacement<2, 8>::PlaneStrain>() );
    const static bool UelCPS8R_isRegistered = BftElementFactory::
        registerElement( "UELCPS8R",
                         userLibrary::ElementCode::UelCPS8R,
                         makeFactoryFunction<UelDisplacement<2, 8>,
                                             ReducedIntegration,
                                             UelDisplacement<2, 8>::PlaneStress>() );
    const static bool UelCPE8R_isRegistered = BftElementFactory::
        registerElement( "UELCPE8R",
                         userLibrary::ElementCode::UelCPE8R,
                         makeFactoryFunction<UelDisplacement<2, 8>,
                                             ReducedIntegration,
                                             UelDisplacement<2, 8>::PlaneStrain>() );

    const static bool UelC3D8_isRegistered = userLibrary::BftElementFactory::
        registerElement( "UELC3D8", userLibrary::ElementCode::UelC3D8, []( int elementID ) -> BftElement* {
            return new UelDisplacement<3, 8>( elementID,
                                              bft::NumIntegration::IntegrationTypes::FullIntegration,
                                              UelDisplacement<3, 8>::SectionType::Solid );
        } );

    const static bool UelC3D20_isRegistered = userLibrary::BftElementFactory::
        registerElement( "UELC3D20", userLibrary::ElementCode::UelC3D20, []( int elementID ) -> BftElement* {
            return new UelDisplacement<3, 20>( elementID,
                                               bft::NumIntegration::IntegrationTypes::FullIntegration,
                                               UelDisplacement<3, 20>::SectionType::Solid );
        } );

    const static bool UelC3D20R_isRegistered = userLibrary::BftElementFactory::
        registerElement( "UELC3D20R", userLibrary::ElementCode::UelC3D20R, []( int elementID ) -> BftElement* {
            return new UelDisplacement<3, 20>( elementID,
                                               bft::NumIntegration::IntegrationTypes::ReducedIntegration,
                                               UelDisplacement<3, 20>::SectionType::Solid );
        } );

    BftElement* generateUelT2D2( int elementID )
    {
        auto uelT2D2 = std::unique_ptr<BftElement>(
            new UelDisplacement<1, 2>( elementID,
                                       bft::NumIntegration::IntegrationTypes::FullIntegration,
                                       UelDisplacement<1, 2>::SectionType::UniaxialStress ) );
        constexpr static int indicesToBeWrapped[] = {0, 1};
        constexpr static int nIndicesToBeWrapped  = 2;
        return new BftElementSpatialWrapper( 2,
                                             1,
                                             2,
                                             2,
                                             indicesToBeWrapped,
                                             nIndicesToBeWrapped,
                                             std::move( uelT2D2 ) );
    }
    const static bool
        UelT2D2_isRegistered = userLibrary::BftElementFactory::registerElement( "UELT2D2",
                                                                                userLibrary::ElementCode::UelT2D2,
                                                                                generateUelT2D2 );
} // namespace UelDisplacementRegistration
