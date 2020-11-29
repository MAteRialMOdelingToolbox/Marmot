#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotFiniteElementSpatialWrapper.h"
#include "Marmot/DisplacementFiniteElement.h"
#include "Marmot/userLibrary.h"

namespace DisplacementFiniteElementRegistration {

    template <class T, Marmot::FiniteElement::Quadrature::IntegrationTypes integrationType, typename T::SectionType sectionType>
    userLibrary::MarmotElementFactory::elementFactoryFunction makeFactoryFunction()
    {
        return []( int elementID ) -> MarmotElement* { return new T( elementID, integrationType, sectionType ); };
    }

    using namespace userLibrary;
    using namespace Marmot::FiniteElement::Quadrature;

    const static bool CPS4_isRegistered = MarmotElementFactory::
        registerElement( "UELCPS4",
                         userLibrary::ElementCode::CPS4,
                         makeFactoryFunction<DisplacementFiniteElement<2, 4>,
                                             FullIntegration,
                                             DisplacementFiniteElement<2, 4>::PlaneStress>() );

    const static bool CPE4_isRegistered = MarmotElementFactory::
        registerElement( "UELCPE4",
                         userLibrary::ElementCode::CPE4,
                         makeFactoryFunction<DisplacementFiniteElement<2, 4>,
                                             FullIntegration,
                                             DisplacementFiniteElement<2, 4>::PlaneStrain>() );

    const static bool CPS8R_isRegistered = MarmotElementFactory::
        registerElement( "UELCPS8R",
                         userLibrary::ElementCode::CPS8R,
                         makeFactoryFunction<DisplacementFiniteElement<2, 8>,
                                             ReducedIntegration,
                                             DisplacementFiniteElement<2, 8>::PlaneStress>() );

    const static bool CPE8R_isRegistered = MarmotElementFactory::
        registerElement( "UELCPE8R",
                         userLibrary::ElementCode::CPE8R,
                         makeFactoryFunction<DisplacementFiniteElement<2, 8>,
                                             ReducedIntegration,
                                             DisplacementFiniteElement<2, 8>::PlaneStrain>() );

    const static bool C3D8_isRegistered = userLibrary::MarmotElementFactory::
        registerElement( "UELC3D8", userLibrary::ElementCode::C3D8, []( int elementID ) -> MarmotElement* {
            return new DisplacementFiniteElement<3, 8>( elementID,
                                              Marmot::FiniteElement::Quadrature::IntegrationTypes::FullIntegration,
                                              DisplacementFiniteElement<3, 8>::SectionType::Solid );
        } );

    const static bool C3D20_isRegistered = userLibrary::MarmotElementFactory::
        registerElement( "UELC3D20", userLibrary::ElementCode::C3D20, []( int elementID ) -> MarmotElement* {
            return new DisplacementFiniteElement<3, 20>( elementID,
                                               Marmot::FiniteElement::Quadrature::IntegrationTypes::FullIntegration,
                                               DisplacementFiniteElement<3, 20>::SectionType::Solid );
        } );

    const static bool C3D20R_isRegistered = userLibrary::MarmotElementFactory::
        registerElement( "UELC3D20R", userLibrary::ElementCode::C3D20R, []( int elementID ) -> MarmotElement* {
            return new DisplacementFiniteElement<3, 20>( elementID,
                                               Marmot::FiniteElement::Quadrature::IntegrationTypes::ReducedIntegration,
                                               DisplacementFiniteElement<3, 20>::SectionType::Solid );
        } );

    MarmotElement* generateT2D2( int elementID )
    {
        auto uelT2D2 = std::unique_ptr<MarmotElement>(
            new DisplacementFiniteElement<1, 2>( elementID,
                                       Marmot::FiniteElement::Quadrature::IntegrationTypes::FullIntegration,
                                       DisplacementFiniteElement<1, 2>::SectionType::UniaxialStress ) );
        constexpr static int indicesToBeWrapped[] = {0, 1};
        constexpr static int nIndicesToBeWrapped  = 2;
        return new MarmotElementSpatialWrapper( 2,
                                             1,
                                             2,
                                             2,
                                             indicesToBeWrapped,
                                             nIndicesToBeWrapped,
                                             std::move( uelT2D2 ) );
    }
    const static bool
        T2D2_isRegistered = userLibrary::MarmotElementFactory::registerElement( "UELT2D2",
                                                                                userLibrary::ElementCode::T2D2,
                                                                                generateT2D2 );
} // namespace DisplacementFiniteElementRegistration
