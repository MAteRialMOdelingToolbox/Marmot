#include "Marmot/DisplacementFiniteElement.h"
#include "Marmot/Marmot.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotFiniteElementSpatialWrapper.h"

namespace Marmot::Elements::Registration {

    template < class T,
               Marmot::FiniteElement::Quadrature::IntegrationTypes integrationType,
               typename T::SectionType                             sectionType >
    MarmotLibrary::MarmotElementFactory::elementFactoryFunction makeFactoryFunction()
    {
        return []( int elementID ) -> MarmotElement* { return new T( elementID, integrationType, sectionType ); };
    }

    using namespace MarmotLibrary;
    using namespace Marmot::FiniteElement::Quadrature;

    const static bool CPS4_isRegistered = MarmotElementFactory::
        registerElement( "CPS4",
                         MarmotLibrary::ElementCode::CPS4,
                         makeFactoryFunction< DisplacementFiniteElement< 2, 4 >,
                                              FullIntegration,
                                              DisplacementFiniteElement< 2, 4 >::PlaneStress >() );

    const static bool CPE4_isRegistered = MarmotElementFactory::
        registerElement( "CPE4",
                         MarmotLibrary::ElementCode::CPE4,
                         makeFactoryFunction< DisplacementFiniteElement< 2, 4 >,
                                              FullIntegration,
                                              DisplacementFiniteElement< 2, 4 >::PlaneStrain >() );

    const static bool CPS8R_isRegistered = MarmotElementFactory::
        registerElement( "CPS8R",
                         MarmotLibrary::ElementCode::CPS8R,
                         makeFactoryFunction< DisplacementFiniteElement< 2, 8 >,
                                              ReducedIntegration,
                                              DisplacementFiniteElement< 2, 8 >::PlaneStress >() );

    const static bool CPE8R_isRegistered = MarmotElementFactory::
        registerElement( "CPE8R",
                         MarmotLibrary::ElementCode::CPE8R,
                         makeFactoryFunction< DisplacementFiniteElement< 2, 8 >,
                                              ReducedIntegration,
                                              DisplacementFiniteElement< 2, 8 >::PlaneStrain >() );

    const static bool C3D8_isRegistered = MarmotLibrary::MarmotElementFactory::
        registerElement( "C3D8", MarmotLibrary::ElementCode::C3D8, []( int elementID ) -> MarmotElement* {
            return new DisplacementFiniteElement< 3, 8 >( elementID,
                                                          Marmot::FiniteElement::Quadrature::IntegrationTypes::
                                                              FullIntegration,
                                                          DisplacementFiniteElement< 3, 8 >::SectionType::Solid );
        } );

    const static bool C3D20_isRegistered = MarmotLibrary::MarmotElementFactory::
        registerElement( "C3D20", MarmotLibrary::ElementCode::C3D20, []( int elementID ) -> MarmotElement* {
            return new DisplacementFiniteElement< 3, 20 >( elementID,
                                                           Marmot::FiniteElement::Quadrature::IntegrationTypes::
                                                               FullIntegration,
                                                           DisplacementFiniteElement< 3, 20 >::SectionType::Solid );
        } );

    const static bool C3D20R_isRegistered = MarmotLibrary::MarmotElementFactory::
        registerElement( "C3D20R", MarmotLibrary::ElementCode::C3D20R, []( int elementID ) -> MarmotElement* {
            return new DisplacementFiniteElement< 3, 20 >( elementID,
                                                           Marmot::FiniteElement::Quadrature::IntegrationTypes::
                                                               ReducedIntegration,
                                                           DisplacementFiniteElement< 3, 20 >::SectionType::Solid );
        } );

    MarmotElement* generateT2D2( int elementID )
    {
        auto uelT2D2 = std::unique_ptr< MarmotElement >(
            new DisplacementFiniteElement< 1, 2 >( elementID,
                                                   Marmot::FiniteElement::Quadrature::IntegrationTypes::FullIntegration,
                                                   DisplacementFiniteElement< 1, 2 >::SectionType::UniaxialStress ) );
        constexpr static int indicesToBeWrapped[] = { 0, 1 };
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
        T2D2_isRegistered = MarmotLibrary::MarmotElementFactory::registerElement( "T2D2",
                                                                                  MarmotLibrary::ElementCode::T2D2,
                                                                                  generateT2D2 );
} // namespace Marmot::Elements::Registration
