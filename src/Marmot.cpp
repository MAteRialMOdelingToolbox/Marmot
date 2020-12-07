#include "Marmot/MarmotElement.h"
#include "Marmot/MarmotMaterial.h"
#include "Marmot/Marmot.h"
#include "Marmot/MarmotJournal.h"
#include <map>
#include <string>
#include <tuple>
#include <cassert>

namespace userLibrary {

    // MaterialFactory

    std::map<std::string, MaterialCode> MarmotMaterialFactory::materialNameToCodeAssociation;
    std::map<MaterialCode, MarmotMaterialFactory::materialFactoryFunction> MarmotMaterialFactory::materialFactoryFunctionByCode;

    bool MarmotMaterialFactory::registerMaterial( MaterialCode            materialCode,
                                               const std::string&      materialName,
                                               materialFactoryFunction factoryFunction )
    {
        assert( materialNameToCodeAssociation.find( materialName ) == materialNameToCodeAssociation.end() );
        assert( materialFactoryFunctionByCode.find( materialCode ) == materialFactoryFunctionByCode.end() );

        materialNameToCodeAssociation[materialName] = materialCode;
        materialFactoryFunctionByCode[materialCode] = factoryFunction;

        return true;
    }

    MaterialCode MarmotMaterialFactory::getMaterialCodeFromName( const std::string& materialName )
    {
        try {
            return materialNameToCodeAssociation.at( materialName );
        }
        catch ( const std::out_of_range& e ) {
            throw std::invalid_argument( MakeString() << "Invalid material " << materialName << " requested!" );
        }
    }

    MarmotMaterial* MarmotMaterialFactory::createMaterial( MaterialCode materialCode, const double* materialProperties, int nMaterialProperties, int materialNumber)
    {
        try {
            return materialFactoryFunctionByCode.at( materialCode )( materialProperties, nMaterialProperties, materialNumber );
        }
        catch ( const std::out_of_range& e ) {
            throw std::invalid_argument( MakeString() << "Invalid material " << materialCode << " requested!" );
        }
    }

    // ElementFactory

    std::map<std::string, ElementCode>                               MarmotElementFactory::elementNameToCodeAssociation;
    std::map<ElementCode, MarmotElementFactory::elementFactoryFunction> MarmotElementFactory::elementFactoryFunctionByCode;

    bool MarmotElementFactory::registerElement( const std::string&     elementName,
                                             ElementCode            elementCode,
                                             elementFactoryFunction factoryFunction )
    {
        assert( elementNameToCodeAssociation.find( elementName ) == elementNameToCodeAssociation.end() );
        assert( elementFactoryFunctionByCode.find( elementCode ) == elementFactoryFunctionByCode.end() );

        elementNameToCodeAssociation[elementName] = elementCode;
        elementFactoryFunctionByCode[elementCode] = factoryFunction;

        return true;
    }

    ElementCode MarmotElementFactory::getElementCodeFromName( const std::string& elementName )
    {
        try {
            return elementNameToCodeAssociation.at( elementName );
        }
        catch ( const std::out_of_range& e ) {
            throw std::invalid_argument( MakeString() << "Invalid element " << elementName << " requested!" );
        }
    }

    MarmotElement* MarmotElementFactory::createElement( ElementCode elementCode, int elementNumber )
    {
        try {
            return elementFactoryFunctionByCode.at( elementCode )( elementNumber );
        }
        catch ( const std::out_of_range& e ) {
            throw std::invalid_argument( MakeString() << "Invalid element " << elementCode << " requested!" );
        }
    }

} // namespace userLibrary
