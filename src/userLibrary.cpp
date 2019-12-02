#include "bftElement.h"
#include "bftMaterial.h"
#include "bftTypedefs.h"
#include "bftVoigt.h"
#include "userLibrary.h"
#include <map>
#include <string>
#include <tuple>

namespace userLibrary {

    // MaterialFactory

    std::map<std::string, MaterialCode> BftMaterialFactory::materialNameToCodeAssociation;
    std::map<MaterialCode, BftMaterialFactory::materialFactoryFunction> BftMaterialFactory::materialFactoryFunctionByCode;

    bool BftMaterialFactory::registerMaterial( MaterialCode            materialCode,
                                               const std::string&      materialName,
                                               materialFactoryFunction factoryFunction )
    {
        assert( materialNameToCodeAssociation.find( materialName ) != materialNameToCodeAssociation.end() );
        assert( materialFactoryFunctionByCode.find( materialCode ) != materialFactoryFunctionByCode.end() );

        materialNameToCodeAssociation[materialName] = materialCode;
        materialFactoryFunctionByCode[materialCode] = factoryFunction;

        return true;
    }

    MaterialCode BftMaterialFactory::getMaterialCodeFromName( const std::string& materialName )
    {
        try {
            return materialNameToCodeAssociation.at( materialName );
        }
        catch ( const std::out_of_range& e ) {
            throw std::invalid_argument( MakeString() << "Invalid material " << materialName << " requested!" );
        }
    }

    BftMaterial* BftMaterialFactory::createMaterial( MaterialCode  materialCode,
                                                     const double* materialProperties,
                                                     int           nMaterialProperties,
                                                     int           element,
                                                     int           gaussPt )
    {
        try {
            return materialFactoryFunctionByCode.at( materialCode )( materialProperties, nMaterialProperties, element, gaussPt );
        }
        catch ( const std::out_of_range& e ) {
            throw std::invalid_argument( MakeString() << "Invalid material " << materialCode << " requested!" );
        }
    }

    // ElementFactory

    std::map<std::string, ElementCode>                               BftElementFactory::elementNameToCodeAssociation;
    std::map<ElementCode, BftElementFactory::elementFactoryFunction> BftElementFactory::elementFactoryFunctionByCode;

    bool BftElementFactory::registerElement( const std::string&     elementName,
                                             ElementCode            elementCode,
                                             elementFactoryFunction factoryFunction )
    {
        assert( elementNameToCodeAssociation.find( elementName ) != elementNameToCodeAssociation.end() );
        assert( elementFactoryFunctionByCode.find( elementCode ) != elementFactoryFunctionByCode.end() );

        elementNameToCodeAssociation[elementName] = elementCode;
        elementFactoryFunctionByCode[elementCode] = factoryFunction;

        return true;
    }

    ElementCode BftElementFactory::getElementCodeFromName( const std::string& elementName )
    {
        try {
            return elementNameToCodeAssociation.at( elementName );
        }
        catch ( const std::out_of_range& e ) {
            throw std::invalid_argument( MakeString() << "Invalid element " << elementName << " requested!" );
        }
    }

    BftElement* BftElementFactory::createElement( ElementCode elementCode, int elementNumber )
    {
        try {
            return elementFactoryFunctionByCode.at( elementCode )( elementNumber );
        }
        catch ( const std::out_of_range& e ) {
            throw std::invalid_argument( MakeString() << "Invalid element " << elementCode << " requested!" );
        }
    }

} // namespace userLibrary
