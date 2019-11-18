#include "bftElement.h"
#include "bftMaterial.h"
#include "bftTypedefs.h"
#include "bftVoigt.h"
#include "userLibrary.h"
#include <map>
#include <string>
#include <tuple>

namespace userLibrary {

    std::map<std::string, MaterialCode> BftMaterialFactory::materialNameToCodeAssociation;
    std::map<MaterialCode, BftMaterialFactory::materialCreationFunction>
        BftMaterialFactory::materialCreationFunctionByCode;

    bool BftMaterialFactory::registerMaterial( const std::string&       materialName,
                                               MaterialCode             materialCode,
                                               materialCreationFunction creationFunction )
    {
        materialNameToCodeAssociation[materialName]  = materialCode;
        materialCreationFunctionByCode[materialCode] = creationFunction;

        return true;
    }

    MaterialCode BftMaterialFactory::getMaterialCodeFromName( const std::string& materialName )
    {
        return materialNameToCodeAssociation[materialName];
    }

    BftMaterial* BftMaterialFactory::createMaterial( MaterialCode  materialCode,
                                                     const double* materialProperties,
                                                     int           nMaterialProperties,
                                                     int           element,
                                                     int           gaussPt )
    {
        auto const it = materialCreationFunctionByCode.find( materialCode );
        return it != materialCreationFunctionByCode.end()
                   ? it->second( materialProperties, nMaterialProperties, element, gaussPt )
                   : nullptr;
    }

    std::map<std::string, ElementCode>                                BftElementFactory::elementNameToCodeAssociation;
    std::map<ElementCode, BftElementFactory::elementCreationFunction> BftElementFactory::elementCreationFunctionByCode;

    bool BftElementFactory::registerElement( const std::string&      elementName,
                                             ElementCode             elementCode,
                                             elementCreationFunction creationFunction )
    {
        elementNameToCodeAssociation[elementName]  = elementCode;
        elementCreationFunctionByCode[elementCode] = creationFunction;

        return true;
    }

    ElementCode BftElementFactory::getElementCodeFromName( const std::string& elementName )
    {
        return elementNameToCodeAssociation[elementName];
    }

    BftElement* BftElementFactory::createElement( ElementCode elementCode, int elementNumber )
    {
        auto const it = elementCreationFunctionByCode.find( elementCode );
        return it != elementCreationFunctionByCode.end() ? it->second( elementNumber ) : nullptr;
    }

} // namespace userLibrary
