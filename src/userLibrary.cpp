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
    std::map<MaterialCode, BftMaterialFactory::materialFactoryFunction>
        BftMaterialFactory::materialFactoryFunctionByCode;

    bool BftMaterialFactory::registerMaterial( MaterialCode            materialCode,
                                               const std::string&      materialName,
                                               materialFactoryFunction factoryFunction )
    {
        materialNameToCodeAssociation[materialName] = materialCode;
        materialFactoryFunctionByCode[materialCode] = factoryFunction;

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
        auto const it = materialFactoryFunctionByCode.find( materialCode );
        return it != materialFactoryFunctionByCode.end()
                   ? it->second( materialProperties, nMaterialProperties, element, gaussPt )
                   : nullptr;
    }

    // ElementFactory

    std::map<std::string, ElementCode>                               BftElementFactory::elementNameToCodeAssociation;
    std::map<ElementCode, BftElementFactory::elementFactoryFunction> BftElementFactory::elementFactoryFunctionByCode;

    bool BftElementFactory::registerElement( const std::string&     elementName,
                                             ElementCode            elementCode,
                                             elementFactoryFunction factoryFunction )
    {
        elementNameToCodeAssociation[elementName] = elementCode;
        elementFactoryFunctionByCode[elementCode] = factoryFunction;

        return true;
    }

    ElementCode BftElementFactory::getElementCodeFromName( const std::string& elementName )
    {
        return elementNameToCodeAssociation[elementName];
    }

    BftElement* BftElementFactory::createElement( ElementCode elementCode, int elementNumber )
    {
        auto const it = elementFactoryFunctionByCode.find( elementCode );
        return it != elementFactoryFunctionByCode.end() ? it->second( elementNumber ) : nullptr;
    }

} // namespace userLibrary
