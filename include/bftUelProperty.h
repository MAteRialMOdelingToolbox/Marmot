#pragma once
//#include "userLibrary.h"
//
namespace userLibrary{ enum MaterialCode : int;}

class BftUelProperty {
    public: 
        enum Type { 
            ElementProperties,
            BftMaterialSection };

        const Type type;

    //protected: // protected since only childs should be instantiated

        BftUelProperty(Type type) : type (type) {};
        virtual ~BftUelProperty() {};
};

class BftMaterialSection : public BftUelProperty {
    public:

        userLibrary::MaterialCode   materialCode; 
        const double*               materialProperties;
        int                         nMaterialProperties;

    BftMaterialSection(userLibrary::MaterialCode materialCode, const double* materialProperties, int nMaterialProperties):
        BftUelProperty( Type:: BftMaterialSection ),
        materialCode(materialCode), 
        materialProperties(materialProperties), 
        nMaterialProperties(nMaterialProperties) {};
};

class ElementProperties : public BftUelProperty {
    public:
        const double*               elementProperties;
        int                         nElementProperties;

    ElementProperties(const double* elementProperties, int nElementProperties):
        BftUelProperty( Type::ElementProperties ),
        elementProperties(elementProperties), 
        nElementProperties(nElementProperties) {};
};
