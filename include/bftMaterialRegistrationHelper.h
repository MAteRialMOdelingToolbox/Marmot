#pragma once
#include "userLibrary.h"

template <typename T>
userLibrary::BftMaterialFactory::materialFactoryFunction makeDefaultBftMaterialFactoryFunction()
{
    return
        []( const double* materialProperties, int nMaterialProperties, int element, int gaussPt ) -> BftMaterial* {
            return new T( materialProperties, nMaterialProperties, element, gaussPt );
        };
}

