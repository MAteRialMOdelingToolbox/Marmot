#pragma once
#include "Marmot/userLibrary.h"

template <typename T>
userLibrary::MarmotMaterialFactory::materialFactoryFunction makeDefaultMarmotMaterialFactoryFunction()
{
    return
        []( const double* materialProperties, int nMaterialProperties, int materialNumber) -> MarmotMaterial* {
            return new T( materialProperties, nMaterialProperties, materialNumber );
        };
}

