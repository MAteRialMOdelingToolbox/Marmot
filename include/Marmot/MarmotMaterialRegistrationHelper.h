#pragma once
#include "Marmot/Marmot.h"

template <typename T>
MarmotLibrary::MarmotMaterialFactory::materialFactoryFunction makeDefaultMarmotMaterialFactoryFunction()
{
    return
        []( const double* materialProperties, int nMaterialProperties, int materialNumber) -> MarmotMaterial* {
            return new T( materialProperties, nMaterialProperties, materialNumber );
        };
}

