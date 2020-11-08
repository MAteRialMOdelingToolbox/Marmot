#pragma once
#include "userLibrary.h"

template <typename T>
userLibrary::MarmotMaterialFactory::materialFactoryFunction makeDefaultMarmotMaterialFactoryFunction()
{
    return
        []( int materialNumber) -> MarmotMaterial* {
            return new T( materialNumber );
        };
}

