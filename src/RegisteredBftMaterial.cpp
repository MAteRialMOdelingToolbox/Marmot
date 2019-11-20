#include "RegisteredBftMaterial.h"
namespace userLibrary
{
    template <typename T>
    bool RegisteredBftMaterial<T>::isRegistered = BftMaterialFactory::registerMaterial(
        T::materialCode,
        T::materialName,
        []( const double* materialProperties, int nMaterialProperties, int element, int gaussPt ) -> BftMaterial* {
            return new T( materialProperties, nMaterialProperties, element, gaussPt );
        }

    );
}
