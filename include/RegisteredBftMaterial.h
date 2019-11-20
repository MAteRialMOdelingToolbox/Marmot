#pragma once
#include "userLibrary.h"

namespace userLibrary {
    template <typename T>
    class RegisteredBftMaterial {
        // Provides a convenient interface to derive from for BftMaterials for AUTOMATIC self-registration in the
        // BftMaterialFactory.
      protected:
        static bool isRegistered;
        RegisteredBftMaterial() {
            if ( !isRegistered )
                throw std::invalid_argument( "Failed to register material" );
        }
    };
} // namespace userLibrary
