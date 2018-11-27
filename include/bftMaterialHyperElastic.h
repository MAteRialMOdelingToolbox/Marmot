#pragma once
#include "bftMaterial.h"

class BftMaterialHyperElastic : public BftMaterial {

    /* Interface for a pure hyperelastic material
     * stress measure: Piola - Kirchhoff II
     * strain measure for algorithmic tangent: Green - Lagrange
     * Deformation gradients @ beginning of the time increment and @ end of time increment are passed
     * */

  public:
    using BftMaterial::BftMaterial;

    // Abstract methods
    virtual void computeStress( double*       S,    // PK2
                                double*       dSdE, // d PK2 d GL_E
                                const double* E,
                                const double* F,
                                const double* timeOld,
                                const double  dT,
                                double&       pNewDT ) = 0;

    virtual void computePlaneStress( double*       S,
                                     double*       dSdE,
                                     const double* E,
                                     const double* F,
                                     const double* timeOld,
                                     const double  dT,
                                     double&       pNewDT );

    virtual void computeUniaxialStress( double*       S,
                                        double*       dSdE,
                                        const double* E,
                                        const double* F,
                                        const double* timeOld,
                                        const double  dT,
                                        double&       pNewDT );
};
