#pragma once
#include "bftMaterialElastic.h"

class BftMaterialHyperElastic : public BftMaterialElastic {

    /* Derived abstract base class for a simple, purely hyperelastic material to be used within TL elements
     * stress measure: Piola - Kirchhoff II .. S
     * strain measure for algorithmic tangent: Green - Lagrange .. E = 1/2 ( F^T * F - I )
     * Algorithmic tangent: dS dE
     * */

  public:
    using BftMaterialElastic::BftMaterialElastic;

    virtual void computeStress( double*       S,    // PK2
                                double*       dSdE, // d PK2 d GL_E
                                const double* FOld,
                                const double* FNew,
                                const double* timeOld,
                                const double  dT,
                                double&       pNewDT ) override;

    // Abstract methods
    virtual void computeStressPK2( double*       S,    // PK2
                                   double*       dSdE, // d PK2 d GL_E
                                   const double* E,
                                   const double* timeOld,
                                   const double  dT,
                                   double&       pNewDT ) = 0;

    virtual void computePlaneStressPK2( double*       S,
                                        double*       dSdE,
                                        double*       E,
                                        const double* timeOld,
                                        const double  dT,
                                        double&       pNewDT );

    virtual void computeUniaxialStressPK2( double*       S,
                                           double*       dSdE,
                                           double*       E,
                                           const double* timeOld,
                                           const double  dT,
                                           double&       pNewDT );
};
