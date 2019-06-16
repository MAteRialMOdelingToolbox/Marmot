#pragma once
#include "bftMaterial.h"

class BftMaterialMechanical : public BftMaterial {

    // Abstract basic class for Elastic materials.
    // 'Elastic' is meant in the 'most general sense', i.e., any material which describes a mechanical (cauchy)
    // stress - deformation relationship (e.g, elastic, elasto-plastic, visco-elastic materials)
    //
    // σ = f (σ, dxdX, t, .. ),
    //
    // formulated incrementally as σ_np = f (σ_n, dxdX_n, dxdX_np, Δt, t_n, .. )
    //
    // Algorithmic tangent: dσdε = d σ d [ sym( ( inv(F_n) * F_np) - I) ] (Abaqus compatible)

  public:
    using BftMaterial::BftMaterial;

    virtual void computeStress( double*       stress,
                                double*       dStressDDStrain,
                                const double* FOld,
                                const double* FNew,
                                const double* timeOld,
                                const double  dT,
                                double&       pNewDT ) = 0;

    virtual void computePlaneStress( double*       stress,
                                     double*       dStressDDStrain,
                                     const double* FOld,
                                     double*       FNew,
                                     const double* timeOld,
                                     const double  dT,
                                     double&       pNewDT );

    virtual void computeUniaxialStress( double*       stress,
                                        double*       dStressDDStrain,
                                        const double* FOld,
                                        double*       FNew,
                                        const double* timeOld,
                                        const double  dT,
                                        double&       pNewDT );
};
