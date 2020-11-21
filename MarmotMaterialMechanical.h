#pragma once
#include "Marmot/MarmotMaterial.h"

class MarmotMaterialMechanical : public MarmotMaterial {

    /*
       Abstract basic class for Mechanical materials with scalar nonlocal interaction.

       Formulated incrementally as σ_np, K_local_np = f (σ_n, dxdX_n, dxdX_np, Δt, t_n, , Kn, ΔK, K_local_n.. )

       Algorithmic tangents: dσdF = d σ_np d (dxdX_np)
                             dK_LocaldF = d K_local_np d (dxdX_np)
                             dσdK = d σ_np d ΔK
    */

    public:

        using MarmotMaterial::MarmotMaterial;

        virtual void computeStress( double*       stress,
                double*       dStressDDFNew,
                const double* FOld,
                const double* FNew,
                const double* timeOld,
                const double  dT,
                double&       pNewDT ) = 0;

        virtual void computePlaneStress( double*       stress,
                double*       dStressDDFNew,
                const double* FOld,
                double*       FNew,
                const double* timeOld,
                const double  dT,
                double&       pNewDT );

        virtual void computeUniaxialStress( double*       stress,
                double*       dStressDDFNew,
                const double* FOld,
                double*       FNew,
                const double* timeOld,
                const double  dT,
                double&       pNewDT );
};
