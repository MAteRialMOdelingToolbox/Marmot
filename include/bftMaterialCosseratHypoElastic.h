#pragma once
#include "bftMaterial.h"

class BftMaterialCosseratHypoElastic : public BftMaterial {

    /*
       Abstract basic class for Cosserat materials

       Formulated incrementally as σ_np, K_local_np = f (σ_n, dxdX_n, dxdX_np, Δt, t_n, , Kn, ΔK, K_local_n.. )

       Algorithmic tangents: dσdF = d σ_np d (dxdX_np)
                             dK_LocaldF = d K_local_np d (dxdX_np)
                             dσdK = d σ_np d ΔK
    */

  public:
    struct ConstitutiveResponse {
        double* stress;
        double* coupleStress;
    };

    struct AlgorithmicModuli {
        double* dStressDDStrain;
        double* dStressDDCurvature;
        double* dCoupleStressDDStrain;
        double* dCoupleStressDDCurvature;
    };

    struct DeformationIncrement {
        const double* dStrain;
        const double* dCurvature;
    };

    struct TimeIncrement {
        const double* timeOld;
        const double  dT;
    };

    using BftMaterial::BftMaterial;

    double characteristicElementLength;
    void   setCharacteristicElementLength( double length ) { characteristicElementLength = length; }

    virtual void computeStress( ConstitutiveResponse*       response,
                                AlgorithmicModuli*          algorithmicModuli,
                                const DeformationIncrement* deformationIncrement,
                                const TimeIncrement*        timeIncrement,
                                double&                     pNewDT ) = 0;

    //virtual void computePlaneStress( ConstitutiveResponse* response,
                                     //AlgorithmicModuli*    algorithmicModuli,
                                     //DeformationIncrement* deformationIncrement,
                                     //const TimeIncrement*  timeIncrement,
                                     //double&               pNewDT );
};
