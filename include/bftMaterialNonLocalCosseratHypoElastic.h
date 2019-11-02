#pragma once
#include "bftMaterial.h"

class BftMaterialNonLocalCosseratHypoElastic : public BftMaterial {

    /*
       Abstract basic class for NonLocalCosserat materials

       Formulated incrementally as σ_np, K_local_np = f (σ_n, dxdX_n, dxdX_np, Δt, t_n, , Kn, ΔK, K_local_n.. )

       Algorithmic tangents: dσdF = d σ_np d (dxdX_np)
                             dK_LocaldF = d K_local_np d (dxdX_np)
                             dσdK = d σ_np d ΔK
    */

  public:
    struct ConstitutiveResponse {
        double* stress;
        double* coupleStress;
        double& localField;
        double& nonLocalRadius;
    };

    struct AlgorithmicModuli {
        double* dStressDDStrain;
        double* dStressDDCurvature;
        double* dStressDNonLocalField;
        double* dCoupleStressDDStrain;
        double* dCoupleStressDDCurvature;
        double* dCoupleStressDNonLocalField;
        double* dLocalFieldDDStrain;
        double* dLocalFieldDDCurvature;
        double& dLocalFieldDNonLocalField;
    };

    struct DeformationIncrement {
        const double* dStrain;
        const double* dCurvature;
        const double& nonLocalField;
    };

    struct TimeIncrement {
        const double* timeOld;
        const double  dT;
    };

    using BftMaterial::BftMaterial;

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
