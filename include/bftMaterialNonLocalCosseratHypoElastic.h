#pragma once
#include "bftMaterial.h"
#include "bftTypedefs.h"

class BftMaterialNonLocalCosseratHypoElastic : public BftMaterial {

    /*
       Abstract basic class for nonlocal Cosserat materials

       Formulated incrementally as σ_np, K_local_np = f (σ_n, dxdX_n, dxdX_np, Δt, t_n, , Kn, ΔK, K_local_n.. )

       Algorithmic tangents: dσdF = d σ_np d (dxdX_np)
                             dK_LocaldF = d K_local_np d (dxdX_np)
                             dσdK = d σ_np d ΔK
    */
        static constexpr int getNumberOfDofForRotation( int nDim )
        {
            if ( nDim == 2 )
                return 1;
            else
                return 3;
        }

  public:

      template <int nDim>
    struct ConstitutiveResponse {
        Eigen::Matrix<double, nDim, nDim> stress;
        Eigen::Matrix<double, nDim, getNumberOfDofForRotation(nDim)> coupleStress;
        double localField;
        double nonLocalRadius;
    };

      template <int nDim>
    struct AlgorithmicModuli {
        Eigen::TensorFixedSize< double, Eigen::Sizes< nDim, nDim, nDim, nDim>> dStressDDStrain;
        Eigen::TensorFixedSize< double, Eigen::Sizes< nDim, nDim, nDim, getNumberOfDofForRotation(nDim)>> dStressDDCurvature;
        Eigen::Matrix<double, nDim, nDim> dStressDNonLocalField;
        Eigen::TensorFixedSize< double, Eigen::Sizes< nDim, getNumberOfDofForRotation(nDim),  nDim, nDim>> dCoupleStressDDStrain;
        Eigen::TensorFixedSize< double, Eigen::Sizes< nDim, getNumberOfDofForRotation(nDim),  nDim,  getNumberOfDofForRotation(nDim) >> dCoupleStressDDCurvature;
        Eigen::Matrix<double, nDim, getNumberOfDofForRotation(nDim)> dCoupleStressDNonLocalField;
        Eigen::Matrix<double, nDim, nDim> dLocalFieldDDStrain;
        Eigen::Matrix<double, nDim, getNumberOfDofForRotation(nDim)> dLocalFieldDDCurvature;
        double dLocalFieldDNonLocalField;
    };

    template <int nDim>
    struct DeformationIncrement {
        Eigen::Matrix<double, nDim, nDim> dStrain;
        Eigen::Matrix<double, nDim, getNumberOfDofForRotation(nDim) > dCurvature;
        double nonLocalField;
    };

    struct TimeIncrement {
        const double* timeOld;
        const double  dT;
    };

    using BftMaterial::BftMaterial;

    virtual void computeStress( ConstitutiveResponse<3>&,
                                AlgorithmicModuli<3>&,
                                const DeformationIncrement<3>&,
                                const TimeIncrement&,
                                double& pNewDT ) = 0;

    virtual void computePlaneStrain( ConstitutiveResponse<2>&,
                                AlgorithmicModuli<2>&,
                                const DeformationIncrement<2>&,
                                const TimeIncrement&,
                                double& pNewDT );

    //virtual void computePlaneStress( ConstitutiveResponse*,
                                //AlgorithmicModuli*,
                                //DeformationIncrement*,
                                //const TimeIncrement*,
                                //double& pNewDT );
};
