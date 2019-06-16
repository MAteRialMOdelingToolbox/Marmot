#pragma once
#include "bftMaterialMechanical.h"

class BftMaterialHypoElastic : public BftMaterialMechanical {

    // Derived abstract base class for elastic materials expressed purely in rate form in terms of stretching rate d, i.e,
    // 'hypoelastic materials':
    //
    // ∇ 
    // σ = f (σ, d, t, .. ), 
    //
    // formulated incrementally as σ_np = f (σ_n, Δε, Δt, t_n, .. ) 
    // with Δε = d * Δt
    //
    // Algorithmic tangent: dσdε = d Δσ d Δε 
    //
    // compatible with Abaqus interface 

  public:
    using BftMaterialMechanical::BftMaterialMechanical;

    double characteristicElementLength;
    void setCharacteristicElementLength( double length );

    virtual void computeStress( double*       stress,
                                double*       dStressDDStrain,
                                const double* FOld,
                                const double* FNew,
                                const double* timeOld,
                                const double  dT,
                                double&       pNewDT ) override;

    // Abstract methods
    virtual void computeStress( double*       stress,
                                double*       dStressDDStrain,
                                const double* dStrain,
                                const double* timeOld,
                                const double  dT,
                                double&       pNewDT ) = 0;

    virtual void computePlaneStress( double*       stress,
                                     double*       dStressDDStrain,
                                     double*       dStrain,
                                     const double* timeOld,
                                     const double  dT,
                                     double&       pNewDT );

    virtual void computeUniaxialStress( double*       stress,
                                        double*       dStressDDStrain,
                                        double*       dStrain,
                                        const double* timeOld,
                                        const double  dT,
                                        double&       pNewDT );
};
