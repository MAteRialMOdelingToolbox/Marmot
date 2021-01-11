#pragma once
#include "Marmot/MarmotMaterialMechanical.h"

class MarmotMaterialHypoElastic : public MarmotMaterialMechanical {

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

    /**
     * \brief Base class for elastic materials expressed purely in rate form in terms of stretching rate i.e 'hypoelastic materials' 
     */

  public:
    using MarmotMaterialMechanical::MarmotMaterialMechanical;
    
    /// \brief Characteristic element length
    double characteristicElementLength;
    /**< #charactersisticElementLength represents the characteristic length of a finite element.
     * It is needed for the regularization of materials with softening behavior to assure mesh independent results.
     */
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

    using MarmotMaterialMechanical::computePlaneStress;
    virtual void computePlaneStress( double*       stress,
                                     double*       dStressDDStrain,
                                     double*       dStrain,
                                     const double* timeOld,
                                     const double  dT,
                                     double&       pNewDT );

    using MarmotMaterialMechanical::computeUniaxialStress;
    virtual void computeUniaxialStress( double*       stress,
                                        double*       dStressDDStrain,
                                        double*       dStrain,
                                        const double* timeOld,
                                        const double  dT,
                                        double&       pNewDT );
};
