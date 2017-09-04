#pragma once 
#include "uelDisplacement.h"
#include "bftAbaqusUmatWrapper.h"

template <int nNodes> 
class UelDisplacementSolid: public UelDisplacement<3, nNodes>
{
    /* Derived Solid Element 
     * */
    static constexpr int nDim = 3;
    typedef UelDisplacement<nDim, nNodes> ParentUelDisplacement;

    public:

        using  UelDisplacement<nDim,nNodes>::UelDisplacement; 
        void computeYourself( const double* QTotal,
                                            const double* dQ,
                                            double* Pe,
                                            double* Ke,
                                            const double* time,
                                            double dT,
                                            double& pNewDT);
};

template <int nNodes>
void UelDisplacementSolid<nNodes>::computeYourself( const double* QTotal_,
                                    const double* dQ_,
                                    double* Pe_,
                                    double* Ke_,
                                    const double* time,
                                    double dT,
                                    double& pNewDT
                                   ){
    using namespace bft;

    Map<const typename ParentUelDisplacement::RhsSized>                  QTotal(QTotal_);				
    Map<const typename ParentUelDisplacement::RhsSized>                  dQ(dQ_);			
    Map<typename ParentUelDisplacement::KeSizedMatrix>                   Ke(Ke_);
    Map<typename ParentUelDisplacement::RhsSized>                        Pe(Pe_);

    Ke.setZero();
    Pe.setZero();

    Vector6 dStrain;
    Matrix6 Cep;

    for(size_t i = 0; i < this->gaussPts.size() ; i++){
        
        typename ParentUelDisplacement::GaussPt& gaussPt = this->gaussPts[i];
        const typename ParentUelDisplacement::ParentGeometryElement::BSized& B = gaussPt.B;

        const double detJ =         gaussPt.detJ;
        const double charElemlen =  std::cbrt(8*detJ);

        gaussPt.material->setCharacteristicElementLength( charElemlen ); // could be moved to a constructor for solid (like plane)

        dStrain = B * dQ;
        Cep.setZero();

        gaussPt.material->computeStress(gaussPt.stress.data(), 
                                                Cep.data(),  
                                                gaussPt.strain.data(),
                                                dStrain.data(), 
                                                time, dT, pNewDT);
        if (pNewDT<1.0)
            return;

        const double vol = detJ * gaussPt.weight; // could be moved to a constructor for solid

        Ke += B.transpose() * Cep * B * vol;
        Pe -= B.transpose() * gaussPt.stress* vol;

        gaussPt.strain += dStrain; 
    }
}
