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

    for(int i = 0; i < this->gaussPointList.rows() ; i++){

        const typename ParentUelDisplacement::ParentGeometryElement::BSized& B =   this->BAtGauss[i]; 
        const double detJ =         this->detJAtGauss[i];
        const double charElemlen =  std::cbrt(8*detJ);

        this->materialAtGauss[i]->setCharacteristicElementLength( charElemlen );

        Ref<Vector6>    stress(             this->stressAtGauss(i));
        Ref<Vector6>    strain(             this->strainAtGauss(i));

        dStrain = B * dQ;
        Cep.setZero();

        this->materialAtGauss[i]->computeStress(stress.data(), 
                                                Cep.data(),  
                                                strain.data(),
                                                dStrain.data(), 
                                                time, dT, pNewDT);

        if (pNewDT<1.0)
            return;

        const double vol = detJ * this->gaussWeights(i);

        Ke += B.transpose() * Cep * B * vol;
        Pe -= B.transpose() * stress* vol;

        strain += dStrain; 
    }
}
