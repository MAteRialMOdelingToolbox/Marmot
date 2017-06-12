#pragma once
#include "uelDisplacement.h"
#include "bftAbaqusUmatWrapper.h"
#include "bftFiniteElement.h"

template <int nNodes> 
class UelDisplacementPlane : public UelDisplacement<2, nNodes>
{
    /* Derived Plane element, capable of plane stress and plane strain
     * */

    static constexpr int nDim = 2;
    typedef UelDisplacement<nDim, nNodes> ParentUelDisplacement;

    public:

    enum SectionType {
        PlaneStress,
        PlaneStrain,
    };

    const SectionType sectionType;

    UelDisplacementPlane(const double* coordinates,
            double* stateVars,
            int nStateVars,
            const double* propertiesElement,
            int nPropertiesElement,
            int noEl,
            const bft::pUmatType umat,
            int nStateVarsUmat,
            const double* propertiesUmat,
            int nPropertiesUmat,
            bft::NumIntegration::IntegrationTypes integrationType,
            SectionType sectionType): 
        ParentUelDisplacement(coordinates, stateVars, nStateVars, propertiesElement, nPropertiesElement, noEl, 
                umat, nStateVarsUmat, propertiesUmat, nPropertiesUmat, integrationType), sectionType(sectionType){};

    void computeDistributedLoad( BftUel::DistributedLoadTypes loadType,
            double* P, 
            const int elementFace, 
            const double* load,
            const double* time,
            double dT);

    void computeYourself( const double* QTotal,
            const double* dQ,
            double* Pe,
            double* Ke,
            const double* time,
            double dT,
            double& pNewdT);

};

template <int nNodes>
void UelDisplacementPlane<nNodes>::computeYourself( const double* QTotal_,
        const double* dQ_,
        double* Pe_,
        double* Ke_,
        const double* time,
        double dT,
        double& pNewdT
        ){

    using namespace bft;

    Map<const typename ParentUelDisplacement::RhsSized>                  QTotal(QTotal_);				
    Map<const typename ParentUelDisplacement::RhsSized>                  dQ(dQ_);			
    Map<typename ParentUelDisplacement::KeSizedMatrix>                   Ke(Ke_);
    Map<typename ParentUelDisplacement::RhsSized>                        Pe(Pe_);

    Ke.setZero();
    Pe.setZero();

    const int nStateVarsUmatPerGaussPt =        this->nStateVarsUmat;
    const int nStateVarsTotalPerGaussPt =       nStateVarsUmatPerGaussPt + 2*Vgt::VoigtSize;// stateVars per each integrationPoint	

    Vector6 dStrain;
    Matrix6 Cep;

    for(int i = 0; i < this->gaussWeights.size(); i++){

        const typename ParentUelDisplacement::ParentGeometryElement::BSized& B =   this->BAtGauss[i]; 
        const double detJ =         this->detJAtGauss[i];
        const double charElemlen =  std::sqrt(4*detJ);

        // mapping for global state vars; sequence: stateVarsLocal[nStateVarsUmatPerGaussPt], stress[nTens], strain[nTens]
        int GaussShiftStateVars =   this->nUelStatVars + i*nStateVarsTotalPerGaussPt;
        Ref<VectorXd> stateVarsUmat(this->stateVars.segment(GaussShiftStateVars,  nStateVarsUmatPerGaussPt));
        Ref<Vector6> stress( 		this->stateVars.segment(GaussShiftStateVars + nStateVarsUmatPerGaussPt,   Vgt::VoigtSize)  );
        Ref<Vector6> strain( 		this->stateVars.segment(GaussShiftStateVars + nStateVarsUmatPerGaussPt +  Vgt::VoigtSize, Vgt::VoigtSize) );	

        Cep.setZero();
        dStrain =       Vgt::planeVoigtToVoigt(B * dQ);
        const double vol =    detJ * this->propertiesElement[0] * this->gaussWeights(i);

        if (sectionType == SectionType::PlaneStress){

            bft::umatPlaneStress(Cep, stress, stateVarsUmat, strain, dStrain, this->propertiesUmat, 
                    pNewdT, charElemlen, time, dT, this->elLabel, i, this->umat);
            if (pNewdT<1.0)
                return;

            Ke += B.transpose() * mechanics::getPlaneStressTangent(Cep) * B * vol;
            Pe -= B.transpose() * Vgt::voigtToPlaneVoigt(stress)* vol;
        }

        else if(sectionType == SectionType::PlaneStrain)
        {
            bft::simpleUmat(Cep, stress, stateVarsUmat, strain, dStrain, this->propertiesUmat, 
                    pNewdT, charElemlen, time, dT, this->elLabel, i, this->umat);

            if (pNewdT<1.0)
                return; 

            Vector3d stressCondensed;
            stressCondensed.head(2) = stress.head(2);
            stressCondensed(2) = stress(3);

            Ke += B.transpose() * bft::mechanics::getPlaneStrainTangent(Cep) * B *  vol;
            Pe -= B.transpose() * stressCondensed * vol;
        }

        strain += dStrain; 

    }

}


template <int nNodes>
void UelDisplacementPlane<nNodes>::computeDistributedLoad( BftUel::DistributedLoadTypes loadType,
        double* P, 
        const int elementFace, 
        const double* load,
        const double* time,
        double dT){

    Map<typename ParentUelDisplacement::RhsSized> fU(P);

    using namespace bft::FiniteElement::BoundaryElementFactory;
    VectorXd boundaryCoordIndices = getBoundaryNodeList(this->shape, elementFace);
    
    VectorXd boundaryCoordinates(boundaryCoordIndices.size());
    for(int i = 0; i < boundaryCoordIndices.size(); i++)
        boundaryCoordinates(i) = this->coordinates( boundaryCoordIndices(i) );

    switch(loadType){

        case BftUel::Pressure: { 
            const double p = load[0];
            
            if (std::abs(p)<bft::Constants::numZeroPos)
                return;

            VectorXd Pk = VectorXd::Zero(boundaryCoordIndices.size());
            MatrixXd gp = getGaussPointList(this->shape); 
            VectorXd gpWeight = getGaussWeights(this->shape);  

            for(int i=0; i<gp.rows(); i++){
                MatrixXd xi = gp.row(i);        // necessary matrix mapping, as factory return type of gauss points is a matrix (with regard to future 3d elements) 
                VectorXd tractionVec = -p * getNormalVector(this->shape, boundaryCoordinates, xi);
                Pk += getIntVol(this->shape, boundaryCoordinates, xi) * this->propertiesElement(0) * gpWeight.row(i) * tractionVec.transpose() * getNB(this->shape, xi);}
            
            for(int i = 0; i < boundaryCoordIndices.size(); i++)
                fU( boundaryCoordIndices(i) ) +=  Pk(i);
            
            break;
        }
        default:  {std::cout << "Load type not supported" << std::endl; exit(-1);}
    }
}

