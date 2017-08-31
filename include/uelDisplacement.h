#pragma once
#include "bftTypedefs.h"
#include "bftVoigt.h"
#include "bftConstants.h"
#include "bftUel.h"
#include "bftFiniteElement.h"
#include "bftFunctions.h"
#include "bftGeometryElement.h"
#include <iostream>
#include <vector>

template <int nDim, int nNodes>
class UelDisplacement: public BftUel, public BftGeometryElement<nDim, nNodes>{
    public:

        typedef BftGeometryElement<nDim, nNodes> ParentGeometryElement;

        static constexpr int sizeLoadVector =   nNodes * nDim;
        static constexpr int nCoordinates =     nNodes * nDim;
        static constexpr int nUelStatVars =     0;

        typedef Matrix<double, sizeLoadVector, 1>                RhsSized;
        typedef Matrix<double, sizeLoadVector, sizeLoadVector>   KeSizedMatrix;

        Map<VectorXd>               stateVars;
        const int                   nStateVars;
        const Map<const VectorXd>   propertiesElement;
        const Map<const VectorXd>   propertiesUmat;
        const int                   elLabel;
        const bft::pUmatType        umat;
        const int                   nStateVarsUmat;

        MatrixXd gaussPointList;
        VectorXd gaussWeights;

        std::vector< double >                                       detJAtGauss;        
        std::vector< typename ParentGeometryElement::JacobianSized> JAtGauss;
        std::vector< typename ParentGeometryElement::JacobianSized> JInvAtGauss;
        std::vector< typename ParentGeometryElement::dNdXiSized   > dNdXiAtGauss;
        std::vector< typename ParentGeometryElement::dNdXiSized   > dNdXAtGauss;
        std::vector< typename ParentGeometryElement::BSized       > BAtGauss;

    public:

        UelDisplacement(const double* coordinates,
                double* stateVars,
                int nStateVars,
                const double* propertiesElement,
                int nPropertiesElement,
                int noEl,
                const bft::pUmatType umat,
                int nStateVarsUmat,
                const double* propertiesUmat,
                int nPropertiesUmat,
                bft::NumIntegration::IntegrationTypes integrationType
                );

        virtual void setInitialConditions(StateTypes state, const double* values);

        virtual void computeDistributedLoad( BftUel::DistributedLoadTypes loadType,
                double* P, 
                const int elementFace, 
                const double* load,
                const double* time,
                double dT)
        {}

        virtual void computeYourself( const double* QTotal,
                const double* dQ,
                double* Pe,
                double* Ke,
                const double* time,
                double dT,
                double& pNewdT) = 0;
};

template <int nDim, int nNodes>
UelDisplacement<nDim, nNodes>::UelDisplacement(const double* coords, double* stateVars,
        int nStateVars,
        const double* properties,
        int nProperties,
        int noEl,
        const bft::pUmatType umat,
        int nStateVarsUmat ,
        const double* propertiesUmat,
        int nPropertiesUmat,
        bft::NumIntegration::IntegrationTypes integrationType
        ):
    ParentGeometryElement(coords),
    stateVars(Map<VectorXd>(stateVars, nStateVars)),
    nStateVars(nStateVars),
    propertiesElement(Map<const VectorXd>(properties, nProperties)),
    propertiesUmat(Map<const VectorXd>(propertiesUmat, nPropertiesUmat)),
    elLabel(noEl),
    umat(umat),
    nStateVarsUmat(nStateVarsUmat)
{
    gaussPointList =    bft::NumIntegration::getGaussPointList(this->shape, integrationType);
    gaussWeights =      bft::NumIntegration::getGaussWeights( this->shape, integrationType);

    for(int i = 0; i < gaussPointList.rows(); i++){
        const typename ParentGeometryElement::XiSized& xiGauss = gaussPointList.row(i);

        dNdXiAtGauss    .push_back		(this->dNdXi(xiGauss));
        JAtGauss		.push_back		(this->Jacobian(dNdXiAtGauss [i] ));
        JInvAtGauss		.push_back		(JAtGauss [i] .inverse());
        dNdXAtGauss		.push_back		(this->dNdX(dNdXiAtGauss [i] , JInvAtGauss [i] ));
        BAtGauss		.push_back		(this->B(dNdXAtGauss [i] ));
        detJAtGauss		.push_back		(JAtGauss [i] .determinant());
    }
}

template <int nDim, int nNodes>
void UelDisplacement<nDim, nNodes>::setInitialConditions(StateTypes state, const double* values)
{
    const int nStateVarsTotalPerGaussPt = nStateVarsUmat + 2*bft::Vgt::VoigtSize;    

    switch(state){

        case BftUel::GeostaticStress: 
            { 
                for(int i = 0; i < this->gaussPointList.rows(); i++) 
                {
                    typename ParentGeometryElement::XiSized coordAtGauss =  this->NB( this->N(this->gaussPointList.row(i))) * this->coordinates; 

                    int GaussShiftStateVars =   nUelStatVars + i*nStateVarsTotalPerGaussPt;
                    Ref<bft::Vector6> stress(stateVars.segment(GaussShiftStateVars + nStateVarsUmat, bft::Vgt::VoigtSize));

                    const double sigY1 = values[0];
                    const double sigY2 = values[2];
                    const double y1    = values[1];
                    const double y2    = values[3];
                    
                    stress(1) = bft::Functions::linearInterpolation(coordAtGauss[1], y1, y2, sigY1, sigY2);  // sigma_y
                    stress(0) = values[4]*stress(1);  // sigma_x
                    stress(2) = values[5]*stress(1);}  // sigma_z

                    break; 
            }

        default: break;
    }
}

