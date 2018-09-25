#include "bftFiniteElementUelSpatialWrapper.h"
#include "bftVoigt.h"
#include "Eigen/Sparse"
#include <iostream>

BftUelSpatialWrapper::BftUelSpatialWrapper(
        int nDim, int nDimChild, 
        int nNodes, 
        int nRhsChild, 
        const int rhsIndicesToBeWrapped[], 
        int nRhsIndicesToBeWrapped,
        std::unique_ptr<BftUel> childElement) :
    nDim(nDim),
    nDimChild(nDimChild),
    nNodes(nNodes),
    nRhsChild(nRhsChild),
    rhsIndicesToBeProjected (rhsIndicesToBeWrapped, nRhsIndicesToBeWrapped),
    projectedSize (nRhsChild),
    unprojectedSize (nRhsChild + rhsIndicesToBeProjected.size() * (nDim - nDimChild)),
    childElement ( std::move( childElement ) )
{
}

int BftUelSpatialWrapper::getNumberOfRequiredStateVars(){
    return childElement->getNumberOfRequiredStateVars();
}

std::vector< std::vector<std::string>> BftUelSpatialWrapper::getNodeFields()
{
    return childElement->getNodeFields();
}
int BftUelSpatialWrapper::getNNodes(){
    return nNodes;}

int BftUelSpatialWrapper::getNDofPerElement (){
    return unprojectedSize;}

std::string BftUelSpatialWrapper::getElementShape(){
    return childElement->getElementShape();}

void BftUelSpatialWrapper::assignProperty(BftUel::PropertyTypes property, int propertyInfo, const double* propertyValues, int nProperties) {
    childElement->assignProperty(property, propertyInfo, propertyValues, nProperties);}

std::vector<int> BftUelSpatialWrapper::getDofIndicesPermutationPattern()
{
    auto childPermutationPattern = childElement->getDofIndicesPermutationPattern();
    std::vector<int> permutationPattern;

    int i=0, j=0, k=0; //i: projected, j: unprojected, k: indicesToProject
    int indexChild;
    while( i < projectedSize)
    {
        indexChild = childPermutationPattern[i];
        if ( k < rhsIndicesToBeProjected.size() && i == rhsIndicesToBeProjected(k) )
        {
            // copy the Transformation block wise ( if current DOF is projected )
            for(int l = 0; l < nDim; l ++)
                permutationPattern.push_back( indexChild + j +  l );

            i += nDimChild;
            j += (nDim-nDimChild);
            k += 1;
        }
        else
        { 
            permutationPattern.push_back( indexChild + j );
            i++;
            j++;
        }
    }
   
    return permutationPattern;
}


void BftUelSpatialWrapper::assignStateVars(double *stateVars, int nStateVars)
{
    childElement->assignStateVars(stateVars, nStateVars);
}

void BftUelSpatialWrapper::initializeYourself(const double *coordinates) {

    Map < const MatrixXd > unprojectedCoordinates (coordinates, nDim, nNodes);


    if(nDimChild == 1)
    {
        // take the (normalized) directional vector based on the first 2 nodes (valid for truss2, truss3);
        VectorXd n = unprojectedCoordinates.col(1) - unprojectedCoordinates.col(0);
        n /= n.norm();
        T = n.transpose ();
    }

    P = MatrixXd::Zero(projectedSize, unprojectedSize);
    // Alternatively: P as a Sparse matrix; currently not necessary...

    int i=0, j=0, k=0; //i: projected, j: unprojected, k: indicesToProject
    while( i < projectedSize)
    {
        if ( k < rhsIndicesToBeProjected.size() && i == rhsIndicesToBeProjected(k) )
        {
            // copy the Transformation block wise ( if current DOF is projected )
            P.block(i, j, nDimChild, nDim) = T;

            i += nDimChild;
            j += nDim;
            k += 1;
        }
        else
        { 
            P(i, j) = 1.0;
            i++;
            j++;
        }
    }

    // Projection of node coordinates
    projectedCoordinates = MatrixXd::Zero(nDimChild, nNodes);
    for (int i = 0; i < nNodes; i++)
        projectedCoordinates.col(i) = T * unprojectedCoordinates.col(i);

    childElement->initializeYourself(projectedCoordinates.data());
}

void BftUelSpatialWrapper::computeYourself( const double* Q,
        const double* dQ,
        double* Pe_,
        double* Ke_,
        const double* time,
        double dT,
        double& pNewDT) 
{
    Map<const VectorXd> Q_Unprojected (  Q, unprojectedSize);
    Map<const VectorXd> dQ_Unprojected ( dQ, unprojectedSize);

    VectorXd  Q_Projected = P * Q_Unprojected;
    VectorXd dQ_Projected = P * dQ_Unprojected;

    VectorXd Pe_Projected(projectedSize);
    MatrixXd Ke_Projected(projectedSize, projectedSize);

    childElement->computeYourself(Q_Projected.data(),
            dQ_Projected.data(),
            Pe_Projected.data(),
            Ke_Projected.data(),
            time, dT, pNewDT);

    if(pNewDT < 1.0)
        return;

    Map<VectorXd> Pe_Unprojected(Pe_, unprojectedSize);
    Map<MatrixXd> Ke_Unprojected(Ke_, unprojectedSize, unprojectedSize);

    Ke_Unprojected = P.transpose() * Ke_Projected * P ;
    Pe_Unprojected = P.transpose() * Pe_Projected;
}

void BftUelSpatialWrapper::setInitialConditions(StateTypes state, const double* values)
{
    childElement->setInitialConditions(state, values); 
}

void BftUelSpatialWrapper::computeDistributedLoad(
        DistributedLoadTypes loadType,
        double* P_, 
        int elementFace, 
        const double* load,
        const double* time,
        double dT) 
{
    VectorXd P_Projected = VectorXd::Zero(projectedSize);

    childElement->computeDistributedLoad(loadType, P_Projected.data(), elementFace, load, time, dT);

    Map<VectorXd> P_Unprojected(P_, unprojectedSize);
    P_Unprojected = P.transpose() * P_Projected;

}

double* BftUelSpatialWrapper::getPermanentResultPointer(const std::string& resultName, int gaussPt, int& resultLength) 
{
    if (resultName == "BftUelSpatialWrapper.T")
    {
        resultLength = T.size();
        return T.data();
    }
    else
    {
        return childElement->getPermanentResultPointer(resultName, gaussPt, resultLength);
    }
}
