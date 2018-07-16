#include "bftFiniteElementUelSpatialWrapper.h"
#include "bftVoigt.h"
#include "Eigen/Sparse"
#include <iostream>

BftUelSpatialWrapper::BftUelSpatialWrapper(int nDim, int nDimChild, 
        const double* nodeCoordinates, 
        int nNodes, 
        int nRhsChild, 
        const int rhsIndicesToBeWrapped[], int nRhsIndicesToBeWrapped,
        const std::function< BftUel* (const double *reducedNodeCoordinates)>& childGenerator)
{
    Map < const MatrixXd > unprojectedCoordinates (nodeCoordinates, nDim, nNodes);
    Map < const VectorXi > rhsIndicesToBeProjected (rhsIndicesToBeWrapped, nRhsIndicesToBeWrapped);

    int nDimDiff = nDim - nDimChild;

    projectedSize =   nRhsChild;
    unprojectedSize = nRhsChild + nRhsIndicesToBeWrapped * nDimDiff;

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
    MatrixXd projectedCoordinates(nDimChild, nNodes);
    for (int i = 0; i < nNodes; i++)
        projectedCoordinates.col(i) = T * unprojectedCoordinates.col(i);

    childElement = std::unique_ptr<BftUel> ( childGenerator ( projectedCoordinates.data() ) );
}

void BftUelSpatialWrapper::computeYourself( const double* Q,
        const double* dQ,
        double* Pe_,
        double* Ke_,
        const double* time,
        double dT,
        double& pNewDT) 
{
    VectorXd  Q_Projected = P * Map<const VectorXd>(  Q, unprojectedSize);
    VectorXd dQ_Projected = P * Map<const VectorXd>( dQ, unprojectedSize);

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
