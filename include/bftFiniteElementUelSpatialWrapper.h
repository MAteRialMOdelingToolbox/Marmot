#pragma once
#include "bftUel.h"
#include "Eigen/Sparse"
#include <functional>
#include <memory>

class BftUelSpatialWrapper : public BftUel
{
    /* Wrapper for Reduced Dimension Elements (e.g. Truss elements) to be used in higher order dimensions (2D, 3D).
     * The Projected is computed automatically based on the provided node coordinates, and the actual (child) element is created through a 
     * provided generator functor (e.g, function pointer)
     * */
    public:

        int projectedSize, unprojectedSize;

        std::unique_ptr<BftUel> childElement;
        Eigen::MatrixXd T;
        Eigen::MatrixXd P;

        BftUelSpatialWrapper(int nDim, int nChildDim, 
                const double* nodeCoordinates, 
                int nNodes, 
                int sizeRhsChild, 
                const int rhsIndicesToBeWrapped_[], int nRhsIndicesToBeWrapped,
                const std::function< BftUel* (const double *reducedNodeCoordinates)>& childGenerator);

        void computeYourself( const double* QTotal,
                                            const double* dQ,
                                            double* Pe,
                                            double* Ke,
                                            const double* time,
                                            double dT,
                                            double& pNewdT);

        void setInitialConditions(StateTypes state, const double* values);

        void computeDistributedLoad(
                                    DistributedLoadTypes loadType,
                                    double* P, 
                                    int elementFace, 
                                    const double* load,
                                    const double* time,
                                    double dT);

        double* getPermanentResultPointer(const std::string& resultName, int gaussPt, int& resultLength);
};
