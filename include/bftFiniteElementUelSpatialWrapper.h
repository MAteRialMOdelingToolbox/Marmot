#pragma once
#include "Eigen/Sparse"
#include "bftUel.h"
#include "bftUelProperty.h"
#include <functional>
#include <memory>

class BftUelSpatialWrapper : public BftUel {
    /* Wrapper for Reduced Dimension Elements (e.g. Truss elements) to be used in higher order
     * dimensions (2D, 3D). The Projected is computed automatically based on the provided node
     * coordinates, and the actual (child) element is created through a provided generator functor
     * (e.g, function pointer)
     * */

  public:
    const int                               nDim;
    const int                               nDimChild;
    const int                               nNodes;
    const int                               nRhsChild;
    const Eigen::Map<const Eigen::VectorXi> rhsIndicesToBeProjected;
    const int                               projectedSize, unprojectedSize;

    std::unique_ptr<BftUel> childElement;
    Eigen::MatrixXd         T;
    Eigen::MatrixXd         P;
    Eigen::MatrixXd         projectedCoordinates;

    BftUelSpatialWrapper( int                     nDim,
                          int                     nChildDim,
                          int                     nNodes,
                          int                     sizeRhsChild,
                          const int               rhsIndicesToBeWrapped_[],
                          int                     nRhsIndicesToBeWrapped,
                          std::unique_ptr<BftUel> childElement );

    int getNumberOfRequiredStateVars();

    std::vector<std::vector<std::string>> getNodeFields();

    std::vector<int> getDofIndicesPermutationPattern();

    int getNNodes();

    int getNDofPerElement();

    std::string getElementShape();

    void assignStateVars( double* stateVars, int nStateVars );

    void assignProperty( const BftUelProperty& property );

    void initializeYourself( const double* coordinates );

    void computeYourself( const double* QTotal,
                          const double* dQ,
                          double*       Pe,
                          double*       Ke,
                          const double* time,
                          double        dT,
                          double&       pNewdT );

    void setInitialConditions( StateTypes state, const double* values );

    void computeDistributedLoad( DistributedLoadTypes loadType,
                                 double*              P,
                                 int                  elementFace,
                                 const double*        load,
                                 const double*        time,
                                 double               dT );

    void computeBodyForce( double*              P,
                                 const double*        load,
                                 const double*        time,
                                 double               dT );

    double* getPermanentResultPointer( const std::string& resultName, int gaussPt, int& resultLength );
};
