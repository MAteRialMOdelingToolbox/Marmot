#pragma once
#include "Eigen/Core"
#include "Eigen/Dense"

using namespace Eigen;

namespace bft{
    typedef Matrix<double,6,6> Matrix6;
	typedef Matrix<double,3,6> Matrix36;
    typedef Map<Matrix6> mMatrix6;

    typedef Matrix<double,6,1> Vector6;
	typedef Matrix<double,1,6> RowVector6d;
    typedef Map<Vector6> mVector6;
    typedef Map<const Vector6> mConstVector6;
}
