#pragma once
#include "Eigen/Core"
#include "Eigen/Dense"

using namespace Eigen;

namespace bft{
    typedef Matrix<double,6,6> Matrix6;
    typedef Map<Matrix6> mMatrix6;

    typedef Matrix<double,6,1> Vector6;
    typedef Map<Vector6> mVector6;
    typedef Map<const Vector6> mConstVector6;

    template <std::size_t size> 
        using YieldSurfFlagArr = Array<bool, 1, size>;
    template <std::size_t size>
        using YieldSurfResArr = Array<double, 1, size>;
};
