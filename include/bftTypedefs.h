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

    template <int size> 
		struct YieldSurfFlagArr
		{
			typedef Array<bool, 1, size> nSurfs;
		};
		
    template <int size>
		struct YieldSurfResArr
		{
			typedef Array<double, 1, size> nSurfs;
		};
};
