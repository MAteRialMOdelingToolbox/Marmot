#pragma once
#include "Eigen/Core"
#include "Eigen/Dense"

// using namespace Eigen;

namespace bft {
    typedef Eigen::Matrix<double, 6, 6> Matrix6;
    typedef Eigen::Matrix<double, 7, 7> Matrix7;
    typedef Eigen::Map<Matrix6>         mMatrix6;

    typedef Eigen::Matrix<double, 6, 1>    Vector6;
    typedef Eigen::Matrix<double, 7, 1>    Vector7d;
    typedef Eigen::Matrix<double, 8, 1>    Vector8d;
    typedef Eigen::Matrix<int, 8, 1>       Vector8i;
    typedef Eigen::Matrix<double, 1, 6>    RowVector6d;
    typedef Eigen::Map<Vector6>            mVector6;
    typedef Eigen::Map<Eigen::VectorXd>    mVectorXd;
    typedef Eigen::Map<const bft::Vector6> mConstVector6;

    typedef Eigen::Matrix<double, 3, 6> Matrix36;

    typedef void ( *pUmatType )( double[],
                                 double[],
                                 double[],
                                 double&,
                                 double&,
                                 double&,
                                 double&,
                                 double[],
                                 double[],
                                 double&,
                                 const double[],
                                 const double[],
                                 const double[2],
                                 const double&,
                                 const double&,
                                 const double&,
                                 const double[],
                                 const double[],
                                 const char[80],
                                 const int&,
                                 const int&,
                                 const int&,
                                 const int&,
                                 const double[],
                                 const int&,
                                 const double[3],
                                 const double[9],
                                 double&,
                                 const double&,
                                 const double[9],
                                 const double[9],
                                 const int&,
                                 const int&,
                                 const int&,
                                 const int&,
                                 const int[4],
                                 const int&,
                                 const int );

} // namespace bft
