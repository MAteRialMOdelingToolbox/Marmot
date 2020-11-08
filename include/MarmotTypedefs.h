#pragma once
#include "Eigen/Core"
#include "Eigen/Dense"
#include "unsupported/Eigen/CXX11/Tensor"

namespace marmot {
    typedef Eigen::Matrix<double, 6, 6> Matrix6;
    typedef Eigen::Matrix<double, 6, 9> Matrix69d;
    typedef Eigen::Matrix<double, 9, 9> Matrix99d;
    typedef Eigen::Matrix<double, 3, 4> Matrix34d;
    typedef Eigen::Matrix<double, 7, 7> Matrix7;
    typedef Eigen::Map<Matrix6>         mMatrix6;

    typedef Eigen::Matrix<double, 6, 1>    Vector6;
    typedef Eigen::Matrix<double, 7, 1>    Vector7d;
    typedef Eigen::Matrix<double, 8, 1>    Vector8d;
    typedef Eigen::Matrix<double, 9, 1>    Vector9d;
    typedef Eigen::Matrix<int, 8, 1>       Vector8i;
    typedef Eigen::Matrix<double, 1, 6>    RowVector6d;
    typedef Eigen::Map<Vector6>            mVector6;
    typedef Eigen::Map<Eigen::VectorXd>    mVectorXd;
    typedef Eigen::Map<const marmot::Vector6> mConstVector6;

    typedef Eigen::Matrix<double, 3, 6> Matrix36d;
    typedef Eigen::Matrix<double, 3, 6> Matrix36;
    typedef Eigen::Matrix<double, 9, 9> Matrix9d;

    typedef Eigen::TensorFixedSize<double, Eigen::Sizes<6,3,3>> Tensor633d; 
    typedef Eigen::TensorFixedSize<double, Eigen::Sizes<3,2,2>> Tensor322d; 
    typedef Eigen::TensorFixedSize<double, Eigen::Sizes<3,3,3,3>> Tensor3333d; 
    typedef Eigen::TensorFixedSize<double, Eigen::Sizes<3, 3, 3>> Tensor333d;
    typedef Eigen::TensorFixedSize<double, Eigen::Sizes<1, 2, 2>> Tensor122d;
    typedef Eigen::TensorFixedSize<double, Eigen::Sizes<2, 2, 2, 2>> Tensor2222d;
    typedef Eigen::TensorFixedSize<double, Eigen::Sizes<2, 2, 1, 2>> Tensor2212d;
    typedef Eigen::TensorFixedSize<double, Eigen::Sizes<2, 1, 2, 2>> Tensor2122d;
    typedef Eigen::TensorFixedSize<double, Eigen::Sizes<2, 1, 1, 2>> Tensor2112d;

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

} // namespace marmot
