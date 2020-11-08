#pragma once
#include "MarmotVoigt.h"

class HughesWinget {
  public:
    enum Formulation {
        AbaqusLike

    };

    HughesWinget( const Eigen::Matrix3d& FOld, const Eigen::Matrix3d& FNew, Formulation formulation );

    marmot::Vector6    getStrainIncrement();
    Eigen::Matrix3d getRotationIncrement();
    marmot::Vector6    rotateTensor( const marmot::Vector6& tensor );

    marmot::Tensor633d compute_dS_dF( const marmot::Vector6&    stress,
                                   const Eigen::Matrix3d& FInv,
                                   const marmot::Matrix6&    dChauchyDEps );
    Eigen::Matrix3d compute_dScalar_dF( const Eigen::Matrix3d& FInv, const marmot::Vector6& dScalarDEps );

  private:
    Formulation     theFormulation;
    Eigen::Matrix3d l;
    Eigen::Matrix3d dOmega;
    Eigen::Matrix3d dR;
    marmot::Vector6    dEps;
};
