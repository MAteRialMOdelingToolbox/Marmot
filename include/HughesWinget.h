#pragma once
#include "bftVoigt.h"

class HughesWinget {
  public:
    enum Formulation {
        AbaqusLike

    };

    HughesWinget( const Eigen::Matrix3d& FOld, const Eigen::Matrix3d& FNew, Formulation formulation );

    bft::Vector6    getStrainIncrement();
    Eigen::Matrix3d getRotationIncrement();
    bft::Vector6    rotateTensor( const bft::Vector6& tensor );

    bft::Tensor633d compute_dS_dF( const bft::Vector6&    stress,
                                   const Eigen::Matrix3d& FInv,
                                   const bft::Matrix6&    dChauchyDEps );
    Eigen::Matrix3d compute_dScalar_dF( const Eigen::Matrix3d& FInv, const bft::Vector6& dScalarDEps );

  private:
    Formulation     theFormulation;
    Eigen::Matrix3d l;
    Eigen::Matrix3d dOmega;
    Eigen::Matrix3d dR;
    bft::Vector6    dEps;
};
