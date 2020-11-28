#pragma once
#include "Marmot/MarmotVoigt.h"

class HughesWinget {
  public:
    enum Formulation {
        AbaqusLike

    };

    HughesWinget( const Eigen::Matrix3d& FOld, const Eigen::Matrix3d& FNew, Formulation formulation );

    Marmot::Vector6d    getStrainIncrement();
    Eigen::Matrix3d getRotationIncrement();
    Marmot::Vector6d    rotateTensor( const Marmot::Vector6d& tensor );

    Marmot::Tensor633d compute_dS_dF( const Marmot::Vector6d&    stress,
                                   const Eigen::Matrix3d& FInv,
                                   const Marmot::Matrix6&    dChauchyDEps );
    Eigen::Matrix3d compute_dScalar_dF( const Eigen::Matrix3d& FInv, const Marmot::Vector6d& dScalarDEps );

  private:
    Formulation     theFormulation;
    Eigen::Matrix3d l;
    Eigen::Matrix3d dOmega;
    Eigen::Matrix3d dR;
    Marmot::Vector6d    dEps;
};
