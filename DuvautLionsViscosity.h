#pragma once
#include "Marmot/MarmotTypedefs.h"

namespace Marmot {
    template <int nMatTangentSize>
    class DuvautLionsViscosity {
      private:
        const double viscosity;

      public:
        typedef Eigen::Matrix<double, nMatTangentSize, nMatTangentSize> TangentSizedMatrix;

        DuvautLionsViscosity( double viscosity );
        double             applyViscosityOnStateVar( double stateVarTrial, double StateVarInf, double dT );
        Marmot::Vector6d       applyViscosityOnStress( const Marmot::Vector6d& trialStress,
                                                   const Marmot::Vector6d& stressInf,
                                                   double              dT );
        TangentSizedMatrix applyViscosityOnMatTangent( const TangentSizedMatrix& matTangentInv, double dT );
    };
} // namespace Marmot

namespace Marmot {
    template <int s>
    DuvautLionsViscosity<s>::DuvautLionsViscosity( double viscosity ) : viscosity( viscosity )
    {
    }

    template <int s>
    double DuvautLionsViscosity<s>::applyViscosityOnStateVar( double stateVarTrial, double StateVarInf, double dT )
    {
        return ( stateVarTrial + ( dT / viscosity ) * StateVarInf ) / ( dT / viscosity + 1 );
    }

    template <int s>
    Marmot::Vector6d DuvautLionsViscosity<s>::applyViscosityOnStress( const Marmot::Vector6d& trialStress,
                                                                  const Marmot::Vector6d& stressInf,
                                                                  double              dT )
    {
        return ( trialStress + ( dT / viscosity ) * stressInf ) / ( dT / viscosity + 1 );
    }

    template <int s>
    typename DuvautLionsViscosity<s>::TangentSizedMatrix DuvautLionsViscosity<s>::applyViscosityOnMatTangent(
        const TangentSizedMatrix& matTangentInv,
        double                    dT )
    {
        return ( 1 / ( 1 + dT / viscosity ) ) * ( TangentSizedMatrix::Identity() + dT / viscosity * matTangentInv );
    }
} // namespace Marmot
