#pragma once
#include "bftTypedefs.h"

namespace bft{
    template <int nMatTangentSize>
    class DuvautLionsViscosity
    {
        private:

            const double viscosity;
        public:

            using TangentSizedMatrix= Matrix<double, nMatTangentSize, nMatTangentSize>;
            DuvautLionsViscosity(double viscosity);
            double applyViscosityOnStateVar(double stateVarTrial, double StateVarInf, double dT);
            Vector6 applyViscosityOnStress(const Vector6& trialStress, const Vector6& stressInf, double dT);
            TangentSizedMatrix applyViscosityOnMatTangent(
                    const TangentSizedMatrix& matTangentInv, 
                    double dT); 
    };
}

namespace bft
{
    template <int s>
    DuvautLionsViscosity<s>::DuvautLionsViscosity(double viscosity):
            viscosity(viscosity)
    {

    }

    template <int s>
    double DuvautLionsViscosity<s>::applyViscosityOnStateVar(double stateVarTrial, double StateVarInf,  double dT)
    {
            return (stateVarTrial + (dT/viscosity) * StateVarInf)/(dT/viscosity + 1); 
    }

    template <int s>
    Vector6 DuvautLionsViscosity<s>::applyViscosityOnStress(const Vector6& trialStress, const Vector6& stressInf, double dT)
    {
            return (trialStress + (dT/viscosity) * stressInf)/(dT/viscosity + 1); 
    }
            
    template <int s>
    typename DuvautLionsViscosity<s>::TangentSizedMatrix 
        DuvautLionsViscosity<s>::applyViscosityOnMatTangent(const TangentSizedMatrix& matTangentInv, double dT)
    {
            return (1/(1+dT/viscosity)) * (TangentSizedMatrix::Identity() + dT/viscosity * matTangentInv);
    }
}
