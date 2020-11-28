#include "Marmot/MarmotAbaqusUtility.h"
#include "Marmot/MarmotFunctions.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"
/*
 * convenience functions for umats in Abaqus
 * */

using namespace Eigen;

namespace Marmot {

    void backToAbaqus( const Matrix6& jacobian,
                       Map<MatrixXd>& ABQJacobian,
                       const Vector6d& stress,
                       Map<VectorXd>& ABQStress,
                       int            nTensor )
    {
        ABQStress   = stress.head( nTensor );
        ABQJacobian = jacobian.topLeftCorner( nTensor, nTensor );
        return;
    }

    void backToAbaqusPlaneStress( const Matrix6& jacobian,
                                  Map<MatrixXd>& ABQJacobian,
                                  const Vector6d& stress,
                                  Map<VectorXd>& ABQStress )
    {
        ABQStress( 0 ) = stress( 0 );
        ABQStress( 1 ) = stress( 1 );
        ABQStress( 2 ) = stress( 3 );
        ABQJacobian    = Marmot::mechanics::getPlaneStressTangent( jacobian );
        return;
    }

    void backToAbaqusNonLocal( const Matrix6& dStressdStrain,
                               Ref<MatrixXd>  ABQdStressDStrain,
                               const Vector6d& stress,
                               Ref<VectorXd>  ABQStress,
                               double         intParameterLocal,
                               double&        ABQParameterLocal,
                               const Vector6d& dStressDIntParamNonLocal,
                               Ref<VectorXd>  ABQDStressDIntParamNonLocal,
                               const Vector6d& dIntParamLocalDStrain,
                               Ref<VectorXd>  ABQDIntParameterLocalDStrain,
                               double         nonLocalRadius,
                               double&        ABQNonLocalRadius,
                               int            nTensor )
    {
        ABQdStressDStrain            = dStressdStrain.topLeftCorner( nTensor, nTensor );
        ABQStress                    = stress.head( nTensor );
        ABQParameterLocal            = intParameterLocal;
        ABQDStressDIntParamNonLocal  = dStressDIntParamNonLocal.head( nTensor );
        ABQDIntParameterLocalDStrain = dIntParamLocalDStrain.head( nTensor );
        ABQNonLocalRadius            = nonLocalRadius;
    }

    void discardTheIncrement( double& pNewDT, double value, const std::string& message )
    {
        pNewDT = value;
        MarmotJournal::warningToMSG( message );
        return;
    }

} // namespace Marmot
