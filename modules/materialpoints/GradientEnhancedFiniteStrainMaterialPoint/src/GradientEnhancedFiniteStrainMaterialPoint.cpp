#include "Marmot/GradientEnhancedFiniteStrainMaterialPoint.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotUtils.h"
#include <stdexcept>

namespace Marmot::MaterialPoints {

  // =========================================================================================
  //  plane strain
  // =========================================================================================

  void GradientEnhancedFiniteStrainMaterialPoint2D::incrementDeformation( const TensorD&  du,
                                                                          const TensorDD& du_dY,
                                                                          double          dn )
  {
    mapEigenToFastor( state->du ) += mapEigenToFastor( expandTo3D( du ) );
    mapEigenToFastor( state->dx_dY ) += mapEigenToFastor( expandTo3D( du_dY ) );
    state->nonLocalDamage += dn;
  };

  void GradientEnhancedFiniteStrainMaterialPoint2D::computeYourself( double timeNew, double dT )
  {
    const Material::TimeIncrement timeIncrement{ timeNew, dT };

    // total deformation gradient: the increment since the last accepted state, composed onto it
    const FastorStandardTensors::Tensor33d dx_dX = state->dx_dY % state->dY_dX;

    Material::ConstitutiveResponse< 3 > response3D{ FastorStandardTensors::Tensor33d( 0.0 ), 0.0, 0.0, 0.0, 0.0,
                                                    state->materialState.data() };
    Material::AlgorithmicModuli< 3 >    algorithmicModuli3D;
    Material::Deformation< 3 >          deformation3D{ dx_dX, state->nonLocalDamage };

    if ( hasEigenDeformation )
      material->computePlaneStrain( response3D, algorithmicModuli3D, deformation3D, timeIncrement,
                                    { state->F0_XX, state->F0_YY, state->F0_ZZ } );
    else
      material->computePlaneStrain( response3D, algorithmicModuli3D, deformation3D, timeIncrement );

    // clang-format off
    this->response = { .S              = reduceTo2D< U, U >( response3D.tau ),
                       .dL             = response3D.L - state->localDamage,
                       .nonLocalRadius = response3D.nonLocalRadius };
    // clang-format on

    state->localDamage = response3D.L;
    state->stress      = response3D.tau;

    using namespace FastorStandardTensors;
    using namespace FastorIndices;

    const auto& I = Marmot::FastorStandardTensors::Spatial3D::I;

    // chain rule from d/dF to d/d(DeltaF), since the unknown of the increment is DeltaF = dx_dY
    const Tensor3333d dF_dDeltaF_3D = einsum< ij, JI, to_iIjJ >( I, state->dY_dX );
    const TensorDDDD  dF_dDeltaF    = reduceTo2D< U, U, U, U >( dF_dDeltaF_3D );

    this->tangents = {
      reduceTo2D< U, U, U, U >( algorithmicModuli3D.dTau_dF ),
      reduceTo2D< U, U >( algorithmicModuli3D.dTau_dN ),
      reduceTo2D< U, U >( algorithmicModuli3D.dL_dF ),
      algorithmicModuli3D.dL_dN,
    };

    tangents.dS_dDeltaF = einsum< ijmn, mnKL >( tangents.dS_dDeltaF, dF_dDeltaF );
    tangents.dL_dDeltaF = einsum< mn, mnKL >( tangents.dL_dDeltaF, dF_dDeltaF );

    _density = material->getDensity( state->materialState.data() );
  }

  // =========================================================================================
  //  3D
  // =========================================================================================

  void GradientEnhancedFiniteStrainMaterialPoint3D::incrementDeformation( const TensorD&  du,
                                                                          const TensorDD& du_dY,
                                                                          double          dn )
  {
    mapEigenToFastor( state->du ) += mapEigenToFastor( du );
    mapEigenToFastor( state->dx_dY ) += mapEigenToFastor( du_dY );
    state->nonLocalDamage += dn;
  };

  void GradientEnhancedFiniteStrainMaterialPoint3D::computeYourself( double timeNew, double dT )
  {
    const Material::TimeIncrement timeIncrement{ timeNew, dT };

    const FastorStandardTensors::Tensor33d dx_dX = state->dx_dY % state->dY_dX;

    Material::ConstitutiveResponse< 3 > response3D{ FastorStandardTensors::Tensor33d( 0.0 ), 0.0, 0.0, 0.0, 0.0,
                                                    state->materialState.data() };
    Material::AlgorithmicModuli< 3 >    algorithmicModuli3D;
    Material::Deformation< 3 >          deformation3D{ dx_dX, state->nonLocalDamage };

    if ( hasEigenDeformation )
      material->computeStress( response3D, algorithmicModuli3D, deformation3D, timeIncrement,
                               { state->F0_XX, state->F0_YY, state->F0_ZZ } );
    else
      material->computeStress( response3D, algorithmicModuli3D, deformation3D, timeIncrement );

    // clang-format off
    this->response = { .S              = ( response3D.tau ),
                       .dL             = response3D.L - state->localDamage,
                       .nonLocalRadius = response3D.nonLocalRadius };
    // clang-format on

    state->localDamage = response3D.L;
    state->stress      = response3D.tau;

    using namespace FastorStandardTensors;
    using namespace FastorIndices;

    const auto& I = Marmot::FastorStandardTensors::Spatial3D::I;

    const Tensor3333d dF_dDeltaF = einsum< ij, JI, to_iIjJ >( I, state->dY_dX );

    this->tangents = {
      ( algorithmicModuli3D.dTau_dF ),
      ( algorithmicModuli3D.dTau_dN ),
      ( algorithmicModuli3D.dL_dF ),
      algorithmicModuli3D.dL_dN,
    };

    tangents.dS_dDeltaF = einsum< ijmn, mnKL >( tangents.dS_dDeltaF, dF_dDeltaF );
    tangents.dL_dDeltaF = einsum< mn, mnKL >( tangents.dL_dDeltaF, dF_dDeltaF );

    _density = material->getDensity( state->materialState.data() );
  }

} // namespace Marmot::MaterialPoints
