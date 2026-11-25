/* ---------------------------------------------------------------------
 *                                       _
 *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_
 * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
 * | | | | | | (_| | |  | | | | | | (_) | |_
 * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
 *
 * Unit of Strength of Materials and Structural Analysis
 * University of Innsbruck
 * 2020 - today
 *
 * festigkeitslehre@uibk.ac.at
 *
 * This file is part of the MAteRialMOdellingToolbox (marmot).
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of marmot.
 * ---------------------------------------------------------------------
 */
#include "Marmot/GradientEnhancedFiniteStrainMaterialPoint.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotUtils.h"
#include <stdexcept>

namespace Marmot::Discretizations::MaterialPoints {

  template < int nDim, int nNonlocal >
  void GradientEnhancedFiniteStrainMaterialPoint< nDim, nNonlocal >::assignStateVars( double* stateVars,
                                                                                      int     nStateVars )
  {
    state                    = std::make_unique< MPStateVarManager >( stateVars, nStateVars );
    this->_assignedStateVars = { stateVars, nStateVars };
  }

  template < int nDim, int nNonlocal >
  Marmot::Utilities::StateView GradientEnhancedFiniteStrainMaterialPoint< nDim, nNonlocal >::getStateView(
    const std::string& stateName ) const
  {
    if ( state->contains( stateName ) )
      return state->getStateView( stateName );
    else
      return material->getStateView( stateName, state->materialState.data() );
  }

  template < int nDim, int nNonlocal >
  void GradientEnhancedFiniteStrainMaterialPoint< nDim, nNonlocal >::assignMaterial(
    const Marmot::Discretizations::MaterialSection& section )
  {
    material.reset( Marmot::Registry::MaterialGradientEnhancedFiniteStrainFactory<
                    nNonlocal >::createMaterial( section.materialName,
                                                 section.materialProperties,
                                                 section.nMaterialProperties,
                                                 _mpNumber ) );
  }

  template < int nDim, int nNonlocal >
  void GradientEnhancedFiniteStrainMaterialPoint< nDim, nNonlocal >::initializeYourself()
  {
    state->dY_dX.eye();
    state->F0_XX = 1.0;
    state->F0_YY = 1.0;
    state->F0_ZZ = 1.0;
    this->prepareYourself( 0, 0 );
    material->initializeYourself( state->materialState.data(), static_cast< int >( state->materialState.size() ) );
    // known from the start, so that inertia can be assembled before the first computation
    _density = material->getDensity( state->materialState.data() );
  }

  template < int nDim, int nNonlocal >
  void GradientEnhancedFiniteStrainMaterialPoint< nDim, nNonlocal >::prepareYourself( double timeNew, double dT )
  {
    state->du.zeros();
    state->dx_dY.eye();
  }

  template < int nDim, int nNonlocal >
  void GradientEnhancedFiniteStrainMaterialPoint< nDim, nNonlocal >::incrementDeformation(
    const TensorD&                          du,
    const TensorDD&                         du_dY,
    const NonlocalVector&                   dn,
    const std::array< TensorD, nNonlocal >& dn_dY )
  {
    if constexpr ( nDim == 2 ) {
      Marmot::Math::Tensors::mapEigenToFastor( state->du ) += Marmot::Math::Tensors::mapEigenToFastor(
        Marmot::Math::Tensors::expandTo3D( du ) );
      Marmot::Math::Tensors::mapEigenToFastor( state->dx_dY ) += Marmot::Math::Tensors::mapEigenToFastor(
        Marmot::Math::Tensors::expandTo3D( du_dY ) );
    }
    else {
      Marmot::Math::Tensors::mapEigenToFastor( state->du ) += Marmot::Math::Tensors::mapEigenToFastor( du );
      Marmot::Math::Tensors::mapEigenToFastor( state->dx_dY ) += Marmot::Math::Tensors::mapEigenToFastor( du_dY );
    }
    state->nonlocalField += dn;

    // grad_X dn = F_n^T grad_Y dn, with the deformation gradient of the last accepted state
    for ( int K = 0; K < nNonlocal; K++ )
      for ( int i = 0; i < nDim; i++ )
        for ( int j = 0; j < nDim; j++ )
          state->nonlocalFieldGradient( i, K ) += state->dY_dX( j, i ) * dn_dY[K]( j );
  }

  template < int nDim, int nNonlocal >
  Marmot::Math::Tensors::FastorTensors::Tensor33d GradientEnhancedFiniteStrainMaterialPoint< nDim, nNonlocal >::
    deformationGradient() const
  {
#if defined( __GNUC__ ) && !defined( __clang__ )
// GCC's -Wmaybe-uninitialized misfires on Fastor's AVX copies of the mapped state tensors: the
// mapped memory is initialized by the host, which GCC cannot see (observed with GCC 14 at x86-64-v3).
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
    const Marmot::Math::Tensors::FastorTensors::Tensor33d dx_dX = state->dx_dY % state->dY_dX;
#if defined( __GNUC__ ) && !defined( __clang__ )
#  pragma GCC diagnostic pop
#endif
    return dx_dX;
  }

  template < int nDim, int nNonlocal >
  void GradientEnhancedFiniteStrainMaterialPoint< nDim, nNonlocal >::computeYourself( double timeNew, double dT )
  {
    using namespace Marmot::Math::Tensors::FastorTensors;
    using namespace Marmot::Math::Tensors::FastorIndices;
    using Marmot::Math::Tensors::reduceTo2D;
    using Marmot::Math::Tensors::U;

    const typename Material::TimeIncrement timeIncrement{ timeNew, dT };

    // total deformation gradient: the increment since the last accepted state, composed onto it
    const Tensor33d dx_dX = deformationGradient();

    typename Material::template ConstitutiveResponse< 3 > response3D{ Tensor33d( 0.0 ),
                                                                      NonlocalVector::Zero(),
                                                                      NonlocalVector::Zero(),
                                                                      0.0,
                                                                      0.0,
                                                                      state->materialState.data() };
    typename Material::template AlgorithmicModuli< 3 >    algorithmicModuli3D;
    typename Material::template Deformation< 3 >          deformation3D{ dx_dX, state->nonlocalField };

    const std::tuple< double, double, double > eigenDeformation = { state->F0_XX, state->F0_YY, state->F0_ZZ };

    if constexpr ( nDim == 2 ) {
      if ( hasEigenDeformation )
        material->computePlaneStrain( response3D, algorithmicModuli3D, deformation3D, timeIncrement, eigenDeformation );
      else
        material->computePlaneStrain( response3D, algorithmicModuli3D, deformation3D, timeIncrement );
    }
    else {
      if ( hasEigenDeformation )
        material->computeStress( response3D, algorithmicModuli3D, deformation3D, timeIncrement, eigenDeformation );
      else
        material->computeStress( response3D, algorithmicModuli3D, deformation3D, timeIncrement );
    }

    // the in-plane part of a 3D tensor in plane strain, the tensor itself in 3D
    auto inPlane = []( const auto& tensor3D ) {
      if constexpr ( nDim == 2 ) {
        if constexpr ( std::decay_t< decltype( tensor3D ) >::rank() == 2 )
          return reduceTo2D< U, U >( tensor3D );
        else
          return reduceTo2D< U, U, U, U >( tensor3D );
      }
      else
        return tensor3D;
    };

    this->response.S              = inPlane( response3D.tau );
    this->response.dL             = response3D.L - state->localField;
    this->response.nonLocalRadius = response3D.nonLocalRadius;

    state->localField = response3D.L;
    state->stress     = response3D.tau;

    const auto& I = Spatial3D::I;

    // chain rule from d/dF to d/d(DeltaF), since the unknown of the increment is DeltaF = dx_dY
    const Tensor3333d dF_dDeltaF_3D = einsum< ij, JI, to_iIjJ >( I, state->dY_dX );
    const TensorDDDD  dF_dDeltaF    = inPlane( dF_dDeltaF_3D );

    tangents.dS_dDeltaF = inPlane( algorithmicModuli3D.dTau_dF );
    tangents.dS_dDeltaF = einsum< ijmn, mnKL >( tangents.dS_dDeltaF, dF_dDeltaF );

    for ( int K = 0; K < nNonlocal; K++ ) {
      tangents.dS_dN[K]      = inPlane( algorithmicModuli3D.dTau_dN[K] );
      tangents.dL_dDeltaF[K] = inPlane( algorithmicModuli3D.dL_dF[K] );
      tangents.dL_dDeltaF[K] = einsum< mn, mnKL >( tangents.dL_dDeltaF[K], dF_dDeltaF );
    }
    tangents.dL_dN = algorithmicModuli3D.dL_dN;

    // the out-of-plane parts in 2D, for the hoop terms of an axisymmetric particle (see OutOfPlaneTangents)
    if constexpr ( nDim == 2 ) {
      const auto& F_n = state->dY_dX;

      this->response.tau33 = response3D.tau( 2, 2 );
      tangents.outOfPlane  = computeOutOfPlaneTangents( algorithmicModuli3D.dTau_dF, F_n );
      for ( int K = 0; K < nNonlocal; K++ ) {
        tangents.dTau33_dN( K )    = algorithmicModuli3D.dTau_dN[K]( 2, 2 );
        tangents.dL_dDeltaF33( K ) = 0.0;
        for ( int I = 0; I < 3; I++ )
          tangents.dL_dDeltaF33( K ) += algorithmicModuli3D.dL_dF[K]( 2, I ) * F_n( 2, I );
      }
    }
  }

  template class GradientEnhancedFiniteStrainMaterialPoint< 2, 1 >;
  template class GradientEnhancedFiniteStrainMaterialPoint< 3, 1 >;
  template class GradientEnhancedFiniteStrainMaterialPoint< 3, 6 >;

} // namespace Marmot::Discretizations::MaterialPoints
