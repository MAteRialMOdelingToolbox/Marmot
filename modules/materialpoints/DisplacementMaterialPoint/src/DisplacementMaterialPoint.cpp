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
#include "Marmot/DisplacementMaterialPoint.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include <cstring>
#include <ctime>

namespace Marmot::MaterialPoints {

  void DisplacementMaterialPoint2D::incrementDeformation( const TensorD& du, const TensorDD& du_dY )
  {
    Marmot::Math::mapEigenToFastor( state->du ) += Marmot::Math::mapEigenToFastor( Marmot::Math::expandTo3D( du ) );

    Marmot::Math::mapEigenToFastor( state->dx_dY ) += Marmot::Math::mapEigenToFastor( Marmot::Math::expandTo3D( du_dY ) );
  };

  void DisplacementMaterialPoint2D::computeYourself( double timeNew, double dT )
  {

    const Material ::TimeIncrement timeIncrement{ timeNew, dT };

    // TODO: make auto&
    const Marmot::Math::FastorTensors::Tensor33d dx_dX = state->dx_dY % state->dY_dX;

    using namespace Marmot;

    Material::ConstitutiveResponse< 3 > response3D{ 0, 0, 0, state->materialState.data() };

    Material::AlgorithmicModuli< 3 > algorithmicModuli3D;

    Material::Deformation< 3 > deformation3D{ dx_dX };

    material->computePlaneStrain( response3D, algorithmicModuli3D, deformation3D, timeIncrement );

    // clang-format off
    this->response = { .S = Marmot::Math::reduceTo2D< Marmot::Math::U, Marmot::Math::U >( response3D.tau ),
    };

    state->stress = response3D.tau;

    using namespace Marmot::Math::FastorTensors;
    using namespace Marmot::Math::Indices;

    const auto& I = Marmot::Math::FastorTensors::Spatial3D::I;

    const Tensor3333d dF_dDeltaF_3D = einsum< ij, JI, to_iIjJ >(I, state->dY_dX);
    const TensorDDDD dF_dDeltaF = Marmot::Math::reduceTo2D<Marmot::Math::U,Marmot::Math::U,Marmot::Math::U,Marmot::Math::U>(dF_dDeltaF_3D);

    this->tangents = {
        Marmot::Math::reduceTo2D< Marmot::Math::U, Marmot::Math::U, Marmot::Math::U, Marmot::Math::U >    ( algorithmicModuli3D.dTau_dF),
    };

    tangents.dS_dDeltaF = einsum<ijmn, mnKL>(tangents.dS_dDeltaF, dF_dDeltaF);
    // clang-format on

    _density = material->getDensity( state->materialState.data() );
  }

  void DisplacementMaterialPoint3D::incrementDeformation( const TensorD& du, const TensorDD& du_dY )
  {
    Marmot::Math::mapEigenToFastor( state->du ) += Marmot::Math::mapEigenToFastor( ( du ) );
    Marmot::Math::mapEigenToFastor( state->dx_dY ) += Marmot::Math::mapEigenToFastor( ( du_dY ) );
  };

  void DisplacementMaterialPoint3D::computeYourself( double timeNew, double dT )
  {

    const Material ::TimeIncrement timeIncrement{ timeNew, dT };

    /* const Marmot::Math::FastorTensors::Tensor33d dx_dX_n  = state->dY_dX; */
    const Marmot::Math::FastorTensors::Tensor33d dx_dX = state->dx_dY % state->dY_dX;

    using namespace Marmot;

    Material::ConstitutiveResponse< 3 > response3D{ 0, 0, 0, state->materialState.data() };

    Material::AlgorithmicModuli< 3 > algorithmicModuli3D;

    Material::Deformation< 3 > deformationIncrement3D{ dx_dX };

    material->computeStress( response3D, algorithmicModuli3D, deformationIncrement3D, timeIncrement );

    // clang-format off
    this->response = { .S = ( response3D.tau ),
    };

    state->stress = response3D.tau;

    using namespace Marmot::Math::FastorTensors;
    using namespace Marmot::Math::Indices;

    const auto& I = Marmot::Math::FastorTensors::Spatial3D::I;

    const Tensor3333d dF_dDeltaF_3D = einsum< ij, JI, to_iIjJ >(I, state->dY_dX);
    const TensorDDDD dF_dDeltaF = (dF_dDeltaF_3D);

    this->tangents = { ( algorithmicModuli3D.dTau_dF ) };

    tangents.dS_dDeltaF = einsum<ijmn, mnKL>(tangents.dS_dDeltaF, dF_dDeltaF);
    // clang-format on

    _density = material->getDensity( state->materialState.data() );
  }

} // namespace Marmot::MaterialPoints
