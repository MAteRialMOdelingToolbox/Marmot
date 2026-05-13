/* ---------------------------------------------------------------------
 *                                       _
 *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_
 * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
 * | | | | | | (_| | |  | | | | | | (_) | |_
 * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
 *
 * Unit of Strength of Materials and Structural Analysis
 * University of Innsbruck,
 * 2020 - today
 *
 * festigkeitslehre@uibk.ac.at
 *
 * Alexandros Stathas alexandros.stathas@boku.ac.at
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
#pragma once
#include "Marmot/MarmotFastorTensorBasics.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/src/Core/Matrix.h>
#include <Fastor/Fastor.h>
#include <Fastor/tensor_algebra/indicial.h>
#include <cmath>
#include <tuple>
#include <unsupported/Eigen/CXX11/Tensor>

using namespace Eigen;
using namespace Fastor;

using Marmot::FastorStandardTensors::Tensor3333d;
using Marmot::FastorStandardTensors::Tensor333d;
using Marmot::FastorStandardTensors::Tensor33d;
using Marmot::FastorStandardTensors::Tensor3d;

namespace Marmot::Materials {

  namespace InterfaceMaterialHelperFunctions {

    Tensor33d compute_inv( const Tensor33d& Q );

    std::tuple< Tensor3333d, const Tensor3333d, Tensor3333d, Tensor33d > interfaceGeometrySystemCouplings(
      const Tensor33d&   N,
      const Tensor33d&   T,
      const Tensor3333d& L );

    std::tuple< Tensor3333d, Tensor33d, Tensor333d, Tensor3333d > calculateMaterialMatrices(
      const Tensor3d&    normal,
      const Tensor33d&   N,
      const Tensor33d&   T,
      const Tensor3333d& C_nu_aibj );

    Eigen::Matrix< double, 9, 9 > convert4thOrderTensorToMatrix_9x9( const Tensor3333d& tensor );
    Eigen::Matrix< double, 9, 3 > convert3rdOrderTensorToMatrix_9x3( const Tensor333d& tensor );
    Eigen::Matrix< double, 3, 9 > convert3rdOrderTensorToMatrix_3x9( const Tensor333d& tensor );
    Eigen::Matrix< double, 3, 3 > convert2ndOrderTensorToMatrix_3x3( const Tensor33d& tensor );

    std::tuple< Tensor3333d, Tensor3333d, Tensor3333d, Tensor3333d, Tensor33d, Tensor3333d > calculateFY(
      const Tensor33d&   N,
      const Tensor33d&   T,
      const Tensor3333d& C_nu_aibj );

    std::tuple< Tensor3333d, Tensor33d, Tensor333d, Tensor3333d > calculateInterfaceMaterialParameters(
      const Tensor3d& normal,
      const double&   nu_0 );

    std::tuple< Tensor3333d, Tensor33d, Tensor333d, Tensor3333d > calculateInterfaceMaterialParameters(
      const Tensor3d&                      normal,
      const Eigen::Matrix< double, 6, 6 >& C_ep_voigt );
  } // namespace InterfaceMaterialHelperFunctions
} // namespace Marmot::Materials
