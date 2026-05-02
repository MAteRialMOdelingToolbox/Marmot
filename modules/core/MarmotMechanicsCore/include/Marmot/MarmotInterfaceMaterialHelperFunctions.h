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

using Tensor1D = Fastor::Tensor< double, 3 >;
using Tensor2D = Fastor::Tensor< double, 3, 3 >;
using Tensor3D = Fastor::Tensor< double, 3, 3, 3 >;
using Tensor4D = Fastor::Tensor< double, 3, 3, 3, 3 >;

namespace Marmot::Materials {

  namespace InterfaceMaterialHelperFunctions {

    Tensor2D compute_inv( const Tensor2D& I, Tensor2D& Q );

    std::tuple< Tensor4D, const Tensor4D, Tensor4D, Tensor2D > interfaceGeometrySystemCouplings( const Tensor2D& I,
                                                                                                 const Tensor2D& N,
                                                                                                 const Tensor2D& T,
                                                                                                 const Tensor4D& L );

    std::tuple< Tensor4D, Tensor2D, Tensor3D, Tensor4D > calculateMaterialMatrices( const Tensor1D& normal,
                                                                                    const Tensor2D& I,
                                                                                    const Tensor2D& N,
                                                                                    const Tensor2D& T,
                                                                                    const Tensor4D& C_nu_aibj );

    Tensor4D voigtToStiffness( const Eigen::Matrix< double, 6, 6 >& voigtStiffness );

    Eigen::Matrix< double, 9, 9 > convert4thOrderTensorToMatrix_9x9( const Tensor4D& tensor );
    Eigen::Matrix< double, 9, 3 > convert3rdOrderTensorToMatrix_9x3( const Tensor3D& tensor );
    Eigen::Matrix< double, 3, 9 > convert3rdOrderTensorToMatrix_3x9( const Tensor3D& tensor );
    Eigen::Matrix< double, 3, 3 > convert2ndOrderTensorToMatrix_3x3( const Tensor2D& tensor );

    std::tuple< Tensor4D, Tensor4D, Tensor4D, Tensor4D, Tensor2D, Tensor4D > calculateFY( const Tensor2D& I,
                                                                                          const Tensor2D& N,
                                                                                          const Tensor2D& T,
                                                                                          const Tensor4D& C_nu_aibj );

    std::tuple< Tensor4D, Tensor2D, Tensor3D, Tensor4D > calculateInterfaceMaterialParameters( const Tensor1D& normal,
                                                                                               const double&   nu_0 );

    std::tuple< Tensor4D, Tensor2D, Tensor3D, Tensor4D > calculateInterfaceMaterialParameters(
      const Tensor1D&                      normal,
      const Eigen::Matrix< double, 6, 6 >& C_ep_voigt );
  } // namespace InterfaceMaterialHelperFunctions
} // namespace Marmot::Materials
