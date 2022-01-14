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
 * Matthias Neuner matthias.neuner@uibk.ac.at
 * Magdalena Schreter magdalena.schreter@uibk.ac.at
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
#include "Eigen/Core"
#include "Eigen/Dense"
#include "unsupported/Eigen/CXX11/Tensor"

namespace Marmot {
  typedef Eigen::Matrix< double, 6, 6 > Matrix6d;
  typedef Eigen::Matrix< double, 6, 9 > Matrix69d;
  typedef Eigen::Matrix< double, 9, 9 > Matrix99d;
  typedef Eigen::Matrix< double, 3, 4 > Matrix34d;
  typedef Eigen::Map< Matrix6d >        mMatrix6d;
  typedef Eigen::Matrix< double, 3, 3 > Matrix3d;

  typedef Eigen::Matrix< double, 3, 1 >        Vector3d;
  typedef Eigen::Matrix< double, 6, 1 >        Vector6d;
  typedef Eigen::Matrix< double, 7, 1 >        Vector7d;
  typedef Eigen::Matrix< double, 8, 1 >        Vector8d;
  typedef Eigen::Matrix< double, 9, 1 >        Vector9d;
  typedef Eigen::Matrix< int, 8, 1 >           Vector8i;
  typedef Eigen::Matrix< double, 1, 6 >        RowVector6d;
  typedef Eigen::Map< Vector6d >               mVector6d;
  typedef Eigen::Map< Eigen::VectorXd >        mVectorXd;
  typedef Eigen::Map< const Marmot::Vector6d > mConstVector6d;

  typedef Eigen::Matrix< double, 3, 6 > Matrix36d;
  typedef Eigen::Matrix< double, 3, 6 > Matrix36;
  typedef Eigen::Matrix< double, 6, 3 > Matrix63d;
  typedef Eigen::Matrix< double, 9, 9 > Matrix9d;

  // complex matrix definitions
  typedef std::complex< double >               complexDouble;
  typedef Eigen::Matrix< complexDouble, 6, 1 > Vector6cd;

  // definitions for dual numbers

  namespace EigenTensors {

    typedef Eigen::TensorFixedSize< double, Eigen::Sizes< 6, 3, 3 > >    Tensor633d;
    typedef Eigen::TensorFixedSize< double, Eigen::Sizes< 3, 2, 2 > >    Tensor322d;
    typedef Eigen::TensorFixedSize< double, Eigen::Sizes< 3, 3, 3, 3 > > Tensor3333d;
    typedef Eigen::TensorFixedSize< double, Eigen::Sizes< 3, 3, 3 > >    Tensor333d;
    typedef Eigen::TensorFixedSize< double, Eigen::Sizes< 1, 2, 2 > >    Tensor122d;
    typedef Eigen::TensorFixedSize< double, Eigen::Sizes< 2, 2, 2, 2 > > Tensor2222d;
    typedef Eigen::TensorFixedSize< double, Eigen::Sizes< 2, 2, 1, 2 > > Tensor2212d;
    typedef Eigen::TensorFixedSize< double, Eigen::Sizes< 2, 1, 2, 2 > > Tensor2122d;
    typedef Eigen::TensorFixedSize< double, Eigen::Sizes< 2, 1, 1, 2 > > Tensor2112d;

  } // namespace EigenTensors

} // namespace Marmot
