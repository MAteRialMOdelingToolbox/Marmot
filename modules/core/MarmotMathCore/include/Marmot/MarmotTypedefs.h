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
#include "Eigen/Dense"
#include "autodiff/forward/dual/dual.hpp"
#include "autodiff/forward/dual/eigen.hpp"
#include "unsupported/Eigen/CXX11/Tensor"

/**
 * @file MarmotTypedefs.h
 * @brief Common Eigen matrix/vector type aliases used throughout Marmot.
 *
 * This header collects all fixed-size and map-based Eigen types that are shared
 * across the Marmot library, so that downstream code can rely on a single,
 * consistent set of names (e.g. @c Marmot::Vector6d, @c Marmot::Matrix6d).
 */

/** @brief Macro to compute the Voigt notation index from the dimension */
#define VOIGTFROMDIM( x ) ( ( ( x * x ) + x ) >> 1 )

namespace Marmot {
  typedef Eigen::Matrix< double, 6, 6 > Matrix6d;              ///< 6×6 double matrix (Voigt tangent)
  typedef Eigen::Matrix< double, 6, 9 > Matrix69d;             ///< 6×9 double matrix
  typedef Eigen::Matrix< double, 9, 9 > Matrix99d;             ///< 9×9 double matrix
  typedef Eigen::Matrix< double, 3, 4 > Matrix34d;             ///< 3×4 double matrix
  typedef Eigen::Map< Matrix6d >        mMatrix6d;             ///< Non-owning map over a Matrix6d
  typedef Eigen::Matrix< double, 3, 3 > Matrix3d;              ///< 3×3 double matrix

  typedef Eigen::Matrix< double, 3, 1 >        Vector3d;       ///< 3-component double column vector
  typedef Eigen::Matrix< double, 4, 1 >        Vector4d;       ///< 4-component double column vector
  typedef Eigen::Matrix< double, 6, 1 >        Vector6d;       ///< 6-component double column vector (Voigt notation)
  typedef Eigen::Matrix< double, 7, 1 >        Vector7d;       ///< 7-component double column vector
  typedef Eigen::Matrix< double, 8, 1 >        Vector8d;       ///< 8-component double column vector
  typedef Eigen::Matrix< double, 9, 1 >        Vector9d;       ///< 9-component double column vector
  typedef Eigen::Matrix< int, 8, 1 >           Vector8i;       ///< 8-component integer column vector
  typedef Eigen::Matrix< double, 1, 6 >        RowVector6d;    ///< 6-component double row vector
  typedef Eigen::Map< Vector6d >               mVector6d;      ///< Non-owning map over a Vector6d
  typedef Eigen::Map< Eigen::VectorXd >        mVectorXd;      ///< Non-owning map over a dynamic double vector
  typedef Eigen::Map< const Marmot::Vector6d > mConstVector6d; ///< Non-owning const map over a Vector6d

  typedef Eigen::Matrix< double, 3, 6 > Matrix36d;             ///< 3×6 double matrix
  typedef Eigen::Matrix< double, 3, 6 > Matrix36;              ///< 3×6 double matrix (alias for Matrix36d)
  typedef Eigen::Matrix< double, 6, 3 > Matrix63d;             ///< 6×3 double matrix
  typedef Eigen::Matrix< double, 9, 9 > Matrix9d;              ///< 9×9 double matrix (alias for Matrix99d)

  // complex matrix definitions
  typedef std::complex< double >               complexDouble; ///< Standard complex double scalar
  typedef Eigen::Matrix< complexDouble, 6, 1 > Vector6cd;     ///< 6-component complex double column vector

  // autodiff dual matrix definitions
  typedef Eigen::Matrix< autodiff::dual, 6, 1 > Vector6dual;       ///< 6-component autodiff::dual column vector
  typedef Eigen::Map< Vector6dual >             mVector6dual;      ///< Non-owning map over a Vector6dual
  typedef Eigen::Map< const Vector6dual >       mVector6dualConst; ///< Non-owning const map over a Vector6dual

  // definitions for templated scalar type
  template < typename T >
  using Vector6t = Eigen::Matrix< T, 6, 1 >; ///< 6-component column vector with scalar type @p T

  template < typename T >
  using VectorXt = Eigen::Matrix< T, -1, 1 >; ///< Dynamic column vector with scalar type @p T

  namespace TensorUtility::EigenTensors {

    typedef Eigen::TensorFixedSize< double, Eigen::Sizes< 6, 3, 3 > > Tensor633d; ///< Fixed-size 6×3×3 double tensor
    typedef Eigen::TensorFixedSize< double, Eigen::Sizes< 3, 2, 2 > > Tensor322d; ///< Fixed-size 3×2×2 double tensor
    typedef Eigen::TensorFixedSize< double, Eigen::Sizes< 3, 3, 3, 3 > >
      Tensor3333d; ///< Fixed-size 3×3×3×3 double tensor
    typedef Eigen::TensorFixedSize< double, Eigen::Sizes< 3, 3, 3 > > Tensor333d; ///< Fixed-size 3×3×3 double tensor
    typedef Eigen::TensorFixedSize< double, Eigen::Sizes< 1, 2, 2 > > Tensor122d; ///< Fixed-size 1×2×2 double tensor
    typedef Eigen::TensorFixedSize< double, Eigen::Sizes< 2, 2, 2, 2 > >
      Tensor2222d;                                                              ///< Fixed-size 2×2×2×2 double tensor
    typedef Eigen::TensorFixedSize< double, Eigen::Sizes< 2, 2, 1, 2 > >
      Tensor2212d;                                                              ///< Fixed-size 2×2×1×2 double tensor
    typedef Eigen::TensorFixedSize< double, Eigen::Sizes< 2, 1, 2, 2 > >
      Tensor2122d;                                                              ///< Fixed-size 2×1×2×2 double tensor
    typedef Eigen::TensorFixedSize< double, Eigen::Sizes< 2, 1, 1, 2 > >
                                                                   Tensor2112d; ///< Fixed-size 2×1×1×2 double tensor
    typedef Eigen::TensorFixedSize< double, Eigen::Sizes< 3, 3 > > Tensor33d;   ///< Fixed-size 3×3 double tensor

  }                                                                             // namespace TensorUtility::EigenTensors

} // namespace Marmot
