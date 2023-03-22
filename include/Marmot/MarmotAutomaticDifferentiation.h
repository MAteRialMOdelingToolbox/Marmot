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
 * Alexander Dummer alexander.dummer@uibk.ac.at
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
#include "autodiff/forward/dual.hpp"
#include "autodiff/forward/dual/eigen.hpp"
#include <functional>

namespace Marmot {
  namespace AutomaticDifferentiation {
    using function_type = std::function< autodiff::VectorXdual( const autodiff::VectorXdual& X ) >;
    Eigen::MatrixXd forwardMode( const function_type& F, const Eigen::VectorXd& X );

    std::tuple< Eigen::VectorXd, Eigen::MatrixXd > dF_dX( const function_type& F, const Eigen::VectorXd& X );
  } // namespace AutomaticDifferentiation

} // namespace Marmot
