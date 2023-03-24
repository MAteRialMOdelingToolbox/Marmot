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
#include "Marmot/MarmotTensor.h"
#include "autodiff/forward/dual.hpp"
#include "autodiff/forward/dual/eigen.hpp"
#include <autodiff/forward/dual/dual.hpp>
#include <functional>

namespace Marmot {

  namespace AutomaticDifferentiation {

    using namespace autodiff;
    using namespace Eigen;

    dual2nd shiftTo2ndOrderDual( const dual& x );

    using scalar_to_scalar_function_type = std::function< dual( const dual& ) >;
    double df_dx( const scalar_to_scalar_function_type& f, const double& x );

    using scalar_to_scalar_function_type_2nd = std::function< dual2nd( const dual2nd& ) >;
    dual df_dx( const scalar_to_scalar_function_type_2nd& f, const dual& x );

    using vector_to_vector_function_type = std::function< VectorXdual( const VectorXdual& X ) >;
    MatrixXd forwardMode( const vector_to_vector_function_type& F, const VectorXd& X );

    using vector_to_vector_function_type_dual = std::function< VectorXdual( const VectorXdual& X ) >;
    std::pair< VectorXd, MatrixXd > jacobian( const vector_to_vector_function_type_dual& F, const VectorXd& X );

  } // namespace AutomaticDifferentiation

} // namespace Marmot
