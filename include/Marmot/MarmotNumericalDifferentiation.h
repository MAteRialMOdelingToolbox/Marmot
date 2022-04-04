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
#include "Marmot/MarmotTypedefs.h"
#include <functional>

namespace Marmot {
  namespace NumericalAlgorithms::Differentiation {

    using function_type = std::function< Eigen::VectorXd( const Eigen::VectorXd& X ) >;
    Eigen::MatrixXd forwardDifference( const function_type& F, const Eigen::VectorXd& X );
    Eigen::MatrixXd centralDifference( const function_type& F, const Eigen::VectorXd& X );

    namespace Complex {

      using function_type_complex = std::function< Eigen::VectorXcd( const Eigen::VectorXcd& X ) >;

      Eigen::MatrixXd forwardDifference( const function_type_complex& F, const Eigen::VectorXd& X );

      Eigen::MatrixXd centralDifference( const function_type_complex& F, const Eigen::VectorXd& X );

      Eigen::MatrixXd fourthOrderAccurateDerivative( const function_type_complex& F, const Eigen::VectorXd& X );

    } // namespace Complex
  }   // namespace NumericalAlgorithms::Differentiation
} // namespace Marmot
