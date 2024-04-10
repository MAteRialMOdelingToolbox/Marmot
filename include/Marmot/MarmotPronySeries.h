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

namespace Marmot::Materials {

  namespace PronySeries {
    using namespace Marmot;
    using namespace Eigen;

    struct Properties {
      size_t                  nPronyTerms;
      Matrix6d                ultimateStiffnessMatrix;
      Matrix< double, 6, -1 > pronyStiffnesses;
      Matrix< double, 6, -1 > pronyRelaxationTimes;
    };

    typedef Eigen::Matrix< double, 6, Eigen::Dynamic > StateVarMatrix;
    typedef Eigen::Map< StateVarMatrix >               mapStateVarMatrix;

    void evaluatePronySeries( const Properties&               props,
                              Vector6d&                       stress,
                              Matrix6d&                       stiffness,
                              Eigen::Ref< mapStateVarMatrix > stateVars,
                              const Vector6d&                 dStrain,
                              const double                    dT,
                              const bool                      updateStateVars = false );

    void updateStateVars( const Properties&               props,
                          Eigen::Ref< mapStateVarMatrix > stateVars,
                          const Vector6d&                 dStrain,
                          const double                    dT );

  } // namespace PronySeries
} // namespace Marmot::Materials
