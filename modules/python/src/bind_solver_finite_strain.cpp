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
#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include "Marmot/MarmotMaterialPointSolverFiniteStrain.h"
#include "fastor_nanobind_caster.h"

namespace nb = nanobind;
using Solver = Marmot::Solvers::MarmotMaterialPointSolverFiniteStrain;

void bind_finite_strain_solver( nb::module_& m )
{
  nb::class_< Solver > solver( m, "FiniteStrainSolver" );

  nb::class_< Solver::Step >( solver, "Step" )
    .def( nb::init<>() )
    .def_rw( "gradUIncrementTarget", &Solver::Step::gradUIncrementTarget )
    .def_rw( "stressIncrementTarget", &Solver::Step::stressIncrementTarget )
    .def_rw( "isGradUComponentControlled", &Solver::Step::isGradUComponentControlled )
    .def_rw( "isStressComponentControlled", &Solver::Step::isStressComponentControlled )
    .def_rw( "timeStart", &Solver::Step::timeStart )
    .def_rw( "timeEnd", &Solver::Step::timeEnd )
    .def_rw( "dTStart", &Solver::Step::dTStart )
    .def_rw( "dTMin", &Solver::Step::dTMin )
    .def_rw( "dTMax", &Solver::Step::dTMax )
    .def_rw( "maxIncrements", &Solver::Step::maxIncrements )
    .def( "checkControl", &Solver::Step::checkControl );

  nb::class_< Solver::Increment >( solver, "Increment" )
    .def( nb::init<>() )
    .def_rw( "gradUIncrement", &Solver::Increment::gradUIncrement )
    .def_rw( "stressIncrement", &Solver::Increment::stressIncrement )
    .def_rw( "isGradUComponentControlled", &Solver::Increment::isGradUComponentControlled )
    .def_rw( "isStressComponentControlled", &Solver::Increment::isStressComponentControlled )
    .def_rw( "timeOld", &Solver::Increment::timeOld )
    .def_rw( "dT", &Solver::Increment::dT );

  nb::class_< Solver::SolverOptions >( solver, "SolverOptions" )
    .def( nb::init<>() )
    .def_rw( "maxIterations", &Solver::SolverOptions::maxIterations )
    .def_rw( "residualTolerance", &Solver::SolverOptions::residualTolerance )
    .def_rw( "correctionTolerance", &Solver::SolverOptions::correctionTolerance )
    .def( "__repr__", []( const Solver::SolverOptions& o ) {
      return "<SolverOptions maxIterations=" + std::to_string( o.maxIterations ) +
             " residualTolerance=" + std::to_string( o.residualTolerance ) +
             " correctionTolerance=" + std::to_string( o.correctionTolerance ) + ">";
    } );

  nb::class_< Solver::HistoryEntry >( solver, "HistoryEntry" )
    .def( nb::init<>() )
    .def_ro( "time", &Solver::HistoryEntry::time )
    .def_ro( "stress", &Solver::HistoryEntry::stress )
    .def_ro( "F", &Solver::HistoryEntry::F )
    .def_ro( "dTau_dF", &Solver::HistoryEntry::dTau_dF )
    .def_ro( "stateVars", &Solver::HistoryEntry::stateVars )
    .def( "print", &Solver::HistoryEntry::print )
    .def( "__repr__", []( const Solver::HistoryEntry& h ) {
      return "<HistoryEntry time=" + std::to_string( h.time ) + " stress[0,0]=" + std::to_string( h.stress( 0, 0 ) ) +
             ">";
    } );

  solver
    .def(
      "__init__",
      []( Solver*                                            t,
          const std::string&                                 name,
          nb::ndarray< double, nb::ndim< 1 >, nb::c_contig > props,
          const Solver::SolverOptions&                       opts ) {
        new ( t ) Solver( name, props.data(), static_cast< int >( props.size() ), opts );
      },
      nb::arg( "name" ),
      nb::arg( "props" ),
      nb::arg( "opts" ) = Solver::SolverOptions{} )
    .def( "addStep", &Solver::addStep, nb::arg( "step" ) )
    .def( "getSteps", &Solver::getSteps )
    .def( "clearSteps", &Solver::clearSteps )
    .def( "setInitialState", &Solver::setInitialState, nb::arg( "initialStress" ), nb::arg( "initialStateVars" ) )
    .def( "getNumberOfStateVariables", &Solver::getNumberOfStateVariables )
    .def( "resetToInitialState", &Solver::resetToInitialState )
    .def( "solve", &Solver::solve )
    .def( "getHistory", &Solver::getHistory )
    .def( "clearHistory", &Solver::clearHistory )
    .def( "printHistory", &Solver::printHistory )
    .def( "exportHistoryToCSV", &Solver::exportHistoryToCSV, nb::arg( "filename" ) );
}
