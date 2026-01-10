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

#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotTypedefs.h"
#include <iostream>
#pragma once

namespace Marmot {
  namespace Solvers {
    /**
     * @brief Solver for material point problems with hypo-elastic materials
     * @details This class implements a solver for material point problems
     * using hypo-elastic material models. It supports loading steps with
     * controlled strain and stress components, adaptive time stepping,
     * and history recording.
     */
    class MarmotMaterialPointSolverHypoElastic {

    public:
      /**
       * @struct Step
       * @brief Struct to define a loading step
       * @details Each step contains target strain and stress states,
       * time information and time step control parameters.
       */
      struct Step {
        Marmot::Vector6d strainIncrementTarget; ///< Target strain increment for the step
        Marmot::Vector6d stressIncrementTarget; ///< Target stress increment for the step
        Eigen::Vector< bool, 6 >
          isStrainComponentControlled;          ///< Flags to indicate which strain components are controlled
        Eigen::Vector< bool, 6 >
               isStressComponentControlled;     ///< Flags to indicate which stress components are controlled
        double timeStart     = 0.0;             ///< Start time of the step
        double timeEnd       = 1.0;             ///< End time of the step
        double dTStart       = 0.1;             ///< Initial time step size
        double dTMin         = 1e-6;            ///< Minimum time step size
        double dTMax         = 0.5;             ///< Maximum time step size
        int    maxIncrements = 100;             ///< Maximum number of increments in the step

        /**
         * @brief Check that for each component, either strain or stress is controlled
         * @throws std::runtime_error if the condition is not met
         */
        void checkControl() const
        {
          for ( int i = 0; i < 6; i++ ) {
            if ( isStrainComponentControlled[i] == isStressComponentControlled[i] ) {
              throw std::runtime_error(
                "exactly one of strain or stress component must be controlled for each component." );
            }
          }
        }
      };
      /**
       * @struct Increment
       * @brief Struct to define a loading increment
       * @details Each increment contains strain and stress increments,
       * control flags, time information, and iteration limits.
       */
      struct Increment {
        Marmot::Vector6d strainIncrement;   ///< Strain increment for the increment
        Marmot::Vector6d stressIncrement;   ///< Stress increment for the increment
        Eigen::Vector< bool, 6 >
          isStrainComponentControlled;      ///< Flags to indicate which strain components are controlled
        Eigen::Vector< bool, 6 >
               isStressComponentControlled; ///< Flags to indicate which stress components are controlled
        double timeOld;                     ///< Old time at the beginning of the increment
        double dT;                          ///< Time step size for the increment
      };
      /**
       * @struct HistoryEntry
       * @brief Struct to record the history of the simulation
       * @details Each entry contains time, stress, strain, and state variables.
       */
      struct HistoryEntry {
        double           time;           ///< Time at the history entry
        Marmot::Vector6d stress;         ///< Stress at the history entry
        Marmot::Vector6d strain;         ///< Strain at the history entry
        Marmot::Matrix6d dStressdStrain; ///< Material tangent at the history entry
        Eigen::VectorXd  stateVars;      ///< State variables at the history entry

        void print() const
        {
          std::cout.precision( 6 );
          std::cout << std::scientific << "Time: " << time << std::endl;
          std::cout << "  Stress:     " << stress.transpose() << std::endl;
          std::cout << "  Strain:     " << strain.transpose() << std::endl;
          std::cout << "  dStress_dStrain:\n" << dStressdStrain.transpose() << std::endl;
          std::cout << "  State Vars: " << stateVars.transpose() << std::endl;
        }
      };

      /**
       * @struct SolverOptions
       * @brief Struct to define solver options
       * @details Contains parameters for controlling the solver's behavior.
       */
      struct SolverOptions {
        int    maxIterations       = 25;    ///< Maximum number of iterations per increment
        double residualTolerance   = 1e-10; ///< Convergence tolerance
        double correctionTolerance = 1e-10; ///< Correction tolerance
      };

      /**
       * @brief Constructor for the MarmotMaterialPointSolverHypoElastic class
       * @param materialName Name of the hypo-elastic material model
       * @param materialProperties Array of material properties
       * @param nMaterialProperties Number of material properties
       */
      MarmotMaterialPointSolverHypoElastic( std::string&         materialName,
                                            double*              materialProperties,
                                            int                  nMaterialProperties,
                                            const SolverOptions& options );

      /**
       * @brief Add a loading step to the solver
       * @param step The Step to be added
       */
      void addStep( const Step& step );

      /**
       * @brief Get the list of added loading steps
       * @return A vector of Step containing the added steps
       */
      std::vector< Step > getSteps() const { return steps; }

      /**
       * @brief Clear all added loading steps
       */
      void clearSteps() { steps.clear(); }

      /**
       * @brief Set the initial state of the material model
       * @param initialStress The initial stress in Voigt notation
       * @param initialStateVars The initial state variables
       */
      void setInitialState( const Marmot::Vector6d& initialStress, const Eigen::VectorXd& initialStateVars );

      /**
       * @brief Get the number of state variables in the material model
       * @return The number of state variables
       */
      int getNumberOfStateVariables() const { return nStateVars; }

      /**
       * @brief Reset the solver to the initial state
       * @details This function resets the initial stress
       * and state variables of the material model.
       */
      void resetToInitialState();

      /**
       * @brief Solve the material point problem for all added steps
       */
      void solve();

      /**
       * @brief Get the recorded history of the simulation
       * @return A vector of HistoryEntry containing the recorded history
       */
      std::vector< HistoryEntry > getHistory() const { return history; }

      /**
       * @brief Clear the recorded history
       */
      void clearHistory() { history.clear(); }

      /**
       * @brief Print the recorded history to the console
       */
      void printHistory();

      /**
       * @brief Export the recorded history to a CSV file
       * @param filename The name of the CSV file to export to
       */
      void exportHistoryToCSV( const std::string& filename );

    private:
      /**
       * @brief Solve a single loading step
       * @param step The Step to be solved
       *
       * @details This function iterates over the increments
       * defined in the step, calling solveIncrement for each increment.
       * It manages time stepping and ensures that the entire step
       * is covered.
       */
      void solveStep( const Step& step );

      /**
       * @brief Solve a single increment within a loading step
       * @param increment The Increment to be solved
       * @throws std::runtime_error if the solver does not converge
       *
       * @details This function implements a Newton-Raphson iterative
       * solver to compute the stress and strain state for the given increment.
       * It updates the material state variables and records the history
       * after convergence.
       *
       */
      void solveIncrement( const Increment& increment );

      /**
       * @brief Compute the residual for the current increment
       * @param stressIncrement The computed stress increment
       * @param target The target (mixed) stress/strain increment
       * @param increment The Increment containing control information
       * @return The computed residual vector
       *
       * @details This function calculates the residual vector
       * based on the difference between the computed stress increment
       * and the target increment, taking into account which components
       * are controlled by strain or stress.
       */
      Marmot::Vector6d computeResidual( const Marmot::Vector6d& stressIncrement,
                                        const Marmot::Vector6d& target,
                                        const Increment&        increment );

      /**
       * @brief Modify the material tangent matrix based on control type
       * @param tangent The material tangent matrix to be modified
       * @param increment The Increment containing control information
       *
       * @details This function adjusts the tangent matrix to account for
       * the components that are controlled by strain or stress, ensuring
       * that the solver correctly handles mixed control scenarios.
       * This is done by zeroing out rows corresponding to strain-controlled
       * components and setting their diagonal entries to one.
       */
      void modifyTangent( Eigen::Matrix< double, 6, 6 >& tangent, const Increment& increment );

      /// @brief The hypo-elastic material model
      MarmotMaterialHypoElastic* material;

      /// @brief Number of state variables in the material model
      int nStateVars;

      /// @brief Current state variables
      Eigen::VectorXd stateVars;

      /// @brief Initial state variables
      Eigen::VectorXd _initialStateVars;

      /// @brief Temporary state variables for computations
      Eigen::VectorXd stateVarsTemp;

      /// @brief The stress in Voigt notation
      Marmot::Vector6d stress = Marmot::Vector6d::Zero();

      /// @brief The initial stress in Voigt notation
      Marmot::Vector6d _initialStress = Marmot::Vector6d::Zero();

      /// @brief The strain in Voigt notation
      Marmot::Vector6d strain = Marmot::Vector6d::Zero();

      /// @brief The material tangent matrix in Voigt notation
      Marmot::Matrix6d dStressDStrain = Marmot::Matrix6d::Zero();

      /// @brief List of loading steps
      std::vector< Step > steps;

      /// @brief History of the simulation
      std::vector< HistoryEntry > history;

      /// @brief Solver options
      const SolverOptions options;
    };
  } // namespace Solvers
} // namespace Marmot
