Material Point Solver for Hypo-Elastic Materials
================================================

Overview
--------

The **Material Point Solver (Hypo-Elastic)** provides a standalone driver for simulating
the response of a single material point using **hypo-elastic constitutive models**
from the *MAteRialMOdellingToolbox (Marmot)* framework.

It is designed for:

- Testing and verifying constitutive models,
- Performing mixed strain/stress-controlled loading paths,
- Recording full stress–strain histories including tangent matrices and internal variables.

The solver automatically handles time stepping, convergence checks, and output formatting.

Main Features
-------------

- Supports **mixed control** (independent strain or stress control for each Voigt component)
- **Adaptive time stepping** within each load step
- **Newton–Raphson** iteration with configurable convergence tolerances
- **History recording** for post-processing
- **CSV export** for analysis or plotting

Solver Structure
----------------

The solver is organized around a few key data structures:

- **SolverOptions** — defines numerical tolerances and iteration limits.
- **Step** — represents a loading phase with target increments, control flags, and time control.
- **Increment** — represents a sub-step automatically generated within each Step.
- **HistoryEntry** — records time, stress, strain, tangent, and state variables at each increment.

Usage Workflow
--------------

A typical workflow consists of:

1. **Initialize** the solver with a material model name and its properties.
2. **Set** the initial state (stress and internal variables).
3. **Define** one or more loading Steps with control flags and targets.
4. **Add** the steps to the solver.
5. **Solve** the problem.
6. **Inspect or export** the resulting history.

Example Usage
-------------

.. code-block:: c++

   #include "Marmot/MarmotMaterialPointSolverHypoElastic.h"
   #include <Eigen/Dense>
   #include <iostream>

   int main() {
     using Vec6 = Marmot::Vector6d;

     // 1) Define the material model (must be registered in Marmot)
     std::string materialName = "LINEARELASTIC";
     double properties[] = { 210e9, 0.3 }; // Example: Young's modulus, Poisson's ratio
     int nProps = 2;

     // 2) Configure solver options
     MarmotMaterialPointSolverHypoElastic::SolverOptions options;
     options.maxIterations       = 25;
     options.residualTolerance   = 1e-10;
     options.correctionTolerance = 1e-10;

     // 3) Create solver instance
     MarmotMaterialPointSolverHypoElastic solver(materialName, properties, nProps, options);

     // 4) Set the initial state
     int nSV = 0;
     solver.getNumberOfStateVariables(nSV);
     Eigen::VectorXd initialSV = Eigen::VectorXd::Zero(nSV);
     solver.setInitialState(Vec6::Zero(), initialSV);

     // 5) Define a loading step (uniaxial strain-controlled on ε_11)
     MarmotMaterialPointSolverHypoElastic::Step step;
     step.timeStart = 0.0;
     step.timeEnd   = 1.0;
     step.dTStart   = 0.2;
     step.dTMin     = 1e-6;
     step.dTMax     = 0.5;
     step.maxIncrements = 50;

     step.strainIncrementTarget = Vec6::Zero();
     step.stressIncrementTarget = Vec6::Zero();

     step.strainIncrementTarget[0] = 1e-3; // Total applied strain over the step

     // Control flags: exactly one of each must be true
     step.isStrainComponentControlled = Eigen::Vector<bool,6>::Constant(false);
     step.isStressComponentControlled = Eigen::Vector<bool,6>::Constant(false);
     step.isStrainComponentControlled[0] = true;  // ε_11 controlled
     for (int i = 1; i < 6; ++i)
       step.isStressComponentControlled[i] = true; // others stress-controlled (zero target)

     // 6) Run the solver
     solver.addStep(step);
     solver.solve();

     // 7) Output the results
     solver.printHistory();
     solver.exportHistoryToCSV("mp_history.csv");

     return 0;
   }

Typical Output
--------------

Running the above example will print information about each load step and increment, for example:

.. code-block:: none

   Solving step from 0 to 1
     Solving increment 1, time: 0 to 0.2, dT: 0.2
       Iteration 0, ||ddE||: 1.0e+00, ||R||: 2.5e-05
       Iteration 1, ||ddE||: 3.1e-07, ||R||: 1.2e-10
     Converged after 2 iterations.
   Material Point History:
   Time: 1.000000e+00
     Stress: 2.100000e+08  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00
     Strain: 1.000000e-03  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00
     State Vars: ...

The history can be exported as a CSV file (``mp_history.csv``) for visualization or analysis.

Guidelines and Notes
--------------------

- Each of the six Voigt components must be either **strain-controlled** or **stress-controlled**.
  ``Step::checkControl()`` is used to verify consistency when adding a step.
- The total target increments are distributed automatically over time according to ``dT``.
- If convergence fails, the solver halves ``dT`` until reaching ``dTMin``.
- The exported CSV includes:
  ``time``, ``stress[6]``, ``strain[6]``, and all state variables.

Practical Applications
----------------------

- Material verification and constitutive model calibration
- Generating reference stress–strain curves
- Investigating convergence and stability of material models
- Educational demonstrations of hypo-elastic material response


.. doxygenclass:: MarmotMaterialPointSolverHypoElastic
   :allow-dot-graphs:
