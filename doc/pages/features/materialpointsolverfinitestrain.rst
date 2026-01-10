Material Point Solver for Finite Strain Materials
=================================================

Overview
--------

The **Material Point Solver (Finite Strain)** provides a standalone driver for simulating
the response of a single material point using **finite strain constitutive models**
from the *MAteRialMOdellingToolbox (marmot)* framework.

Unlike the hypo-elastic solver, this class handles large deformations using full :math:`3\times3`
tensor kinematics. It solves for the **Displacement Gradient** (:math:`\nabla \mathbf{u}`)
and the **Kirchhoff Stress** (:math:`\boldsymbol{\tau}`).

It is designed for:

- Testing and verifying finite strain material models (e.g., hyperelasticity, finite plasticity),
- Performing mixed displacement-gradient/stress-controlled loading paths,
- Recording full tensor histories including material tangents and internal state variables.

The solver automatically handles time stepping, adaptive cutbacks upon non-convergence,
and output formatting.

Main Features
-------------

- **Finite Strain Kinematics**: Solves for :math:`\mathbf{F} = \mathbf{I} + \nabla \mathbf{u}`.
- **Mixed Control**: Supports independent control of displacement gradient or stress for each of the 9 tensor components.
- **Adaptive Time Stepping**: Automatically reduces time step size (:math:`\Delta t`) if convergence fails.
- **Newton–Raphson Solver**: robust iterative solver with configurable residual and correction tolerances.
- **History Recording**: Stores time, stress, deformation gradient, tangent moduli, and state variables.
- **CSV Export**: Exports simulation history for post-processing.

Solver Structure
----------------

The solver is organized around a few key data structures:

- **SolverOptions** — Defines numerical tolerances (residual/correction) and iteration limits.
- **Step** — Represents a loading phase. It defines target increments for :math:`\nabla \mathbf{u}` and ::math:`\boldsymbol{\tau}`, along with boolean flags determining which variable controls which tensor component.
- **Increment** — A sub-step automatically generated within each Step based on time discretization.
- **HistoryEntry** — Records the converged state (Time, :math:`\boldsymbol{\tau}`, :math:`\mathbf{F}`, :math:`\mathbb{C}`, state variables) at specific points.

Usage Workflow
--------------

A typical workflow consists of:

1. **Initialize** the solver with a material model name and properties.
2. **Set** the initial state (initial stress and internal variables).
3. **Define** one or more loading ``Step`` objects with control flags and targets.
4. **Add** the steps to the solver.
5. **Solve** the problem.
6. **Inspect or export** the resulting history.

Example Usage
-------------

Here is a simple example demonstrating how to set up and run the finite strain material point solver using a compressible Neo-Hookean material model.
It can be found in the `examples/finite-strain-mp-solver/CompressibleNeoHooke` directory of the Marmot repository together with the appropriate ``Makefile``.

.. literalinclude:: ../../../examples/finite-strain-mp-solver/CompressibleNeoHooke/example.cpp
   :language: c++
   :linenos:
   :caption: Finite Strain Solver Example

To compile and run the example, you can use the following ``Makefile``. Ensure that the paths to Marmot, Eigen, and Fastor match your system configuration.

.. literalinclude:: ../../../examples/finite-strain-mp-solver/CompressibleNeoHooke/Makefile
   :language: make
   :linenos:
   :caption: Makefile for Finite Strain Material Point Solver Example


Typical Output
--------------

Running the solver prints detailed convergence information to the console:

.. code-block:: none

   Solving step from 0 to 1
   +------------------------------------------------------------------------------+
     Solving increment 1, time: 0 to 0.1, dT: 0.1
       Iteration 0, ||ddGradU||: 0.00e+00, ||R||: 5.00e-02
       Iteration 1, ||ddGradU||: 1.25e-03, ||R||: 4.12e-05
       Iteration 2, ||ddGradU||: 2.10e-07, ||R||: 1.10e-11
       Converged after 3 iterations.
     Material state for time: 1.000000e-01
     tau:
      [2.500000e+01, 0.000000e+00, 0.000000e+00
       0.000000e+00, 0.000000e+00, 0.000000e+00
       0.000000e+00, 0.000000e+00, 0.000000e+00]
     ...

Guidelines and Notes
--------------------

- **Tensors**: The solver uses the **Fastor** library. Inputs usually expect ``Tensor33d`` (3x3 double) or ``Tensor9d`` (flattened 9x1 double).
- **Control Flags**: For every component :math:`(i,j)` of the :math:`3\times3` tensor, exactly **one** of ``isGradUComponentControlled(i,j)`` or ``isStressComponentControlled(i,j)`` must be set to true. The ``Step::checkControl()`` method enforces this.
- **Singularity Handling**: When using mixed control, the tangent matrix is modified to ensure invertibility. Rows corresponding to displacement-controlled components are zeroed and the diagonal set to 1.
- **Stress Measure**: The solver works with **Kirchhoff stress** (:math:`\boldsymbol{\tau}`). If the material model computes Cauchy stress, it must be converted internally or by the model.

CSV Export
----------

The history can be exported to CSV. The columns are formatted as follows:

- **Time**: Simulation time.
- **tau[9]**: Kirchhoff stress components (11, 12, 13, 21, 22, 23, 31, 32, 33).
- **F[9]**: Deformation gradient components (11, 12, 13, 21, 22, 23, 31, 32, 33).
- **SV[N]**: Internal state variables (SV1, SV2, ...).

.. doxygenclass:: MarmotMaterialPointSolverFiniteStrain
   :allow-dot-graphs:
