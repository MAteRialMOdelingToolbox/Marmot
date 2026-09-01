import marmot
import marmot.testing
import numpy as np

print("Running example for FiniteStrain material: BERGSTROMBOYCE")

# Material properties: muA, kappaA, muB, kappaB, c1, c2, c3, implementationType (0 = CSDA)
properties = np.array([100.0, 1000.0, 50.0, 1000.0, 0.05, 1.0, 1.0, 0.0], dtype=np.float64)

# Setup solver
options = marmot.solvers.FiniteStrainSolver.SolverOptions()
solver = marmot.solvers.FiniteStrainSolver("BERGSTROMBOYCE", properties, options)

# Setup a loading step: fully prescribed (strain-controlled) arbitrary deformation
step = marmot.solvers.FiniteStrainSolver.Step()
step.timeStart = 0.0
step.timeEnd = 1.0
step.dTStart = 0.1
step.gradUIncrementTarget = np.array([[0.05, 0.02, 0.0], [0.0, -0.01, 0.0], [0.0, 0.0, -0.01]], dtype=np.float64)

step.isGradUComponentControlled = np.array([[True, True, True], [True, True, True], [True, True, True]])
step.isStressComponentControlled = np.logical_not(step.isGradUComponentControlled)
step.stressIncrementTarget = np.zeros((3, 3), dtype=np.float64)

solver.addStep(step)
solver.solve()

history = solver.getHistory()
marmot.testing.run_test(history, __file__)
