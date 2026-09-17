import marmot
import marmot.testing
import numpy as np

print("Running example for FiniteStrain material: COMPRESSIBLENEOHOOKE")

# Material properties extracted from C++ tests
properties = np.array([1.0, 0.3], dtype=np.float64)

# Setup solver
options = marmot.solvers.FiniteStrainSolver.SolverOptions()
solver = marmot.solvers.FiniteStrainSolver("COMPRESSIBLENEOHOOKE", properties, options)

# Setup a loading step
step = marmot.solvers.FiniteStrainSolver.Step()
step.timeStart = 0.0
step.timeEnd = 1.0
step.dTStart = 0.1
step.gradUIncrementTarget = np.zeros((3, 3), dtype=np.float64)
step.isGradUComponentControlled = np.array([[False, True, True], [True, False, True], [True, True, False]])
step.isStressComponentControlled = np.logical_not(step.isGradUComponentControlled)
s = np.zeros((3, 3), dtype=np.float64)
s[0, 0] = 10.0
step.stressIncrementTarget = s

solver.addStep(step)
solver.solve()

history = solver.getHistory()
marmot.testing.run_test(history, __file__)
