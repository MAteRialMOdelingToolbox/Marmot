import marmot
import marmot.testing
import numpy as np

print("Running example for FiniteStrain material: FINITESTRAINORTHOTROPICBIOTVISCOELASTICITY")

# Material properties extracted from C++ tests
properties = np.array(
    [20000.0, 10000.0, 15000.0, 0.25, 0.35, 0.30, 4000.0, 6000.0, 5000.0, 1.0, 0.3, 10.0], dtype=np.float64
)

# Setup solver
options = marmot.solvers.FiniteStrainSolver.SolverOptions()
solver = marmot.solvers.FiniteStrainSolver("FINITESTRAINORTHOTROPICBIOTVISCOELASTICITY", properties, options)

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
