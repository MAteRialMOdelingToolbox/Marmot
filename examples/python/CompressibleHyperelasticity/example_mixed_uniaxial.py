import marmot
import marmot.testing
import numpy as np

print("Running example for FiniteStrain material: COMPRESSIBLEHYPERELASTICITY")

# Material properties extracted from C++ tests
# isochoric base 0 (NeoHooke), volumetric base 6 (DollSchweizerhof), mu, kappa, alpha, beta
properties = np.array([0.0, 6.0, 0.3, 1.0, 1.0, 3.0], dtype=np.float64)

# Setup solver
options = marmot.solvers.FiniteStrainSolver.SolverOptions()
solver = marmot.solvers.FiniteStrainSolver("COMPRESSIBLEHYPERELASTICITY", properties, options)

# Setup a loading step
step = marmot.solvers.FiniteStrainSolver.Step()
step.timeStart = 0.0
step.timeEnd = 1.0
step.dTStart = 0.1

step.isGradUComponentControlled = np.array([[True, True, True], [True, False, True], [True, True, False]], dtype=bool)
step.isStressComponentControlled = np.array(
    [[False, False, False], [False, True, False], [False, False, True]], dtype=bool
)

gradu_target = np.zeros((3, 3), dtype=np.float64)
gradu_target[0, 0] = 0.1
step.gradUIncrementTarget = gradu_target

stress_target = np.zeros((3, 3), dtype=np.float64)
step.stressIncrementTarget = stress_target

solver.addStep(step)
solver.solve()

history = solver.getHistory()
marmot.testing.run_test(history, __file__)
