import marmot
import marmot.testing
import numpy as np

print("Running example for HypoElastic material: LINEARVISCOELASTICPOWERLAW")

# Material properties extracted from C++ tests
properties = np.array([2e5, 0.2, 0.5, 0.1, 10.0, 0.0001, 1.0], dtype=np.float64)

# Setup solver
options = marmot.solvers.HypoElasticSolver.SolverOptions()
solver = marmot.solvers.HypoElasticSolver("LINEARVISCOELASTICPOWERLAW", properties, options)

# Setup a loading step
step = marmot.solvers.HypoElasticSolver.Step()
step.timeStart = 0.0
step.timeEnd = 1.0
step.dTStart = 0.1
step.strainIncrementTarget = np.zeros(6, dtype=np.float64)
step.isStrainComponentControlled = np.array([False, False, False, True, True, True])
step.isStressComponentControlled = np.logical_not(step.isStrainComponentControlled)
s = np.zeros(6, dtype=np.float64)
s[0] = 10.0
step.stressIncrementTarget = s

solver.addStep(step)
solver.solve()

history = solver.getHistory()
marmot.testing.run_test(history, __file__)
