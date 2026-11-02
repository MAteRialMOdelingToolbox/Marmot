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

step.isStrainComponentControlled = np.array([True, False, False, True, True, True], dtype=bool)
step.isStressComponentControlled = np.array([False, True, True, False, False, False], dtype=bool)

strain_target = np.zeros(6, dtype=np.float64)
strain_target[0] = 0.1
step.strainIncrementTarget = strain_target

stress_target = np.zeros(6, dtype=np.float64)
step.stressIncrementTarget = stress_target

solver.addStep(step)
solver.solve()

history = solver.getHistory()
marmot.testing.run_test(history, __file__)
