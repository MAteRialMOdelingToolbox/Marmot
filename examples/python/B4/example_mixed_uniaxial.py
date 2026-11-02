import marmot
import marmot.testing
import numpy as np

print("Running example for HypoElastic material: B4")

# Material properties extracted from C++ tests
properties = np.array(
    [
        0.2,
        20.6,
        91.0,
        4.8,
        5.9,
        0.1,
        0.5,
        12.0,
        1e-5,
        0.0,
        3.0,
        1.45,
        -4.5,
        -0.0015,
        90.0,
        7.0,
        1.0,
        400.0,
        11.0,
        2e-4,
        -100.0,
        1.0,
    ],
    dtype=np.float64,
)

# Setup solver
options = marmot.solvers.HypoElasticSolver.SolverOptions()
solver = marmot.solvers.HypoElasticSolver("B4", properties, options)

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
