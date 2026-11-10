import marmot
import marmot.testing
import numpy as np

print("Running example for FiniteStrain material: COMPRESSIBLEFINITESTRAINFRACTIONALSCHAPERY")

# Material properties: [hyperelasticBase=0 (NeoHooke), onlyShearCreep=1.0, mu=100.0, kappa=300.0, alpha=0.3, tau0=5.0, g_inf=0.2, a1=0.1, a2=0.05, b_S=0.02]
properties = np.array([0.0, 1.0, 100.0, 300.0, 0.3, 5.0, 0.2, 0.1, 0.05, 0.02], dtype=np.float64)

# Setup solver
options = marmot.solvers.FiniteStrainSolver.SolverOptions()
solver = marmot.solvers.FiniteStrainSolver("COMPRESSIBLEFINITESTRAINFRACTIONALSCHAPERY", properties, options)

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
gradu_target[0, 0] = 0.05
step.gradUIncrementTarget = gradu_target

stress_target = np.zeros((3, 3), dtype=np.float64)
step.stressIncrementTarget = stress_target

solver.addStep(step)
solver.solve()

history = solver.getHistory()
marmot.testing.run_test(history, __file__)
