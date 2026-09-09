import marmot
import marmot.testing
import numpy as np

print("Running example for FiniteStrain material: BERGSTROMBOYCE")

# Material properties: hyperelasticBase (0=NeoHooke), kappaA, kappaB, A1, A2, A3, B1, B2, B3, c1, c2, c3, implementationType (0 = CSDA)
properties = np.array([0.0, 1000.0, 1000.0, 100.0, 0.0, 0.0, 50.0, 0.0, 0.0, 0.05, 1.0, 1.0, 0.0], dtype=np.float64)

# Setup solver
options = marmot.solvers.FiniteStrainSolver.SolverOptions()
solver = marmot.solvers.FiniteStrainSolver("BERGSTROMBOYCE", properties, options)

# Setup a loading step: diagonal components stress-controlled, off-diagonal
# components gradU-controlled (held at zero, i.e. no shear).
step = marmot.solvers.FiniteStrainSolver.Step()
step.timeStart = 0.0
step.timeEnd = 1.0
step.dTStart = 0.1
step.gradUIncrementTarget = np.zeros((3, 3), dtype=np.float64)
step.isGradUComponentControlled = np.array([[False, True, True], [True, False, True], [True, True, False]])
step.isStressComponentControlled = np.logical_not(step.isGradUComponentControlled)
s = np.zeros((3, 3), dtype=np.float64)
s[0, 0] = 5.0
step.stressIncrementTarget = s

solver.addStep(step)
solver.solve()

history = solver.getHistory()
marmot.testing.run_test(history, __file__)
