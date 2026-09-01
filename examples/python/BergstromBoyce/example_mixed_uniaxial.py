import marmot
import marmot.testing
import numpy as np

print("Running example for FiniteStrain material: BERGSTROMBOYCE")

# Material properties: muA, kappaA, muB, kappaB, c1, c2, c3, implementationType (0 = CSDA)
properties = np.array([100.0, 1000.0, 50.0, 1000.0, 0.05, 1.0, 1.0, 0.0], dtype=np.float64)

# Setup solver
options = marmot.solvers.FiniteStrainSolver.SolverOptions()
solver = marmot.solvers.FiniteStrainSolver("BERGSTROMBOYCE", properties, options)

# Uniaxial stress-relaxation test: axial component gradU-controlled (ramp, then
# held constant), transverse components stress-free throughout.
is_gradU_controlled = np.array([[True, True, True], [True, False, True], [True, True, False]])
is_stress_controlled = np.logical_not(is_gradU_controlled)

ramp = marmot.solvers.FiniteStrainSolver.Step()
ramp.timeStart = 0.0
ramp.timeEnd = 1.0
ramp.dTStart = 0.1
ramp.isGradUComponentControlled = is_gradU_controlled
ramp.isStressComponentControlled = is_stress_controlled
gradu_target = np.zeros((3, 3), dtype=np.float64)
gradu_target[0, 0] = 0.1
ramp.gradUIncrementTarget = gradu_target
ramp.stressIncrementTarget = np.zeros((3, 3), dtype=np.float64)
solver.addStep(ramp)

hold = marmot.solvers.FiniteStrainSolver.Step()
hold.timeStart = 1.0
hold.timeEnd = 100.0
hold.dTStart = 1.0
hold.isGradUComponentControlled = is_gradU_controlled
hold.isStressComponentControlled = is_stress_controlled
hold.gradUIncrementTarget = np.zeros((3, 3), dtype=np.float64)
hold.stressIncrementTarget = np.zeros((3, 3), dtype=np.float64)
solver.addStep(hold)

solver.solve()

history = solver.getHistory()
marmot.testing.run_test(history, __file__)
