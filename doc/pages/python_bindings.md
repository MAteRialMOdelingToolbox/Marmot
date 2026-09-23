# Python Bindings

`Marmot` provides high-performance Python bindings built with [nanobind](https://github.com/wjakob/nanobind).
These bindings expose Marmot's material point solvers directly in Python, enabling quick material testing, parameter studies, unit response verification, and interactive prototyping with NumPy arrays.

---

## Installation & Setup

Building the Python bindings is optional and requires a Python 3.8+ development environment.
Nanobind is automatically fetched at CMake configuration time.

To build and install Marmot with Python bindings:

```bash
mkdir -p build && cd build
cmake -DMARMOT_BUILD_PYTHON_BINDINGS=ON -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX ..
make install
```

When building inside a build tree without installing globally, make sure your `PYTHONPATH` includes the build's python folder:

```bash
export PYTHONPATH=/path/to/marmot/build/python:$PYTHONPATH
```

---

## Module Overview

The Python package is imported as `marmot`:

```python
import marmot
```

The main solvers are available under the `marmot.solvers` submodule:

- `marmot.solvers.HypoElasticSolver`: Material point solver for small-strain / hypo-elastic constitutive laws.
- `marmot.solvers.FiniteStrainSolver`: Material point solver for finite-strain constitutive models formulated with the deformation gradient.

---

## Quickstart Examples

### 1. Small-Strain Hypo-Elastic Solver

In small-strain problems, strains and stresses are represented in 6D Voigt notation:
`[11, 22, 33, 23, 13, 12]` (or `[xx, yy, zz, yz, xz, xy]`).

```python
import numpy as np
import marmot

# 1. Define material properties (e.g. Young's modulus E=20000, Poisson's ratio nu=0.25)
properties = np.array([20000.0, 0.25], dtype=np.float64)

# 2. Instantiate solver with material model name and properties
solver = marmot.solvers.HypoElasticSolver("LINEARELASTIC", properties)

# 3. Create and configure a loading step
step = marmot.solvers.HypoElasticSolver.Step()
step.timeStart = 0.0
step.timeEnd = 1.0
step.dTStart = 0.1

# Strain-controlled uniaxial stretch in 11 direction
step.isStrainComponentControlled = np.array([True, True, True, True, True, True])
step.isStressComponentControlled = np.logical_not(step.isStrainComponentControlled)
step.strainIncrementTarget = np.array([0.01, 0.0, 0.0, 0.0, 0.0, 0.0], dtype=np.float64)
step.stressIncrementTarget = np.zeros(6, dtype=np.float64)

# 4. Add step and solve
solver.addStep(step)
solver.solve()

# 5. Extract results
history = solver.getHistory()
for entry in history:
    print(f"Time: {entry.time:.2f}, Stress[0]: {entry.stress[0]:.4f}, Strain[0]: {entry.strain[0]:.4f}")
```

### 2. Finite-Strain Solver

In finite-strain problems, deformation measures (displacement gradient $\nabla \mathbf{u}$, deformation gradient $\mathbf{F}$) and stresses (Kirchhoff / Cauchy stresses) are 3x3 matrices:

```python
import numpy as np
import marmot

# 1. Material properties (e.g. Compressible Neo-Hooke: mu=1.0, nu=0.3)
properties = np.array([1.0, 0.3], dtype=np.float64)

# 2. Instantiate FiniteStrainSolver
solver = marmot.solvers.FiniteStrainSolver("COMPRESSIBLENEOHOOKE", properties)

# 3. Define loading step
step = marmot.solvers.FiniteStrainSolver.Step()
step.timeStart = 0.0
step.timeEnd = 1.0
step.dTStart = 0.1

# Mixed control: prescribe displacement gradient component [0,0], stress free in transverse directions
step.isGradUComponentControlled = np.array([
    [True,  True,  True ],
    [True,  False, True ],
    [True,  True,  False]
], dtype=bool)
step.isStressComponentControlled = np.array([
    [False, False, False],
    [False, True,  False],
    [False, False, True ]
], dtype=bool)

gradu_target = np.zeros((3, 3), dtype=np.float64)
gradu_target[0, 0] = 0.1
step.gradUIncrementTarget = gradu_target
step.stressIncrementTarget = np.zeros((3, 3), dtype=np.float64)

# 4. Solve
solver.addStep(step)
solver.solve()

# 5. Inspect history
history = solver.getHistory()
print(f"Recorded {len(history)} increments.")
last_entry = history[-1]
print("Final Deformation Gradient F:\n", last_entry.F)
print("Final Stress:\n", last_entry.stress)
```

---

## API Reference

### `marmot.solvers.HypoElasticSolver`

#### Constructor
```python
HypoElasticSolver(name: str, props: np.ndarray, opts: SolverOptions = SolverOptions())
```
- `name`: String name of the registered hypo-elastic material model (e.g., `"LINEARELASTIC"`, `"VONMISES"`, `"B4"`).
- `props`: 1D continuous NumPy array of material parameters.
- `opts`: Optional `SolverOptions` configuration.

#### Methods
- `addStep(step: Step)`: Add a loading step.
- `getSteps() -> list[Step]`: Return list of configured steps.
- `clearSteps()`: Clear all configured steps.
- `setInitialState(initialStress: np.ndarray, initialStateVars: np.ndarray)`: Set initial 6D stress vector and internal state variable vector.
- `resetToInitialState()`: Reset solver to initial state.
- `getNumberOfStateVariables() -> int`: Get count of required internal state variables.
- `solve()`: Execute simulation for all added steps.
- `getHistory() -> list[HistoryEntry]`: Retrieve output history sequence.
- `clearHistory()`: Clear history buffer.
- `printHistory()`: Print recorded history to standard output.
- `exportHistoryToCSV(filename: str)`: Export history records to a CSV file.

#### Nested Classes
- **`Step`**: Configures step duration and boundary conditions.
  - `timeStart`, `timeEnd`, `dTStart`, `dTMin`, `dTMax`: Time integration bounds.
  - `strainIncrementTarget`: Target 6D strain increment array.
  - `stressIncrementTarget`: Target 6D stress increment array.
  - `isStrainComponentControlled`: 6-element boolean mask for strain-controlled components.
  - `isStressComponentControlled`: 6-element boolean mask for stress-controlled components.
  - `maxIncrements`: Increment count limit.
  - `checkControl()`: Validates that every component is uniquely controlled.
- **`SolverOptions`**:
  - `maxIterations`: Maximum Newton-Raphson iterations per increment (default: `25`).
  - `residualTolerance`: Residual convergence tolerance (default: `1e-10`).
  - `correctionTolerance`: Correction norm tolerance (default: `1e-10`).
- **`HistoryEntry`**:
  - `time`: Time value at increment end.
  - `stress`: 6D stress vector.
  - `strain`: 6D strain vector.
  - `dStressdStrain`: 6x6 material tangent matrix.
  - `stateVars`: Internal state variables vector.
  - `print()`: Print entry details.

---

### `marmot.solvers.FiniteStrainSolver`

#### Constructor
```python
FiniteStrainSolver(name: str, props: np.ndarray, opts: SolverOptions = SolverOptions())
```

#### Nested Classes
- **`Step`**:
  - `timeStart`, `timeEnd`, `dTStart`, `dTMin`, `dTMax`: Time stepping controls.
  - `gradUIncrementTarget`: 3x3 displacement gradient increment matrix $\Delta \nabla \mathbf{u}$.
  - `stressIncrementTarget`: 3x3 stress increment target matrix.
  - `isGradUComponentControlled`: 3x3 boolean mask.
  - `isStressComponentControlled`: 3x3 boolean mask.
  - `maxIncrements`: Max increment count.
  - `checkControl()`: Validates component control flags.
- **`HistoryEntry`**:
  - `time`: Increment timestamp.
  - `stress`: 3x3 stress tensor.
  - `F`: 3x3 deformation gradient tensor $\mathbf{F}$.
  - `dTau_dF`: 3x3x3x3 material tangent tensor $\partial \boldsymbol{\tau} / \partial \mathbf{F}$.
  - `stateVars`: 1D internal state variable array.

---

## Testing Utilities (`marmot.testing`)

For automated testing and continuous integration, `marmot.testing` provides helpers to compare and generate `.npz` reference outputs:

```python
import marmot.testing

# Automatically checks against <script_name>_reference.npz or writes reference with --writeReferenceSolution
marmot.testing.run_test(history, __file__)
```
