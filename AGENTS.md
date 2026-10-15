# AGENTS.md

Guidance for AI Agents working on Marmot (MAteRialMOdellingToolbox).

## What this is

Marmot is a high-performance C++20 shared library (`libMarmot`) providing finite elements and constitutive material models for quasi-brittle and engineering materials. Consumed by external FE frameworks ([EdelweissFE](https://github.com/MAteRialMOdelingToolbox/EdelweissFE), Abaqus, MOOSE, OpenSees) via lightweight interfaces (`Marmot::Materials::Material*` base classes, `Marmot::Discretizations::Elements::Element`, created through the `Marmot::Registry` factories).

Optional Python bindings (`modules/python`, built with [nanobind](https://github.com/wjakob/nanobind)) expose material point solvers (`marmot.solvers.HypoElasticSolver`, `marmot.solvers.FiniteStrainSolver`) for quick material testing and prototyping — see `doc/pages/python_bindings.md`.

Dependencies (located by `cmake/MarmotDependencies.cmake` as imported targets): Eigen (>= 3.3, incl. Eigen 5), autodiff and Fastor (CMake package, or a plain header search; no version enforced — CI pins autodiff v1.1.0 and Fastor V0.6.4).

## Build & Test

Always build out-of-source from `build/`:
```bash
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release && make -j$(nproc)
make install                                       # installs libMarmot and public headers
ctest --output-on-failure                          # run all tests
ctest -R TestLinearElastic                         # run specific test
```

To build and install the optional Python bindings (nanobind is fetched automatically at configure time; requires a Python 3.8+ dev environment):
```bash
cmake .. -DMARMOT_BUILD_PYTHON_BINDINGS=ON -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX
make install                                       # also installs the `marmot` python package
ctest -R Python_                                   # run the python solver example tests only
```

To measure code coverage locally (GCC/Clang only, requires `lcov`):
```bash
cmake .. -DCMAKE_BUILD_TYPE=Debug -DMARMOT_ENABLE_COVERAGE=ON && make -j$(nproc)
ctest --output-on-failure
lcov --directory . --capture --output-file coverage.info
lcov --remove coverage.info '/usr/*' '*/eigen/*' '*/autodiff/*' '*/Fastor/*' '*/build/*' --output-file coverage.filtered.info
genhtml coverage.filtered.info --output-directory coverage_html      # open coverage_html/index.html
```
CI runs the same instrumentation on every push/PR (`.github/workflows/coverage.yml`) and uploads results to Codecov.

## Namespace Convention

All code lives in the global `Marmot` namespace. The second level is a **closed set** — new second-level namespaces require maintainer approval:

| Namespace | Contents |
|---|---|
| `Marmot::Utilities` | Journal, MakeString, exceptions, state bookkeeping (`StateMapper`, `StateView`, `StateLayoutDynamic`, `StateVarVectorManager`) |
| `Marmot::Math` | Pure numerics: `Constants`, `Tensors` (all tensor types and functions, with `::FastorTensors`, `::FastorIndices` for einsum index symbols, `::TensorExponential`), `Differentiation` (`::Automatic`, `::Numerical`), `Integration`, plus root-level linear-algebra typedefs (`Vector6d`, `Matrix6d`, …), general math functions, eigen systems and basis functions |
| `Marmot::ContinuumMechanics` | Constitutive theory in descriptive topic sub-namespaces: `VoigtNotation`, `Elasticity`, `Kinematics`, `HaighWestergaard`, `CommonConstitutiveModels`, `Viscoelasticity`, `DeformationMeasures`, `StressMeasures`, `EnergyDensityFunctions`, `Plasticity`, `GradientDamage`, `PhaseField`, `LocalizationAnalysis`, `GeostaticStress`, `TimeIntegration`, `Algorithms` (substeppers/return mapping) |
| `Marmot::Materials` | Material base classes (`MaterialHypoElastic`, `MaterialFiniteStrain`, …), all concrete materials, material helpers (`KelvinChain`, `Wiechert`, `SolidificationTheory`), material-point solvers |
| `Marmot::Discretizations` | Everything spatial: `::Elements` (FE bases, shape functions, quadrature, concrete elements), `::Meshfree` (RKPM, particles), `::MaterialPoints`, `::Cells`. Root holds `MaterialSection`, `ElementProperties` |
| `Marmot::Registry` | All factories, flat — the external entry point (replaces the former `MarmotLibrary`) |
| `Marmot::Testing` | Test helpers only; production code never uses it |

Rules:

1. No classes/functions/variables in the global namespace or directly in `Marmot` (macros like `MARMOT_API` excepted).
2. Deeper nesting is welcome when it makes names more descriptive (up to 6 levels, e.g. `ContinuumMechanics::Plasticity::FlowIntegration::FirstOrder`); never a namespace that adds no information.
3. Derivative-order namespaces are spelled `FirstOrder` / `SecondOrder` / `ThirdOrder` only.
4. No namespace for a single class, and no namespace named after a class.
5. Class names never start with `Marmot` — the namespace says it. Header **filenames** keep the `Marmot` prefix (flat install dir).
6. Base classes at the root of their audience namespace; concrete implementations beside their base; factories in `Registry`.
7. Registration statics go in **unnamed namespaces** in `src/<Name>Registration.cpp`.
8. C++17 nested namespace syntax only (`namespace Marmot::ContinuumMechanics::Elasticity {`); namespaces are PascalCase.
9. No `using namespace` and no using-declarations at global or namespace scope in headers — function bodies and class scope only. (Sanctioned exception: `Marmot::Math::Differentiation::Automatic` imports `autodiff` inside its own leaf namespace.) `Marmot::Math::Tensors::FastorIndices` exists precisely for function-local `using namespace` in einsum-heavy code.

`scripts/checkNamespaceConventions.sh` greps for violations; CI treats findings as errors.

## General Coding Rules & Best Practices

- **Minimal Diffs**: Keep changes focused and strictly scoped; avoid reformatting unrelated files.
- **Zero-Copy & Memory Hygiene**: Avoid variable/tensor copies. Use `const&`, `Eigen::Map` (`mVector6d`, `mMatrix6d`), `std::span`, and modify `stateVars` strictly in place.
- **No Reinvention**: Reuse core math, tensor transformations, return mappers, and shape functions in `modules/core/`.
- **Reuse & Generalize**: Generalize templated formulations (`<int nDim, int nNodes>`) where applicable.
- **State Variable In-Place Writes**: Manage layout via `stateLayout.add("name", size)` / `finalize()`, access via `stateLayout.getAs<double&>(...)`. Elements slice state via `QPStateVarManager` / `MarmotStateVarVectorManager`.
- **Early Exit on Trivial Increments**: Check `if (dE.isZero(1e-14)) { dS_dE = Cel; return; }` to bypass return mapping.
- **Exceptions & Cutbacks**: Throw `MarmotExceptions.h` subclasses (`StressUpdateFailed`, `SolverConvergenceFailed`, `SolverTimestepExhausted`, `SolverIncrementsExhausted`) with `MakeString() << __PRETTY_FUNCTION__`.
- **Framework Logging**: Use `MarmotJournal::warningToMSG(...)` and `notificationToMSG(...)` instead of `std::cout`/`std::cerr`.
- **Test Every Feature**: Register standalone executables in `test.cmake`. Verify against **analytical solutions** and validate tangents via `MarmotMathCore`'s `numdiff`/`autodiff`.
- **Linting & Commits**: Format with `pre-commit run --all-files`. Commits follow [Conventional Commits](https://www.conventionalcommits.org). PRs target active `next_v<YY>.<MM>` branch.

## Workspace Skills (`.agents/skills/`)

- **[`marmot-add-module`](.agents/skills/marmot-add-module/SKILL.md)**: Universal lifecycle for adding/extending modules across categories and routing to specialized skills.
- **[`marmot-add-material`](.agents/skills/marmot-add-material/SKILL.md)**: Implementing small-strain, finite-strain, or nonlocal materials with in-place `stateVars` and factory registration.
- **[`marmot-add-element`](.agents/skills/marmot-add-element/SKILL.md)**: Implementing templated `<int nDim, int nNodes>` finite elements with `MarmotGeometryElement` and `MarmotElementFactory`.
- **[`marmot-create-test`](.agents/skills/marmot-create-test/SKILL.md)**: Constructing CTest suites validated against analytical solutions and `MarmotMathCore` `numdiff`/`autodiff` tangents.
- **[`marmot-documentation`](.agents/skills/marmot-documentation/SKILL.md)**: Authoring Doxygen API docstrings and Sphinx/Breathe feature pages (`python scripts/buildDocumentation.py`).
- **[`marmot-code-review`](.agents/skills/marmot-code-review/SKILL.md)**: QA review checklist covering C++20 standards, memory hygiene, factory registrations, and pre-commit checks.

## Architecture

### Module System & Category Scanning Order

Every `modules/<category>/<Name>/` with a `module.cmake` is auto-discovered (filters: `CORE_MODULES`, `MATERIAL_MODULES`, ... default `all`). Categories:
1. `core` — shared utilities (`MarmotMathCore`, `MarmotMechanicsCore`, `MarmotFiniteElementCore`, `MarmotUtilitiesCore`)
2. `materials` — constitutive models (bases `Marmot::Materials::{MaterialHypoElastic, MaterialFiniteStrain, MaterialGradientEnhancedHypoElastic}`)
3. `elements` — formulations (`Marmot::Discretizations::Elements::{Element, GeometryElement<nDim, nNodes>}`)
4. `particles` | 5. `materialpoints` | 6. `cells` | 7. `cellelements`

`modules/python` is *not* a category in this scan — it is a single, standalone nanobind extension (`modules/python/CMakeLists.txt`) added via `add_subdirectory` only when `MARMOT_BUILD_PYTHON_BINDINGS=ON`, wrapping the material point solvers on top of `libMarmot`. New materials/elements are registered through the C++ factories above; they do not need bindings added by hand.

Standard `module.cmake` (name = directory name; sources default to `src/*.cpp`):
```cmake
marmot_add_module(MyMaterial
    REQUIRES MarmotFiniteStrainMechanicsCore)
```
`REQUIRES` lists every module whose headers the module includes, except modules already required (directly or
transitively) by a listed module. Each module is an OBJECT library `Marmot_<Name>` that sees only its own and its
required modules' headers, so an undeclared dependency is a compile error. Without a `*_MODULES` filter every module
must be buildable: requiring a module that does not exist (misspelled, or its repository not checked out) fails
configuring. With a filter, a module whose requirement is filtered out or does not exist is skipped, unless it was
requested explicitly. All objects are assembled into the one `libMarmot`. Details: `cmake/MarmotModules.cmake`.

Standard `test.cmake`. A test sees its module's headers and `MarmotTesting.h`; `REQUIRES` lists further modules whose
headers it includes or whose materials/elements it creates by name through the factories (the test is skipped if
one is not built):
```cmake
SET(CURR_TEST_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/test")
add_marmot_test("TestMyModule" "${CURR_TEST_SOURCE_DIR}/test.cpp" [REQUIRES OtherModule])
```

### Self-Registering Factory Pattern

Materials and elements register statically in `src/<ModuleName>Registration.cpp`:

```cpp
// Material
#include "Marmot/MyMaterial.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"
namespace {
  using namespace Marmot::Materials;
  const static bool registered = Marmot::Registry::MaterialHypoElasticFactory::registerMaterial< MyMaterial >( "MYMATERIAL" );
}

// Templated Element
#include "Marmot/MyElement.h"
#include "Marmot/MarmotElementFactory.h"
namespace {
  using namespace Marmot::Discretizations::Elements;
  template < class T >
  Marmot::Registry::ElementFactory::elementFactoryFunction makeFactoryFunction() {
    return []( int id ) -> Element* { return new T( id ); };
  }
  const static bool CPS4_reg = Marmot::Registry::ElementFactory::registerElement( "CPS4", makeFactoryFunction< MyElement< 2, 4 > >() );
  const static bool C3D8_reg = Marmot::Registry::ElementFactory::registerElement( "C3D8", makeFactoryFunction< MyElement< 3, 8 > >() );
}
```
