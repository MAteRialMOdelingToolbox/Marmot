# AGENTS.md

Guidance for AI Agents working on Marmot (MAteRialMOdellingToolbox).

## What this is

Marmot is a high-performance C++20 shared library (`libMarmot`) providing finite elements and constitutive material models for quasi-brittle and engineering materials. Consumed by external FE frameworks ([EdelweissFE](https://github.com/MAteRialMOdelingToolbox/EdelweissFE), Abaqus, MOOSE, OpenSees) via lightweight interfaces (`MarmotMaterial*`, `MarmotElement`).

Optional Python bindings (`modules/python`, built with [nanobind](https://github.com/wjakob/nanobind)) expose material point solvers (`marmot.solvers.HypoElasticSolver`, `marmot.solvers.FiniteStrainSolver`) for quick material testing and prototyping — see `doc/pages/python_bindings.md`.

Dependencies (discovered via top-level `CMakeLists.txt`): Eigen (`find_package(Eigen3 3.3 REQUIRED)`), autodiff and Fastor (located via `find_path()`, no version enforced — CI pins autodiff v1.1.0 and Fastor V0.6.4).

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

Categories are auto-discovered in strict dependency order (modules only depend on earlier categories):
1. `core` — shared utilities (`MarmotMathCore`, `MarmotMechanicsCore`, `MarmotFiniteElementCore`, `MarmotUtilitiesCore`)
2. `materials` — constitutive models (`MarmotMaterialHypoElastic`, `MarmotMaterialFiniteStrain`, `MarmotMaterialGeneralGradientEnhancedHypoElastic`)
3. `elements` — formulations (`MarmotElement`, `MarmotGeometryElement<nDim, nNodes>`)
4. `particles` | 5. `materialpoints` | 6. `cells` | 7. `cellelements`

`modules/python` is *not* a category in this scan — it is a single, standalone nanobind extension (`modules/python/CMakeLists.txt`) added via `add_subdirectory` only when `MARMOT_BUILD_PYTHON_BINDINGS=ON`, wrapping the material point solvers on top of `libMarmot`. New materials/elements are registered through the C++ factories above; they do not need bindings added by hand.

Standard `module.cmake`:
```cmake
list(APPEND INSTALLED_MODULE_INCLUDE_DIRS "${CMAKE_CURRENT_LIST_DIR}/include")
file(GLOB module_sources CONFIGURE_DEPENDS "${CMAKE_CURRENT_LIST_DIR}/src/*.cpp")
list(APPEND sources ${module_sources})
```
Standard `test.cmake`:
```cmake
SET(CURR_TEST_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/test")
add_marmot_test("TestMyModule" "${CURR_TEST_SOURCE_DIR}/test.cpp")
```

### Self-Registering Factory Pattern

Materials and elements register statically in `src/<ModuleName>Registration.cpp`:

```cpp
// Material
#include "Marmot/MyMaterial.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"
namespace Marmot::Materials::Registration {
  using namespace MarmotLibrary;
  const static bool registered = MarmotMaterialHypoElasticFactory::registerMaterial< MyMaterial >( "MYMATERIAL" );
}

// Templated Element
#include "Marmot/MyElement.h"
#include "Marmot/MarmotElementFactory.h"
namespace Marmot::Elements::Registration {
  using namespace MarmotLibrary;
  template < class T >
  MarmotElementFactory::elementFactoryFunction makeFactoryFunction() {
    return []( int id ) -> MarmotElement* { return new T( id ); };
  }
  const static bool CPS4_reg = MarmotElementFactory::registerElement( "CPS4", makeFactoryFunction< MyElement< 2, 4 > >() );
  const static bool C3D8_reg = MarmotElementFactory::registerElement( "C3D8", makeFactoryFunction< MyElement< 3, 8 > >() );
}
```
