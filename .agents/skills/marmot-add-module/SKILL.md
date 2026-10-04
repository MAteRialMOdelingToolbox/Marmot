---
name: marmot-add-module
description: >-
  Universal workflow and architectural lifecycle for adding or extending any kind of module in Marmot.
  Use when implementing new materials, finite elements, core mechanical utilities, particles, or when routing to specialized skills.
---

# Module Development Lifecycle in Marmot

## 1. Skill Routing

| Task | Skill |
| :--- | :--- |
| Materials (small/finite strain, damage, nonlocal, AD) | [`marmot-add-material`](../marmot-add-material/SKILL.md) |
| Elements (templated `<nDim, nNodes>`, UL, micropolar) | [`marmot-add-element`](../marmot-add-element/SKILL.md) |
| Tests (`test.cmake` + analytical / tangent verification) | [`marmot-create-test`](../marmot-create-test/SKILL.md) |
| Documentation (Doxygen + Sphinx/Breathe feature page) | [`marmot-documentation`](../marmot-documentation/SKILL.md) |
| Code Review & QA Checklist | [`marmot-code-review`](../marmot-code-review/SKILL.md) |

---

## 2. Category Map & Layout

Scan order: `core` -> `materials` -> `elements` -> `particles` -> `materialpoints` -> `cells` -> `cellelements`.

```
modules/<category>/<ModuleName>/
├── module.cmake                  # marmot_add_module(<Name> REQUIRES ...)
├── test.cmake                    # ctest registration
├── include/Marmot/<Name>.h       # Public Doxygen header
├── src/<Name>.cpp                # Implementation
├── src/<Name>Registration.cpp    # Static factory registration
└── test/test.cpp                 # Standalone test executable
```

---

## 3. Implementation Checklist & Rules

1. **Minimal Changes & Zero Copies**: Modify code in place without unnecessary copies.
2. **Reuse Core Utilities**: Do not re-implement math, kinematics, or shape functions in `modules/core/`.
3. **Module CMake**:
   ```cmake
   marmot_add_module(<ModuleName>
       REQUIRES MarmotFiniteStrainMechanicsCore)
   ```
   `REQUIRES` lists every module whose headers the module includes, except modules already required (directly or transitively) by a listed module. The name must equal the directory name. An include of a module not
   required (directly or transitively) fails to compile. See `cmake/MarmotModules.cmake`.
4. **Registration**: Register with the appropriate factory in `src/<Name>Registration.cpp`.
5. **Testing**: Follow [`marmot-create-test`](../marmot-create-test/SKILL.md) (analytical benchmark + `MarmotMathCore` tangent test).
6. **Documentation**: Follow [`marmot-documentation`](../marmot-documentation/SKILL.md) (Doxygen + Sphinx `.rst`).
7. **QA**: Run [`marmot-code-review`](../marmot-code-review/SKILL.md) (`pre-commit run --all-files`, `ctest`).
