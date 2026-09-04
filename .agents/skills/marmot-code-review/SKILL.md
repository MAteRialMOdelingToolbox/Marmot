---
name: marmot-code-review
description: >-
  Quality assurance, static check, and architectural review checklist for Marmot changes and pull requests.
  Use when reviewing code, checking pre-commit compliance, inspecting C++20 standards, or preparing PR submissions.
---

# Code Review & QA Checklist for Marmot

## 1. General Coding Rules & Performance
- [ ] **Minimal Code Changes**: Diffs are clean, focused, and scoped strictly to the task.
- [ ] **Zero-Copy & Memory Hygiene**: No copying of vectors, tensors, or state arrays. Memory is passed via `const&`, `std::span`, or mapped via `Eigen::Map` (`mVector6d`, `mMatrix6d`).
- [ ] **No Reinvention & Core Reuse**: Existing math, return mappers, Voigt conversions, and shape functions from `modules/core/` are reused.
- [ ] **In-Place State Variables**: Both materials (`stateLayout`) and elements (`QPStateVarManager`) write `stateVars` **strictly in place**.
- [ ] **Minimal State Variables**: Store **only strictly path-dependent** historical variables. Compute derived quantities on the fly.
- [ ] **Templated Elements**: Elements are templated on `< int nDim, int nNodes >` and use compile-time statically sized matrices from `MarmotGeometryElement< nDim, nNodes >`.

## 2. Formatting, Linting & Hygiene
```bash
pre-commit run --all-files
```
- [ ] **clang-format / python linters**: Pass with zero errors.
- [ ] **No Stray Artifacts**: Clean up temporary print statements or scratch debug dumps.

## 3. Tests & Documentation
- [ ] **Analytical Validation**: Behavior verified against **analytical solutions** (or high-precision benchmarks).
- [ ] **Tangent Verification**: Consistent tangent verified using `MarmotMathCore` numerical differentiation (`MarmotNumericalDifferentiation.h`) or automatic differentiation (`autodiff`).
- [ ] **CTest Registration**: Test executable registered via `add_marmot_test` in `test.cmake` and passes `ctest --output-on-failure`.
- [ ] **Doxygen & Sphinx**: Public APIs have Doxygen docstrings; feature page added in `doc/pages/features/` and builds with `python scripts/buildDocumentation.py`.

## 4. Commits & PR Targeting
- [ ] **Conventional Commits**: `<type>(<scope>): <summary>` (e.g. `feat(materials): add von Mises yield surface`).
- [ ] **PR Target Branch**: Targets active integration branch `next_v<YY>.<MM>` (e.g., `next_v26.11`).
- [ ] **Guidance Maintenance**: Update `AGENTS.md` or `.agents/skills/` if conventions change.
