# Contributing to Marmot

We welcome pull requests of all sizes and issues with clear, minimal reproductions. This guide explains how to set up your environment, follow our commit conventions, run tests, and open high‑quality issues and PRs.

---

## Code of Conduct
By participating in this project you agree to uphold a respectful, inclusive environment. Be kind, be constructive, and assume good intent. If you see a problem, please open an issue.

## Ways to Contribute
- Report bugs and propose enhancements via GitHub Issues.
- Improve documentation, examples, and comments.
- Add tests, improve performance, or refactor internals.
- Implement new features aligned with the project scope.

> **Pull requests are welcome!** See [Pull Requests](#pull-requests).

---

## Pre-commit Hooks
We use **pre-commit** to run formatting and static checks before each commit to keep the codebase clean and consistent.

### Install and Run
```bash
# Install the pre-commit framework (one-time)
pip install pre-commit

# Install the repo's hooks
pre-commit install

# (Recommended) Update hooks to the latest pinned versions
pre-commit autoupdate

# Run hooks against all files (useful before pushing)
pre-commit run --all-files
```
> If any hook modifies files (e.g., formatting), re-stage and re-commit.

> See `.pre-commit-config.yaml` at the repo root for the authoritative list of hooks we use.

---

## Conventional Commits
All commits **must** follow the [Conventional Commits](https://www.conventionalcommits.org) specification. This makes history easier to read and enables automated changelogs.

### Types
Use one of the following types (with an optional **scope** in parentheses):
- **feat**: A new feature
- **fix**: A bug fix
- **docs**: Documentation-only changes
- **test**: Adding or updating tests
- **refactor**: Code change that neither fixes a bug nor adds a feature
- **perf**: Performance improvement
- **style**: Formatting, missing semicolons, etc.; no code change
- **build**: Build system or external dependencies
- **ci**: CI configuration or scripts
- **chore**: Maintenance tasks (no src or test changes)
- **revert**: Reverts a previous commit

### Examples
```
feat(core): add consistent tangent operator for large strain model
fix(io): correct VTK export for element-wise results
docs(readme): clarify build badges and usage examples
test(solver): add ctest for Newton convergence edge cases
refactor(mesh): split node utilities into dedicated module
```
> Keep commit messages concise in the subject line (≤72 chars). Use the body to explain **what** and **why**, and reference issues like `Fixes #123`.

---

## Opening Issues
Before filing a new issue, please search to avoid duplicates. When reporting bugs, include:
- **Environment**: OS, compiler & version, CMake version
- **Exact steps** to reproduce
- **Expected vs actual** behavior
- **Logs / error messages** (as text, not screenshots, when possible)
- **Minimal Reproducible Example (MRE)**

---

## Pull Requests
We use GitHub flow: fork → branch → PR → review → merge.

1. **Fork** the repository and create a feature branch from `master`:
   ```bash
   git checkout -b feat/<short-scope>-<concise-topic>
   ```
2. Make changes, commit using [Conventional Commits](#conventional-commits), and ensure the code builds locally.
3. **Run pre-commit** on all files and **run tests** with `ctest` (see below).
4. **Open a PR** with a clear title and description. Link any related issues.

### Feature PRs: Documentation & Tests
Feature PRs **must** include:
- **Documentation** updates (README, docs site, or inline comments) explaining the new behavior and public APIs.
- **Automated tests** registered with **ctest** that cover the new behavior and edge cases.

### PR Checklist
- [ ] My PR title follow Conventional Commits
- [ ] `pre-commit run --all-files` passes locally
- [ ] The project **builds** with CMake
- [ ] All **tests pass** locally via `ctest`
- [ ] New/changed behavior is **documented**
- [ ] I added/updated **ctest** tests for features/bug fixes

---

## Documentation using Doxygen/Sphinx

Marmot uses **Doxygen** for API-level documentation and **Sphinx** for generating the full developer and theory documentation website. Contributors are expected to document **all** code changes thoroughly.

### Doxygen Conventions
- Every **class**, **function**, **structure**, **enum**, and **member variable** must include Doxygen-style comments.
- Use the standard Doxygen tags:
  - `@brief` — short one-line summary
  - `@param` — for each function parameter
  - `@return` — for return values
  - `@tparam` — for template parameters
  - `@details` — for details on equations, theory, etc.
  - `@note`, `@warning`, `@todo` — for relevant annotations

Example for a `void` function:
```cpp
/**
 * @brief Computes something from ...
 * @param[in] input The description of input.
 * @param[inout] something The description of the input/output parameter.
 * @param[in] parameters The description of the parameters.
 */
void computeSomething(const someType& input, otherType& something, const parameterType&  parameters);
```

> Consistent and complete documentation is **mandatory** for all public APIs and newly introduced code.

### Local documentation build

The documentation is built using Doxygen and Sphinx/Breathe following the commands:
```bash
   # install the requirements
   conda install --file doc/requirements.txt
   # build the documentation
   python scripts/buildDocumentation.py
```

### Module Documentation for New Models and Elements
If you add a **new material model** or **finite element**, you must:
1. Add a corresponding `module.rst` file under the appropriate Sphinx module directory (`docs/pages/features/`).
2. Include:
   - Theoretical background
   - Governing equations and assumptions
   - Parameters and usage examples
   - References, if applicable
3. Register it in the corresponding `materials.rst` or `elements.rst`.

> Documentation completeness will be part of PR review for new physical models or elements.

---

## Adding & Registering Tests
If you are adding tests, please:
1. Place sources under the modules' test directory, i.e., under `modules/my-specific-module/tests/`.
2. Each test source file must include a `main` function which executes all tests.
3. To make your life easier for writing tests, look for convenience functions in `include/Marmot/MarmotTesting.h` to be used for testing.
4. Register the test in CMake in the modules' `test.cmake`.
   ```cmake
   # modules/my-specific-module/test.cmake
   add_marmot_test("TestMarmotMyNewFeature" "${CURR_TEST_SOURCE_DIR}/TestMarmotMyNewFeature.cpp")
   ```
5. Ensure tests are **deterministic** and run quickly.
6. Verify with:
   ```bash
   cd build
   ctest --output-on-failure
   ```
---

## Style & Guidelines
- Prefer readable, well‑documented code over clever one‑liners.
- Include comments for non‑obvious logic and derivations (e.g., equations, references).
- Keep functions small and cohesive; prefer self‑documenting names.
- Avoid introducing new external dependencies unless justified.

> Formatting and basic static checks are enforced via **pre-commit** hooks.

---

## License
By contributing, you agree that your contributions will be licensed under the repository’s license. See `LICENSE` for details.
