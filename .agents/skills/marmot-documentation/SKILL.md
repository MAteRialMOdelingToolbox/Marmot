---
name: marmot-documentation
description: >-
  Procedure for writing, updating, structuring, and building documentation in Marmot.
  Use when adding or modifying Doxygen comments in C++ headers, Sphinx/Breathe docs (.rst files), or theory manuals.
---

# Documentation Guide for Marmot

Marmot compiles API docs via **Doxygen** (XML) and website pages via **Sphinx** + **Breathe**.

---

## 1. Doxygen in C++ Headers
Document public classes and functions using `@brief`, `@param[in/out]`, `@return`, `@tparam`, `@details`, and LaTeX math (`\f[ ... \f]` or `\f$ ... \f$`):

```cpp
/**
 * @brief Computes Cauchy stress and algorithmic consistent tangent.
 * @param[in,out] state State instance holding stress, energy, and in-place state variables.
 * @param[out] dStress_dStrain Algorithmic consistent tangent matrix \f$\mathbf{C}\f$.
 * @param[in] dStrain Linearized strain increment \f$\Delta\boldsymbol{\varepsilon}\f$.
 * @param[in] timeInfo Current time and time increment.
 */
void computeStress( state3D& state, Matrix6d& dStress_dStrain, const Vector6d& dStrain, const timeInfo& timeInfo ) const;
```

---

## 2. Sphinx Feature Page (`doc/pages/features/<name>.rst`)

```rst
.. _mymaterial:

My Material Model
=================

Governing equations:

.. math::
   \boldsymbol{\sigma} = \mathbf{C} : (\boldsymbol{\varepsilon} - \boldsymbol{\varepsilon}^p)

Parameters:
- ``E``: Young's modulus (:math:`\text{MPa}`).
- ``nu``: Poisson's ratio (-).

.. doxygenclass:: Marmot::Materials::MyMaterial
   :project: Marmot
   :members:
```

Index the page in `doc/pages/features/materials.rst` or `doc/pages/features/elements.rst`.

---

## 3. Local Build Command
```bash
python scripts/buildDocumentation.py
```
Output is generated at `doc/doc_out/sphinx/index.html`. Ensure zero Doxygen or Sphinx warnings.
