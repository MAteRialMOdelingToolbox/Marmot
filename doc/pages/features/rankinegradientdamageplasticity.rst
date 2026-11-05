.. _rankineGradientDamagePlasticity:

Gradient-enhanced Rankine damage-plasticity model
=================================================

Implicit-gradient version of the :ref:`rankineDamagePlasticity`: multi-surface Rankine plasticity in the effective
stress with isotropic damage driven by a nonlocal variable, for the general gradient-enhanced displacement element.

.. list-table::
   :header-rows: 1
   :align: left

   * - **Index**
     - **Model Parameter**
     - **Description**
   * - 0
     - :math:`E`
     - Young’s modulus
   * - 1
     - :math:`\nu`
     - Poisson’s ratio
   * - 2
     - :math:`f_\mathrm{t}`
     - Uniaxial tensile strength
   * - 3
     - :math:`\varepsilon_\mathrm{f}`
     - Softening parameter of the damage law (positive)
   * - 4
     - :math:`l`
     - Nonlocal length (nonnegative)
   * - 5
     - :math:`m`
     - Over-nonlocal weight (0: local, 1: nonlocal, > 1: over-nonlocal)
   * - 6
     - :math:`\rho`
     - Density (optional)
   * - 7
     - :math:`\eta`
     - Nonlocal viscosity (optional, default 0; explicit dynamics only)

.. list-table::
   :header-rows: 1
   :align: left

   * - **State Variable**
     - **Name**
     - **Description**
   * - :math:`\kappa_\mathrm{d}`
     - ``kappaD``
     - Local damage driving variable (volumetric plastic strain)
   * -
     - ``acYield``
     - Yield flag: 1 for a plastic, 0 for an elastic increment
   * - :math:`\omega`
     - ``omega``
     - Damage
   * - :math:`\bar\kappa_\mathrm{d}`
     - ``kappaDNonlocal``
     - Nonlocal damage driving variable
   * - :math:`\hat\kappa_\mathrm{max}`
     - ``kappaDWeightedMax``
     - History maximum of the over-nonlocal combination

Theory
------

The local damage driving variable :math:`\kappa_\mathrm{d}`, with
:math:`\Delta\kappa_\mathrm{d} = \sum_A \Delta\lambda_A` as in the :ref:`rankineDamagePlasticity`, is the source of the
Helmholtz equation solved by the element,

.. math:: \bar\kappa_\mathrm{d} - l^2\,\Delta\bar\kappa_\mathrm{d} = \kappa_\mathrm{d}, \qquad
   \nabla\bar\kappa_\mathrm{d}\cdot\boldsymbol{n} = 0 \text{ on the boundary}.

The damage follows from the over-nonlocal combination

.. math:: \hat\kappa = m\,\bar\kappa_\mathrm{d} + (1-m)\,\kappa_\mathrm{d}, \qquad
   \omega = 1 - \exp\left(-\hat\kappa_\mathrm{max}/\varepsilon_\mathrm{f}\right), \qquad
   \hat\kappa_\mathrm{max} = \max_{\tau \le t} \hat\kappa(\tau), \qquad
   \boldsymbol{\sigma} = (1-\omega)\,\bar{\boldsymbol{\sigma}}.

The history maximum keeps the damage irreversible also for :math:`m > 1`, where :math:`\hat\kappa` may decrease.
Since the effective stress does not depend on the damage, :math:`\kappa_\mathrm{d}` does not depend on
:math:`\bar\kappa_\mathrm{d}`.

Consistent Algorithmic Tangent Operator
---------------------------------------

.. math:: \frac{\partial\boldsymbol{\sigma}}{\partial\Delta\boldsymbol{\varepsilon}}
   = (1-\omega)\,\bar{\mathbb{C}}^\mathrm{ep} - (1-m)\,\omega'\,\bar{\boldsymbol{\sigma}} \otimes
   \frac{\partial\kappa_\mathrm{d}}{\partial\Delta\boldsymbol{\varepsilon}}, \qquad
   \frac{\partial\boldsymbol{\sigma}}{\partial\bar\kappa_\mathrm{d}} = -m\,\omega'\,\bar{\boldsymbol{\sigma}}, \qquad
   \frac{\partial\kappa_\mathrm{d}}{\partial\bar\kappa_\mathrm{d}} = 0,

with :math:`\omega' = \exp\left(-\hat\kappa/\varepsilon_\mathrm{f}\right)/\varepsilon_\mathrm{f}` while
:math:`\hat\kappa` exceeds its history maximum, and zero otherwise.

The model is registered as ``RANKINEGRADIENTDAMAGEPLASTICITY`` in the gradient-enhanced hypoelastic factory with one
nonlocal field.

.. doxygenclass:: Marmot::Materials::RankineGradientDamagePlasticity
   :allow-dot-graphs:
