.. _rankineDamagePlasticity:

Rankine damage-plasticity model
===============================

Multi-surface Rankine plasticity in the effective stress combined with isotropic damage, following the damage-plastic
model of Grassl and Jirásek (2006).
The effective stress is bounded by the perfectly plastic :ref:`rankineSofteningPlasticity`; softening results from the
damage only.

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
     - :math:`\omega_\mathrm{max}`
     - Damage threshold, the maximum damage :math:`0 < \omega_\mathrm{max} \le 1` (1: unbounded exponential law)
   * - 5
     - :math:`\rho`
     - Density (optional)

.. list-table::
   :header-rows: 1
   :align: left

   * - **State Variable**
     - **Name**
     - **Description**
   * - :math:`\kappa_\mathrm{d}`
     - ``kappaD``
     - Damage driving variable (volumetric plastic strain)
   * -
     - ``acYield``
     - Yield flag: 1 for a plastic, 0 for an elastic increment
   * - :math:`\omega`
     - ``omega``
     - Damage

Theory
------

The effective stress :math:`\bar{\boldsymbol{\sigma}}` follows from perfect Rankine plasticity with the yield
functions :math:`f_A = \bar\sigma_A - f_\mathrm{t} \leq 0` and associated flow.
The nominal stress and the damage read

.. math:: \boldsymbol{\sigma} = (1-\omega)\,\bar{\boldsymbol{\sigma}}, \qquad
   \omega = \omega_\mathrm{max}\left[1 - \exp\left(-\frac{\kappa_\mathrm{d}}{\omega_\mathrm{max}\,
   \varepsilon_\mathrm{f}}\right)\right], \qquad
   \Delta\kappa_\mathrm{d} = \mathrm{tr}\,\Delta\boldsymbol{\varepsilon}^\mathrm{p} = \sum_A \Delta\lambda_A.

The damage threshold :math:`\omega_\mathrm{max}` bounds the damage without changing its initial slope
:math:`1/\varepsilon_\mathrm{f}`; :math:`\omega_\mathrm{max} = 1` gives the exponential law
:math:`\omega = 1 - \exp\left(-\kappa_\mathrm{d}/\varepsilon_\mathrm{f}\right)`.
All principal plastic strain increments are nonnegative for the Rankine surfaces, hence :math:`\kappa_\mathrm{d}` and
:math:`\omega` never decrease, and under uniaxial tension :math:`\sigma = (1-\omega)\,f_\mathrm{t}` with
:math:`\kappa_\mathrm{d} = \varepsilon^\mathrm{p}`, which tends to the residual strength
:math:`(1-\omega_\mathrm{max})\,f_\mathrm{t}`.
:math:`\varepsilon_\mathrm{f}` is a material parameter without regularization by the element size; for regularized
softening use the :ref:`rankineGradientDamagePlasticity`.

The free energy is :math:`\psi = (1-\omega)\,\bar\psi_\mathrm{el}` with
:math:`\bar\psi_\mathrm{el} = \frac12\bar{\boldsymbol{\sigma}}:\mathbb{C}^{-1}:\bar{\boldsymbol{\sigma}}`, and the
dissipation increment is
:math:`(1-\omega)\,\bar{\boldsymbol{\sigma}}:\Delta\boldsymbol{\varepsilon}^\mathrm{p} + \bar\psi_\mathrm{el}\,\Delta\omega`.

Consistent Algorithmic Tangent Operator
---------------------------------------

.. math:: \frac{\partial\boldsymbol{\sigma}}{\partial\Delta\boldsymbol{\varepsilon}}
   = (1-\omega)\,\bar{\mathbb{C}}^\mathrm{ep}
   - \bar{\boldsymbol{\sigma}} \otimes \frac{\partial\omega}{\partial\Delta\boldsymbol{\varepsilon}},
   \qquad
   \frac{\partial\kappa_\mathrm{d}}{\partial\Delta\boldsymbol{\varepsilon}}
   = \boldsymbol{1} - \bar{\mathbb{C}}^{\mathrm{ep},T}:\mathbb{C}^{-1}:\boldsymbol{1},

with the algorithmic tangent :math:`\bar{\mathbb{C}}^\mathrm{ep}` of the effective stress.

The model is registered as ``RANKINEDAMAGEPLASTICITY`` and, with the co-rotational :ref:`hughesWingetWrapper` wrapper, as
``RANKINEDAMAGEPLASTICITY/HUGHES-WINGET``.

.. doxygenclass:: Marmot::Materials::RankineDamagePlasticity
   :allow-dot-graphs:
