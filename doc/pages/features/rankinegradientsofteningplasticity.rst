.. _rankineGradientSofteningPlasticity:

Gradient-enhanced Rankine softening plasticity model
====================================================

Implicit-gradient, over-nonlocal version of the :ref:`rankineSofteningPlasticity`: multi-surface Rankine plasticity
whose yield stress depends on a combination of the local and the nonlocal softening variable (Vermeer and Brinkgreve,
1994), for the general gradient-enhanced displacement element.

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
     - :math:`f_\mathrm{y}^{0}`
     - Initial uniaxial tensile strength
   * - 3
     - :math:`H_\mathrm{iso,lin}`
     - First hardening/softening parameter
   * - 4
     - :math:`\Delta f_\mathrm{y}^{0,\infty}`
     - Second hardening/softening parameter
   * - 5
     - :math:`\delta`
     - Third hardening/softening parameter
   * - 6
     - :math:`l`
     - Nonlocal length (nonnegative)
   * - 7
     - :math:`m`
     - Over-nonlocal weight (0: local, 1: nonlocal, > 1: over-nonlocal)
   * - 8
     - :math:`\rho`
     - Density (optional)
   * - 9
     - :math:`\eta`
     - Nonlocal viscosity (optional, default 0; explicit dynamics only)

.. list-table::
   :header-rows: 1
   :align: left

   * - **State Variable**
     - **Name**
     - **Description**
   * - :math:`\kappa`
     - ``kappa``
     - Local softening variable
   * -
     - ``acYield``
     - Yield flag: 1 for a plastic, 0 for an elastic increment
   * - :math:`\bar\kappa`
     - ``kappaNonlocal``
     - Nonlocal softening variable
   * - :math:`\hat\kappa`
     - ``kappaWeighted``
     - Over-nonlocal combination

Theory
------

The yield functions, the associated flow rule and the local softening variable :math:`\kappa`, with
:math:`\Delta\kappa = \sum_A \Delta\lambda_A`, are those of the :ref:`rankineSofteningPlasticity`. The local softening
variable is the source of the Helmholtz equation solved by the element,

.. math:: \bar\kappa - l^2\,\Delta\bar\kappa = \kappa, \qquad
   \nabla\bar\kappa\cdot\boldsymbol{n} = 0 \text{ on the boundary}.

The yield stress follows from the over-nonlocal combination

.. math:: \hat\kappa = m\,\bar\kappa + (1-m)\,\kappa, \qquad
   f_\mathrm{y}\left(\hat\kappa\right)
   = \max\left(f_\mathrm{y}^0 + H_\mathrm{iso,lin}\hat\kappa + \Delta f_\mathrm{y}^{0,\infty}\left(1-\exp\left(-\delta\hat\kappa\right)\right),\ 0\right).

For :math:`m = 0` the model is the :ref:`rankineSofteningPlasticity`. For a softening law and :math:`m > 1`, the yield
stress increases with :math:`\kappa` at a fixed nonlocal variable: the material point hardens locally, and softening
enters by the nonlocal variable only, which prevents the plastic strain from localizing in a set of measure zero.
The yield stress is floored at zero in terms of :math:`\hat\kappa`, and the dissipated plastic work is
:math:`\boldsymbol{\sigma}:\Delta\boldsymbol{\varepsilon}^\mathrm{p} = f_\mathrm{y}\left(\hat\kappa\right)\Delta\kappa`.

Implementation
--------------

Stress Update Algorithm
^^^^^^^^^^^^^^^^^^^^^^^

The backward Euler return mapping in principal space is that of the :ref:`rankineSofteningPlasticity`, with the
residuals

.. math:: r_A = \sigma_A^\mathrm{trial} - \sum_{B\in\mathcal{A}} C^\mathrm{e}_{AB}\Delta\lambda_B
   - f_\mathrm{y}\left(m\,\bar\kappa + (1-m)\left(\kappa^n + \sum_B\Delta\lambda_B\right)\right) = 0

for all surfaces in the active set :math:`\mathcal{A}`, and the Jacobian

.. math:: J_{AB} = C^\mathrm{e}_{AB} + (1-m)\,\frac{\partial f_\mathrm{y}}{\partial\hat\kappa}, \qquad A,B\in\mathcal{A}.

Since the yield stress depends on the nonlocal variable, an increment without a strain increment may be plastic.

Consistent Algorithmic Tangent Operator
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The tangent :math:`\partial\boldsymbol{\sigma}/\partial\Delta\boldsymbol{\varepsilon}` follows as for the
:ref:`rankineSofteningPlasticity`, with the hardening modulus :math:`(1-m)\,\partial f_\mathrm{y}/\partial\hat\kappa`.
The nonlocal variable shifts the yield stress of the active surfaces, which gives the coupling terms

.. math:: \frac{\partial\Delta\lambda_A}{\partial\bar\kappa}
   = -m\,\frac{\partial f_\mathrm{y}}{\partial\hat\kappa} \sum_{B\in\mathcal{A}} J^{-1}_{AB}, \qquad
   \frac{\partial\boldsymbol{\sigma}}{\partial\bar\kappa}
   = -\sum_{A}\sum_{B\in\mathcal{A}} C^\mathrm{e}_{AB}\,\frac{\partial\Delta\lambda_B}{\partial\bar\kappa}\,
   \boldsymbol{n}_A\otimes\boldsymbol{n}_A, \qquad
   \frac{\partial\kappa}{\partial\bar\kappa} = \sum_{A\in\mathcal{A}} \frac{\partial\Delta\lambda_A}{\partial\bar\kappa}.

Since the flow directions have unit trace, :math:`\Delta\kappa = \mathrm{tr}\left(\Delta\boldsymbol{\varepsilon} -
\mathbb{C}^{-1}:\Delta\boldsymbol{\sigma}\right)`, and

.. math:: \frac{\partial\kappa}{\partial\Delta\boldsymbol{\varepsilon}}
   = \boldsymbol{1} - \left(\frac{\partial\boldsymbol{\sigma}}{\partial\Delta\boldsymbol{\varepsilon}}\right)^T :
   \mathbb{C}^{-1} : \boldsymbol{1}.

The model is registered as ``RANKINEGRADIENTSOFTENINGPLASTICITY`` in the gradient-enhanced hypoelastic factory with one
nonlocal field.

.. doxygenclass:: Marmot::Materials::RankineGradientSofteningPlasticity
   :allow-dot-graphs:
