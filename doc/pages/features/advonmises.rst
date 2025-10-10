Von Mises model (Automatic Differentiation)
===========================================

An automatic differentiation implementation of classical J2 plasticity with isotropic hardening.

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
     - Initial yield stress
   * - 3
     - :math:`H_\mathrm{iso,lin}`
     - First hardening parameter [#f1]_
   * - 4
     - :math:`\Delta f_\mathrm{y}^{0,\infty}`
     - Second hardening parameter [#f1]_
   * - 5
     - :math:`\delta`
     - Third hardening parameter [#f1]_

.. list-table::
   :header-rows: 1
   :align: left

   * - **State Variable**
     - **Name**
     - **Description**
   * - :math:`\kappa`
     - ``kappa``
     - Hardening variable [#f1]_

.. [#f1] see :ref:`vonMisesModel_hardening`.

Theory
------

See :ref:`vonMisesModel`.

Implementation
--------------

Implementation follows the :ref:`vonMisesModel`, however, automatic differentiation is used to compute the consistent algorithmic tangent operator.

.. doxygenclass:: Marmot::Materials::ADVonMises
   :allow-dot-graphs:
