Compressible Neo Hooke model (Automatic Differentiation)
=========================================================

An automatic differentiation implementation of the compressible Neo-Hookean hyperelastic material model.
This material is intended for demonstration and verification purposes, illustrating how to use the capabilities
of the AD (Automatic Differentiation) material framework.

Theory
------

See :ref:`compressibleneohooke` for the underlying theory.

Implementation
--------------

The implementation follows :ref:`compressibleneohooke`; however,
automatic differentiation is used to compute the consistent algorithmic tangent operator.

.. list-table::
   :header-rows: 1
   :align: left

   * - **Index**
     - **Model Parameter**
     - **Description**
   * - 0
     - :math:`K`
     - Bulk modulus
   * - 1
     - :math:`G`
     - Shear modulus

.. list-table::
   :header-rows: 1
   :align: left

   * - **State Variable**
     - **Name**
     - **Description**

   * - ---
     - ---
     - No state variables required.

.. doxygenclass:: Marmot::Materials::ADCompressibleNeoHooke
   :allow-dot-graphs:
