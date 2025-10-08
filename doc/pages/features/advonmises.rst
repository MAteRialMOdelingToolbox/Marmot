Von Mises model (Automatic Differentiation)
===========================================

An automatic differentiation implementation of classical J2 plasticity with isotropic hardening.

Theory
------

See :ref:`vonMisesModel`.

Implementation
--------------

Implementation follows the :ref:`vonMisesModel`, however, automatic differentiation is used to compute the consistent algorithmic tangent operator.

.. doxygenclass:: Marmot::Materials::ADVonMisesModel
   :allow-dot-graphs:
