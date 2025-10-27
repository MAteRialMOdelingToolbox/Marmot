Linear Viscoelastic Orthotropic Power Law model
===================================

Theory
------

This is an orthotropic implementation of the isotropic model described in :ref:`linearviscoelasticpowerlaw`.

Thus, the normalized initial elastic compliance tensor is replaced by an orthotropic one.

.. math::

   \DelNu = E_1\, \CelInv,

with the Young's modulus :math:`E_1` as the
Young's modulus in the :math:`x_1` material direction.
Constant Poisson's ratios are assumed in all material directions.

Implementation
--------------

.. doxygenclass:: Marmot::Materials::LinearViscoelasticOrthotropicPowerLaw
   :allow-dot-graphs:
