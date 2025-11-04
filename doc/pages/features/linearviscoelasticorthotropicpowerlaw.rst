Linear Viscoelastic Orthotropic Power Law model
===================================

Theory
------

This is an orthotropic implementation of the isotropic model described in :ref:`linearviscoelasticpowerlaw`.
By contrast to the isotropic model, an additional scaling factor is used for scaling the stiffness tensor.
This is helpful because the experimentally determined Young's moduli and shear moduli
are in general not equal to the initial elastic moduli of the viscoelastic model.
By using this scaling factor, the experimentally determined Young's moduli and shear moduli
can be used as input parameters for the viscoelastic model and the initial elastic moduli
are automatically adjusted.
In the following, the scaling factor is denoted as :math:`s` which is usually larger than one.

For the formulation of the viscoelastic model in :ref:`linearviscoelasticpowerlaw`,
the initial elastic stiffness tensor :math:`\Cel` is replaced by

.. math::

   \Cel_0 = s\, \Cel( E_1, E_2, E_3, G_{12}, G_{23}, G_{31}, \nu_{12}, \nu_{23}, \nu_{31} ),
where the Young's moduli :math:`E_1`, :math:`E_2`, and :math:`E_3`
and the shear moduli :math:`G_{12}`, :math:`G_{23}`, and :math:`G_{31}`
are given in the material directions :math:`x_1`, :math:`x_2`, and :math:`x_3`.
Additionally, the Poisson's ratios :math:`\nu_{12}`, :math:`\nu_{23}`, and :math:`\nu_{31}` are used.

.. note::
   The scaling factor :math:`s` is applied to all initial elastic moduli.
   This means that the ratios between the Young's moduli and shear moduli remain unchanged.
   Only the absolute values of the initial elastic moduli are scaled.

.. warning::
   The defnition of Poisson's ratio in Marmot differs from the standard definition.
   In Marmot, Poisson's ratio :math:`\nu_{ij}` is defined as the ratio of the lateral strain in direction :math:`x_j`
   to the axial strain in direction :math:`x_i` when a uniaxial stress is applied in direction :math:`x_i`.
   See namespace :ref:`Marmot::ContinuumMechanics::Elasticity` for more details.

Accordingly, the normalized initial elastic compliance tensor now reads

.. math::

   \DelNu = s\, E_1\, \CelInv_0,

with the Young's modulus :math:`E_1` as the
Young's modulus in the :math:`x_1` material direction.
Constant Poisson's ratios are assumed in all material directions.

Implementation
--------------

.. doxygenclass:: Marmot::Materials::LinearViscoelasticOrthotropicPowerLaw
   :allow-dot-graphs:
