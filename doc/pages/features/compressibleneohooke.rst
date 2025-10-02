Compressible Neo Hooke model
============================

Theory
------

Hyperelastic materials are a class of constitutive models used to describe the elastic behavior of solids.
They postulate a strain-energy density function :math:`\Psi` from which stresses are obtained by differentiation with
respect to appropriate strain measures.
In finite strains we take :math:`\Psi=\Psi(\mathbf C)`, where the deformation
gradient is :math:`\mathbf F`, the right Cauchy–Green tensor is :math:`\mathbf C=\mathbf F^{\mathsf T}\mathbf F`, and the Jacobian is
:math:`J=\det\mathbf F`. The first invariant is :math:`I_1=\operatorname{tr}\mathbf C`.

The model implemented here is a compressible Neo-Hookean material with the Pence–Gou potential (variant B).
An isochoric–volumetric split is defined by

.. math::

   \Psi(\mathbf C)
   =
   \Psi_d(I_1,J) + \Psi_h(J)
   =
   \frac{G}{2}\,\big(I_1\,J^{-2/3}-3\big)
   +
   \frac{K}{8}\,\big(J - J^{-1}\big)^2,

where :math:`G` is the shear modulus and :math:`K` is the bulk modulus.

The second Piola-Kirchhoff stress tensor follows from the potential as

.. math::

   \mathbf S = 2\,\frac{\partial \Psi}{\partial \mathbf C}.


The Kirchhoff and Cauchy stress tensors are obtained by push-forward operations:

.. math::

   \boldsymbol{\tau} = \mathbf F\,\mathbf S\,\mathbf F^{\mathsf T},

.. math::

   \boldsymbol{\sigma} = J^{-1}\,\boldsymbol{\tau}.


By using the chain rule, the first derivative of :math:`\Psi` with respect to :math:`\mathbf C` is

.. math::

   \frac{\partial \Psi}{\partial \mathbf C}
   =
   \frac{\partial \Psi}{\partial J}\,\frac{\partial J}{\partial \mathbf C}
   +
   \frac{\partial \Psi}{\partial I_1}\,\frac{\partial I_1}{\partial \mathbf C}
   =
   \Big[
     \frac{K}{4}\,(J - J^{-1})\!\left(1+\frac{1}{J^{2}}\right)
     - \frac{G}{3}\,I_1\,J^{-5/3}
   \Big]\,
   \frac{J}{2}\,\mathbf C^{-{\mathsf T}}
   +
   \frac{G}{2}\,J^{-2/3}\,\mathbf I.

The resulting second Piola–Kirchhoff stress used in the implementation is

.. math::

   \mathbf S
   =
   G\,J^{-2/3}\Big(\mathbf I - \tfrac{1}{3}\,I_1\,\mathbf C^{-1}\Big)
   +
   \frac{K}{4}\,\big(J^{2}-J^{-2}\big)\,\mathbf C^{-1}.

The second derivative of :math:`\Psi` with respect to :math:`\mathbf C` is

.. math::

   \frac{\partial^{2}\Psi}{\partial \mathbf C\,\partial \mathbf C}
   =
   \frac{\partial^{2}\Psi}{\partial J^{2}}\,
   \frac{\partial J}{\partial \mathbf C}\otimes\frac{\partial J}{\partial \mathbf C}
   +
   \frac{\partial \Psi}{\partial J}\,
   \frac{\partial^{2} J}{\partial \mathbf C\,\partial \mathbf C}
   +
   \frac{\partial^{2}\Psi}{\partial J\,\partial I_1}
   \left(
     \frac{\partial J}{\partial \mathbf C}\otimes\frac{\partial I_1}{\partial \mathbf C}
     +
     \frac{\partial I_1}{\partial \mathbf C}\otimes\frac{\partial J}{\partial \mathbf C}
   \right).

For infinitesimal strains this model reduces to linear isotropic elasticity.


Implementation
--------------

.. doxygenclass:: Marmot::Materials::CompressibleNeoHooke
   :allow-dot-graphs:
