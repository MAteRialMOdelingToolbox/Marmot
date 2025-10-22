Compressible Neo Hooke model
============================

Theory
------

Hyperelastic materials are a class of constitutive models used to describe the elastic behavior of solids.
They postulate a strain-energy density function :math:`\Psi` from which stresses are obtained by differentiation with
respect to appropriate strain measures.
In finite strains we take :math:`\Psi=\Psi(\boldsymbol C)`, where the deformation
gradient is :math:`\boldsymbol F`, the right Cauchy–Green tensor is :math:`\boldsymbol C=\boldsymbol F^{\mathsf T}\boldsymbol F`, and the Jacobian is
:math:`J=\det\boldsymbol F`. The first invariant is :math:`I_1=\operatorname{tr}\boldsymbol C`.

The model implemented here is a compressible Neo-Hookean material with the Pence–Gou potential (variant B).
An isochoric–volumetric split is defined by

.. math::

   \Psi(\boldsymbol C)
   =
   \Psi_{\rm d}(I_1,J) + \Psi_{\rm h}(J)
   =
   \frac{G}{2}\,\big(I_1\,J^{-2/3}-3\big)
   +
   \frac{K}{8}\,\big(J - J^{-1}\big)^2,

where :math:`G` is the shear modulus and :math:`K` is the bulk modulus.

The second Piola-Kirchhoff stress tensor follows from the potential as

.. math::

   \boldsymbol S = 2\,\frac{\partial \Psi}{\partial \boldsymbol C}.


The Kirchhoff and Cauchy stress tensors are obtained by push-forward operations:

.. math::

   \boldsymbol{\tau} = \boldsymbol F\,\boldsymbol S\,\boldsymbol F^{\mathsf T},

.. math::

   \boldsymbol{\sigma} = J^{-1}\,\boldsymbol{\tau}.


Using the chain rule, the first derivative of :math:`\Psi` with respect to :math:`\boldsymbol C` is

.. math::

   \frac{\partial \Psi}{\partial \boldsymbol C}
   =
   \frac{\partial \Psi}{\partial J}\,\frac{\partial J}{\partial \boldsymbol C}
   +
   \frac{\partial \Psi}{\partial I_1}\,\frac{\partial I_1}{\partial \boldsymbol C}
   =
   \Big[
     \frac{K}{4}\,(J - J^{-1})\!\left(1+\frac{1}{J^{2}}\right)
     - \frac{G}{3}\,I_1\,J^{-5/3}
   \Big]\,
   \frac{J}{2}\,\boldsymbol C^{-{\mathsf T}}
   +
   \frac{G}{2}\,J^{-2/3}\,\boldsymbol I.

The resulting second Piola–Kirchhoff stress used in the implementation is

.. math::

   \boldsymbol S
   =
   G\,J^{-2/3}\Big(\boldsymbol I - \tfrac{1}{3}\,I_1\,\boldsymbol C^{-1}\Big)
   +
   \frac{K}{4}\,\big(J^{2}-J^{-2}\big)\,\boldsymbol C^{-1},

where :math:`I_1 = \operatorname{tr}\boldsymbol C`.

The resulting Kirchhoff stress, after the push-forward operation :math:`\boldsymbol{\tau} = \boldsymbol F\,\boldsymbol S\,\boldsymbol F^{\mathsf T}`, is:

.. math::

   \boldsymbol{\tau}
   =
   G\,J^{-2/3}\Big(\boldsymbol{b} - \tfrac{1}{3}\,\operatorname{tr}(\boldsymbol{b})\,\boldsymbol{I}\Big)
   +
   \frac{K}{4}\,\big(J^{2}-J^{-2}\big)\,\boldsymbol{I},

where :math:`\boldsymbol{b} = \boldsymbol F\,\boldsymbol F^{\mathsf T}` is the left Cauchy-Green tensor.

The consistent tangent with respect to the deformation gradient :math:`\frac{\partial\boldsymbol{\tau}(\boldsymbol{S}(\boldsymbol{C}(\boldsymbol{F})), \boldsymbol{F})}{\partial\boldsymbol{F}}` is computed via the chain rule as:

.. math::

   \frac{\partial\boldsymbol{\tau}}{\partial\boldsymbol{F}}
   =
   2\,\frac{\partial\boldsymbol{\tau}}{\partial\boldsymbol{S}} : \frac{\partial^{2}\Psi}{\partial\boldsymbol{C}\,\partial\boldsymbol{C}} : \frac{\partial\boldsymbol{C}}{\partial\boldsymbol{F}}
   +
   \frac{\partial}{\partial\boldsymbol{F}}\left(\boldsymbol{F}\,\boldsymbol{S}\,\boldsymbol{F}^{\mathsf{T}}\right).

The individual derivatives are:

.. math::

   \frac{\partial\boldsymbol{\tau}}{\partial\boldsymbol{S}} = \boldsymbol{F}\otimes\boldsymbol{F},

.. math::

   \frac{\partial^{2}\Psi}{\partial \boldsymbol C\,\partial \boldsymbol C}
   =
   \frac{\partial^{2}\Psi}{\partial J^{2}}\,
   \frac{\partial J}{\partial \boldsymbol C}\otimes\frac{\partial J}{\partial \boldsymbol C}
   +
   \frac{\partial \Psi}{\partial J}\,
   \frac{\partial^{2} J}{\partial \boldsymbol C\,\partial \boldsymbol C}
   +
   \frac{\partial^{2}\Psi}{\partial J\,\partial I_1}
   \left(
     \frac{\partial J}{\partial \boldsymbol C}\otimes\frac{\partial I_1}{\partial \boldsymbol C}
     +
     \frac{\partial I_1}{\partial \boldsymbol C}\otimes\frac{\partial J}{\partial \boldsymbol C}
   \right),

.. math::

   \frac{\partial\boldsymbol{C}}{\partial\boldsymbol{F}} = \boldsymbol{F}^{\mathsf{T}}\otimes\boldsymbol{I} + \boldsymbol{I}\otimes\boldsymbol{F}^{\mathsf{T}},

.. math::

   \frac{\partial}{\partial\boldsymbol{F}}\left(\boldsymbol{F}\,\boldsymbol{S}\,\boldsymbol{F}^{\mathsf{T}}\right) = \boldsymbol{I}\otimes(\boldsymbol{S}\boldsymbol{F})^{\mathsf{T}} + \boldsymbol{S}\boldsymbol{F}\otimes\boldsymbol{I}.

For infinitesimal strains, this model reduces to linear isotropic elasticity.


Implementation
--------------

.. doxygenclass:: Marmot::Materials::CompressibleNeoHooke
   :allow-dot-graphs:
