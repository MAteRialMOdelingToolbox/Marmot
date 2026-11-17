.. _compressiblehyperelasticity:

Compressible hyperelasticity
============================

Theory
------

Hyperelastic materials are a class of constitutive models used to describe the elastic behavior of solids.
They postulate a strain-energy density function :math:`\Psi` from which stresses are obtained by differentiation with
respect to appropriate strain measures.
In finite strains we take :math:`\Psi=\Psi(\boldsymbol C)`, where the deformation
gradient is :math:`\boldsymbol F`, the right Cauchy–Green tensor is :math:`\boldsymbol C=\boldsymbol F^{\mathsf T}\boldsymbol F`, and the Jacobian is
:math:`J=\det\boldsymbol F`.

The model uses a decoupled potential, the sum of an isochoric and a volumetric part,

.. math::

   \Psi(\boldsymbol C) = \Psi_{\rm iso}(\bar I_1, \bar I_2) + U(J),
   \qquad \bar I_1 = I_1\,J^{-2/3}, \quad \bar I_2 = I_2\,J^{-4/3},

where :math:`I_1` and :math:`I_2` are the first and second invariant of :math:`\boldsymbol C`. Both parts are selected
by the first two material properties from the shared energy density functions
(``Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions``).

The second Piola-Kirchhoff stress tensor follows from the potential as

.. math::

   \boldsymbol S = 2\,\frac{\partial \Psi}{\partial \boldsymbol C}.

The Kirchhoff and Cauchy stress tensors are obtained by push-forward operations:

.. math::

   \boldsymbol{\tau} = \boldsymbol F\,\boldsymbol S\,\boldsymbol F^{\mathsf T},
   \qquad
   \boldsymbol{\sigma} = J^{-1}\,\boldsymbol{\tau}.

The consistent tangent with respect to the deformation gradient is computed via the chain rule as

.. math::

   \frac{\partial\boldsymbol{\tau}}{\partial\boldsymbol{F}}
   =
   2\,\frac{\partial\boldsymbol{\tau}}{\partial\boldsymbol{S}} : \frac{\partial^{2}\Psi}{\partial\boldsymbol{C}\,\partial\boldsymbol{C}} : \frac{\partial\boldsymbol{C}}{\partial\boldsymbol{F}}
   +
   \frac{\partial}{\partial\boldsymbol{F}}\left(\boldsymbol{F}\,\boldsymbol{S}\,\boldsymbol{F}^{\mathsf{T}}\right).

The first and second derivative of the potential are evaluated in closed form: both parts are scalar functions of the
invariants :math:`(\bar I_1, \bar I_2, J)`, whose scalar derivatives are added before the chain rule to
:math:`\boldsymbol C` is applied once (the single-pass ``DecoupledPotential``).

For infinitesimal strains, the model reduces to linear isotropic elasticity.

Material properties
-------------------

The material properties are provided as a flat vector: the two base selectors, the parameters of the isochoric part,
the parameters of the volumetric part and, optionally, the density.

.. list-table::
   :header-rows: 1
   :align: left

   * - **Index**
     - **Property**
     - **Description**
   * - 0
     - isochoricBase
     - Selector of the isochoric part, see below
   * - 1
     - volumetricBase
     - Selector of the volumetric part, see below
   * - 2 ...
     - isochoric parameters
     - Parameters of the selected isochoric part
   * - ...
     - volumetric parameters
     - Parameters of the selected volumetric part
   * - last
     - :math:`\rho`
     - Density (optional)

.. rubric:: Isochoric part

.. list-table::
   :header-rows: 1
   :align: left

   * - **Base**
     - **selector value**
     - **Parameters (in order)**
     - :math:`\Psi_{\rm iso}`
   * - NeoHooke
     - 0
     - :math:`\mu`
     - :math:`\frac{\mu}{2}(\bar I_1 - 3)`
   * - Yeoh
     - 1
     - :math:`C_{10}`, :math:`C_{20}`, :math:`C_{30}`
     - :math:`\sum_{i=1}^{3} C_{i0}(\bar I_1 - 3)^i`
   * - MooneyRivlin
     - 2
     - :math:`C_{10}`, :math:`C_{01}`
     - :math:`C_{10}(\bar I_1 - 3) + C_{01}(\bar I_2 - 3)`
   * - ArrudaBoyce8Chain
     - 3
     - :math:`\mu`, :math:`\lambda_L` (locking stretch, :math:`> 1`)
     - closed-form 8-chain potential, see ``Isochoric::ArrudaBoyce8ChainPotential``

.. rubric:: Volumetric part

The selector values are the numbers :math:`U_1` to :math:`U_8` of the volumetric strain energy functions in
Doll & Schweizerhof (2000), https://doi.org/10.1115/1.321146. :math:`\kappa` is the bulk modulus.

.. list-table::
   :header-rows: 1
   :align: left

   * - **Base**
     - **selector value**
     - **Parameters (in order)**
     - :math:`U(J)`
   * - SussmanBathe
     - 1
     - :math:`\kappa`
     - :math:`\frac{\kappa}{2}(J-1)^2`
   * - Penalty
     - 2
     - :math:`\kappa`
     - :math:`\frac{\kappa}{2}(\ln J)^2`
   * - SimoTaylor
     - 3
     - :math:`\kappa`
     - :math:`\frac{\kappa}{4}\left[(J-1)^2 + (\ln J)^2\right]`
   * - Ogden
     - 4
     - :math:`\kappa`, :math:`\theta`
     - :math:`\frac{\kappa}{\theta^2}\left(\theta\ln J + J^{-\theta} - 1\right)`
   * - LiuHofstetterMang
     - 5
     - :math:`\kappa`
     - :math:`\kappa\left(J\ln J - J + 1\right)`
   * - DollSchweizerhof
     - 6
     - :math:`\kappa`, :math:`\alpha > 0`, :math:`\beta > 1`
     - :math:`\frac{\kappa}{\alpha+\beta}\left[\frac{J^{\alpha+1}}{\alpha+1} + \frac{J^{-(\beta-1)}}{\beta-1}\right] - \frac{\kappa}{(\alpha+1)(\beta-1)}`
   * - DollSchweizerhofExponential
     - 7
     - :math:`\kappa`
     - :math:`\frac{\kappa}{2}\left(e^{J-1} - \ln J - 1\right)`
   * - DollSchweizerhofLogarithmic
     - 8
     - :math:`\kappa`
     - :math:`\frac{\kappa}{2}(J-1)\ln J`

.. rubric:: Compressible Neo-Hooke model of Pence & Gou

The compressible Neo-Hooke model with the Pence–Gou potential (variant B),

.. math::

   \Psi(\boldsymbol C) = \frac{G}{2}\,\big(I_1\,J^{-2/3}-3\big) + \frac{K}{8}\,\big(J - J^{-1}\big)^2,

with the shear modulus :math:`G` and the bulk modulus :math:`K`, is the pair NeoHooke and DollSchweizerhof with
:math:`\alpha = 1` and :math:`\beta = 3`, i.e. the material properties ``[0, 6, G, K, 1, 3]``. This is the former
material ``COMPRESSIBLENEOHOOKE`` with the properties ``[K, G]``.

Implementation
--------------

.. doxygenclass:: Marmot::Materials::CompressibleHyperelasticity
   :allow-dot-graphs:
