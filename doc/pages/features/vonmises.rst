.. _vonMisesModel:

Von Mises model
===============

An implementation of classical J2 plasticity with isotropic hardening.

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

Elastoplastic Constitutive Relations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The infinitesimal elastoplastic constitutive relations, which relate the stress and linearized strain rate tensors by means of the elastic continuum tangent operator, are given as

.. math:: \dot{\sigma}_{ij}
   = C^\mathrm{e}_{ijkl}\dot{\varepsilon}^\mathrm{e}_{kl}
   = C^\mathrm{e}_{ijkl}\left(\dot{\varepsilon}_{kl}-\dot{\varepsilon}^\mathrm{p}_{kl}\right).

Therein, an additive split of the strain rate tensor in elastic and plastic parts is assumed.

The elastic continuum tangent operator is given in terms of Young's modulus :math:`E` and Poisson's ratio :math:`\nu` as

.. math:: C^\mathrm{e}_{ijkl}
   = \frac{E\nu\,\delta_{ij}\delta_{kl}}{\left(1+\nu\right)\left(1-2\nu\right)} + \frac{E\left(\delta_{ik}\delta_{jl}+\delta_{il}\delta_{jk}\right)}{2\left(1+\nu\right)}

Alternatively, the constitutive relations can be rewritten based on a volumetric-deviatoric split of the stress and strain rate tensors

.. math:: \dot{\sigma}_{ij}
   = \dot{p}\delta_{ij}+\dot{s}_{ij} \qquad\text{with}\quad \dot{p}=\frac{\dot{\sigma}_{ii}}{3}

and

.. = \dot{\varepsilon}_{ij}^\mathrm{e} + \dot{\varepsilon}_{ij}^\mathrm{p}

.. math:: \dot{\varepsilon}_{ij}
   = \frac{\dot{\varepsilon}_\mathrm{vol}}{3}\delta_{ij}+\dot{e}_{ij}
   = \left( \dot{\varepsilon}^\mathrm{e}_\mathrm{vol}+\dot{\varepsilon}^\mathrm{p}_\mathrm{vol} \right)\frac{\delta_{ij}}{3} + \dot{e}^\mathrm{e}_{ij}+\dot{e}^\mathrm{p}_{ij} \qquad\text{with}\quad \dot{\varepsilon}_\mathrm{vol}=\dot{\varepsilon}_{ii}

as

.. math:: \dot{\sigma}_{ij}
   = \dot{p}\delta_{ij}+\dot{s}_{ij}
   = K \dot{\varepsilon}^\mathrm{e}_\mathrm{vol} + 2G \dot{e}^\mathrm{e}_{ij}

Yield Function
^^^^^^^^^^^^^^

.. math:: s_{ij} = \sigma_{ij}-\frac{1}{3}\sigma_{kk},

The yield function, which is formulated in terms of the deviatoric part of the stress tensor is given as

.. math:: f\left(\boldsymbol{\sigma},\kappa\right)
   = \sqrt{s_{ij}s_{ij}}-\sqrt{\frac{2}{3}}f_\mathrm{y}\left(\kappa\right),

where :math:`f_\mathrm{y}\left(\kappa\right)` is the yield stress under uniaxial tension, which may depend on the hardening variable :math:`\kappa`.

Flow Rule
^^^^^^^^^

The model encompasses associated plastic flow, i.e., the plastic strain rate

.. math:: \dot{\varepsilon}^\mathrm{p}_{ij}
   = \dot{\lambda}\frac{\partial{}f}{\partial{}\sigma_{ij}}
   = \dot{\lambda}\frac{s_{ij}}{\sqrt{s_{ij}s_{ij}}}
   = \dot{\lambda}\boldsymbol{n}

is proportional to the partial derivatives of the yield function with respect to the stress tensor.
The plastic multiplier :math:`\dot{\lambda}` is a proportionality constant.

.. _vonMisesModel_hardening:

Hardening Law
^^^^^^^^^^^^^

The model includes non-linear isotropic hardening, which is controlled by the evolution of the yield stress as

.. math:: f_\mathrm{y}\left(\kappa\right)
   = f_\mathrm{y}^0 + H_\mathrm{iso,lin}\kappa + \Delta f_\mathrm{y}^{0,\infty}\left(1-\exp\left(-\delta\kappa\right)\right).

Therein, :math:`f_\mathrm{y}^0`, :math:`H_\mathrm{iso,lin}`, :math:`\Delta{}f_\mathrm{y}^{0,\infty}` and :math:`\delta` are material parameters.

.. image:: ./vonmises_hardening.svg

The evolution of the hardening variable is given as

.. math:: \dot{\kappa}
   = \sqrt{\frac{2}{3}\dot{\varepsilon}^\mathrm{p}_{ij}\dot{\varepsilon}^\mathrm{p}_{ij}}
   = \dot{\lambda}\sqrt{\frac{2}{3}}.

Therein, the factor :math:`\sqrt{\frac{2}{3}}` ensures that :math:`\dot{\kappa}` is equal to :math:`\dot{\varepsilon}^\mathrm{p}_{11}` for uniaxial loading in :math:`x_1` direction.

Loading-Unloading Conditions, Consistency Condition
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The plastic multiplier and the yield function obey the Karush–Kuhn–Tucker conditions

.. math:: \dot{\lambda}f=0, \quad \dot{\lambda}\geq0, \quad f\leq0,

as well as the consistency condition

.. math:: \dot{\lambda}\dot{f}=0.

Implementation
--------------

Stress Update Algorithm
^^^^^^^^^^^^^^^^^^^^^^^

Given the strain increment :math:`\Delta\boldsymbol{\varepsilon}`, an elastic trial stress is calculated as

.. math:: \boldsymbol{\sigma}^\mathrm{trial}
   = \boldsymbol{\sigma}^n + \boldsymbol{C}^\mathrm{e}:\Delta\boldsymbol{\varepsilon}.

If the trial stress does not satisfy the yield condition, i.e., :math:`f\left(\boldsymbol{\boldsymbol{\sigma}^\mathrm{trial}},\kappa_n\right)<0`, the step is treated as elastic, which means that the updated stress is equal to the trial stress

.. math:: \boldsymbol{\sigma}^{n+1} = \boldsymbol{\sigma}^\mathrm{trial}

and the hardening variable remains constant

.. math:: \kappa_{n+1} = \kappa_n.

If the trial stress satisfies the yield condition, i.e., :math:`f\left(\boldsymbol{\boldsymbol{\sigma}^\mathrm{trial}},\kappa_n\right)\geq0`, the step is treated as elastoplastic, triggering the radial return mapping algorithm outlined below.

.. .. math:: \boldsymbol{s}^\mathrm{trial}
..    = \boldsymbol{s}^n + 2G\Delta\boldsymbol{e}
..
.. .. math:: \boldsymbol{\sigma}^\mathrm{trial}_\mathrm{m}
..    = \boldsymbol{\sigma}^{n}_\mathrm{m} + K\Delta{}{\varepsilon}_\mathrm{vol}
..
.. .. math:: \boldsymbol{s}_{n+1}
..    = \boldsymbol{s}^n + 2G\Delta\boldsymbol{e}^\mathrm{e} = \boldsymbol{s}^\mathrm{trial}-2G\Delta\boldsymbol{e}^\mathrm{p}
..
..
..
..

The updated deviatoric stress :math:`\boldsymbol{s}^{n+1}` can be written as

.. math:: \boldsymbol{s}^{n+1}
   = \boldsymbol{s}^{\mathrm{trial}} - 2G\Delta\boldsymbol{e}^\mathrm{p}
   = \boldsymbol{s}^{\mathrm{trial}} - 2G\Delta\lambda\boldsymbol{n}^{n+1}.

Taking the norm on both sides results in

.. math:: \|\boldsymbol{s}^{n+1}\|
   = \|\boldsymbol{s}^{\mathrm{trial}}\| - 2G\Delta\lambda.

Inserting the resulting expression for :math:`\|\boldsymbol{s}^{n+1}\|` into the yield function :math:`f^{n+1}`

.. math:: f^{n+1}
   = \|\boldsymbol{s}^{n+1}\| -\sqrt{\frac{2}{3}}f_\mathrm{y}\left(\kappa^{n+1}\right)
   = 0

and substituting :math:`\Delta\lambda` with :math:`\sqrt{\frac{3}{2}}\Delta\kappa` results in the scalar expression

.. math:: g\left(\Delta\kappa\right)
   = \|\boldsymbol{s}^\mathrm{trial}\| -\sqrt{6}G\Delta\kappa -\sqrt{\frac{2}{3}}f_\mathrm{y}\left(\Delta\kappa\right)
   = 0,

which can be solved numerically to obtain :math:`\Delta\kappa`.

If :math:`\Delta\kappa` is known, the internal state variable is updated as

.. math:: \kappa^{n+1} = \kappa^{n} + \Delta\kappa

and the stress update is performed by exploiting :math:`\Delta\lambda=\sqrt{\frac{3}{2}}\Delta\kappa`

.. math:: \boldsymbol{\sigma}^{n+1}
   = \boldsymbol{\sigma}^\mathrm{trial} - 2G\Delta\lambda\frac{\boldsymbol{s}^{\mathrm{trial}}}{\|\boldsymbol{s}^{\mathrm{trial}}\|}

Consistent Algorithmic Tangent Operator
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The consistent tangent operator is given as

.. math:: \frac{\partial{\boldsymbol{s}^{n+1}}}{\partial\boldsymbol{\varepsilon}^{n+1}}
   = \mathbf{C}^\mathrm{e} - 2G\left( \left( 1+\frac{\partial f(\kappa^{n+1})}{\partial\kappa} \right)^{-1} -\frac{2G\Delta\lambda^{n+1}}{\|\boldsymbol{s}^\mathrm{trial}\|} \right) \boldsymbol{n}^{\mathrm{trial}}\otimes\boldsymbol{n}^{\mathrm{trial}} - \frac{4G^2\Delta\lambda}{\|\boldsymbol{s}^\mathrm{trial}\|}\boldsymbol{I}^{\mathrm{dev}}.

.. doxygenclass:: Marmot::Materials::VonMisesModel
   :allow-dot-graphs:

