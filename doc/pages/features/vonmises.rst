Von Mises model
===============

Theory
------

Elastoplastic Constitutive Relations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The elastoplastic constitutive relations for the infinitesimal model, which relate the stress and strain rate tensors by means of the elastic continuum tangent operator, are given as

.. math:: \dot{\sigma}_{ij}
   = \mathsf{C}^\mathrm{e}_{ijkl}\dot{\varepsilon}^\mathrm{e}_{kl}
   = \mathsf{C}^\mathrm{e}_{ijkl}\left(\dot{\varepsilon}_{kl}-\dot{\varepsilon}^\mathrm{p}_{kl}\right).

Therein, an additive split of the strain rate tensor in elastic and plastic parts is assumed.

The elastic continuum tangent operator is given in terms of Young's modulus :math:`E` and Poisson's ratio :math:`\nu` as

.. math:: \mathsf{C}^\mathrm{e}_{ijkl}
   = \frac{E\nu\,\delta_{ij}\delta_{kl}}{\left(1+\nu\right)\left(1-2\nu\right)} + \frac{E\left(\delta_{ik}\delta_{jl}+\delta_{il}\delta_{jk}\right)}{2\left(1+\nu\right)}

Yield Function
^^^^^^^^^^^^^^

The yield function, which is formulated in terms of the deviatoric stress tensor

.. math:: s_{ij} = \sigma_{ij}-\frac{1}{3}\sigma_{kk},

is given as

.. math:: f\left(\boldsymbol{\sigma},\kappa\right)
   = \sqrt{s_{ij}s_{ij}}-\sqrt{\frac{2}{3}}f_\mathrm{y}\left(\kappa\right),

where :math:`f_\mathrm{y}\left(\kappa\right)` is the yield stress under uniaxial tension, which depends on the hardening variable :math:`\kappa`.

Flow Rule
^^^^^^^^^

The model encompasses associated plastic flow, i.e., the plastic strain rate

.. math:: \dot{\varepsilon}^\mathrm{p}_{ij}=\dot{\lambda}\frac{\partial{}f}{\partial{}\sigma_{ij}}

is proportional to the partial derivatives of the yield function with respect to the stress tensor.
The plastic multiplier :math:`\dot{\lambda}` is a proportionality constant.

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


Implementation
--------------

.. doxygenclass:: Marmot::Materials::VonMisesModel
   :allow-dot-graphs:
