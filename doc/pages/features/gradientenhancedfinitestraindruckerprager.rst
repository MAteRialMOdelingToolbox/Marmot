Gradient-enhanced finite-strain Drucker-Prager
==============================================

A finite-strain Drucker-Prager damage-plasticity model with implicit-gradient (nonlocal) damage, implementing
``MarmotMaterialGradientEnhancedFiniteStrain``. It can be used with every consumer of that interface, e.g. the
gradient-enhanced finite-strain displacement element and the gradient-enhanced meshfree material point and particle.
It is registered as ``GRADIENTENHANCEDFINITESTRAINDRUCKERPRAGER``.

Theory
------

**Kinematics.** The deformation gradient is split multiplicatively, :math:`\boldsymbol{F} = \boldsymbol{F}^{\rm e}
\boldsymbol{F}^{\rm p}`, with the plastic deformation gradient as state variable. The formulation follows
:doc:`finitestrainj2plasticity`: the yield function is evaluated on the Mandel stress
:math:`\boldsymbol{M} = \boldsymbol{C}^{\rm e}\boldsymbol{S}`, and the plastic flow is integrated with the exponential
map,

.. math::

   \boldsymbol{F}^{\rm e,trial} = \boldsymbol{F}\,(\boldsymbol{F}^{\rm p}_n)^{-1}
   = \boldsymbol{F}^{\rm e}\exp\!\left(\Delta\lambda\,\frac{\partial g}{\partial\boldsymbol{M}}\right) .

**Elasticity** is that of :doc:`compressibleneohooke` (Pence-Gou, variant B) in the elastic stretch,

.. math::

   \Psi = \frac{K}{8}\left(J^2 + J^{-2} - 2\right) + \frac{G}{2}\left(I_1 J^{-2/3} - 3\right),

in :math:`J = \det\boldsymbol{F}^{\rm e}` and :math:`I_1 = \operatorname{tr}\boldsymbol{C}^{\rm e}`.

**Plasticity.** A Drucker-Prager yield function on the Mandel stress with linear hardening of the cohesion and a
non-associated plastic potential,

.. math::

   f = \sqrt{J_2} + \eta\,p - \xi\,(c_0 + H\alpha), \qquad g = \sqrt{J_2} + \bar\eta\,p,
   \qquad p = \tfrac13\operatorname{tr}\boldsymbol{M}\ \text{(tension positive)},

where the cone passes through the outer edges of the Mohr-Coulomb pyramid (compressive meridian),

.. math::

   \eta = \frac{6\sin\phi}{\sqrt3\,(3-\sin\phi)},\qquad \xi = \frac{6\cos\phi}{\sqrt3\,(3-\sin\phi)},\qquad
   \bar\eta = \frac{6\sin\psi}{\sqrt3\,(3-\sin\psi)} .

The return to the cone solves the flow rule above, the hardening law :math:`\alpha = \alpha_n + \xi\,\Delta\lambda`
and the consistency condition :math:`f = 0` for :math:`\{\boldsymbol{F}^{\rm e}, \alpha, \Delta\lambda\}` with Newton's
method. Where no solution on the cone exists (a trial state beyond the apex), the state returns to the apex: by
isotropy, :math:`\boldsymbol{F}^{\rm p}` is determined up to a rotation only, so that
:math:`\boldsymbol{F}^{\rm e} = J_{\rm e}^{1/3}\boldsymbol{I}` with the unknowns :math:`\{\ln J_{\rm e}, \alpha\}`,
:math:`\eta\,p = \xi\,(c_0 + H\alpha)` and :math:`\alpha = \alpha_n + (\xi/\bar\eta)\,\Delta\varepsilon^{\rm p}_v`
(de Souza Neto, Peric & Owen, *Computational Methods for Plasticity*, Sec. 8.3).

The Jacobians of both return mappings are computed by the complex-step method, and the algorithmic tangents follow
from the same Jacobians by the implicit function theorem.

**Damage** is an implicit-gradient damage driven by the plastic flow. The local variable grows with the volumetric
plastic logarithmic strain, weighted by the ductility measure of the concrete damage-plasticity model CDPM2 (Grassl,
Xenos, Nyström, Rempling & Gylltoft, *Int. J. Solids Struct.* 50, 2013), which reduces the damage under confined,
compression-dominated flow,

.. math::

   \Delta\alpha_{\rm local} = \frac{\Delta\varepsilon^{\rm p}_v}{x_s(R_s)},\qquad
   R_s = \frac{\sum_a\langle-\Delta\varepsilon^{\rm p}_a\rangle}{\Delta\varepsilon^{\rm p}_v},\qquad
   x_s = \begin{cases} 1 + A_s R_s^2 & R_s < 1 \\ 1 + A_s(4\sqrt{R_s} - 3) & R_s \ge 1 \end{cases},

and is the source :math:`L` of the nonlocal balance :math:`\bar N - l^2\nabla^2\bar N = L`. The damage follows from
the over-nonlocal weighting of the local and the nonlocal measure, with a history maximum that makes it
irreversible,

.. math::

   \kappa = \max_t\left(m\,\bar N + (1-m)\,\alpha_{\rm local}\right),\qquad
   \omega = \min\left(1 - e^{-\kappa/\varepsilon_f},\ \omega_{\max}\right),\qquad
   \boldsymbol{\tau} = (1-\omega)\,\boldsymbol{\tau}_{\rm eff}.

The dissipation is cumulative: the incoming value is incremented by
:math:`(1-\omega)\,\boldsymbol{M}:\Delta\boldsymbol{\varepsilon}^{\rm p} + \Psi_{\rm eff}\,\Delta\omega`.

Material parameters
-------------------

.. list-table::
   :header-rows: 1

   * - Index
     - Symbol
     - Description
   * - 0
     - :math:`K`
     - bulk modulus
   * - 1
     - :math:`G`
     - shear modulus
   * - 2
     - :math:`c_0`
     - cohesion
   * - 3
     - :math:`\phi`
     - friction angle [deg]
   * - 4
     - :math:`\psi`
     - dilatancy angle [deg], :math:`0 \le \psi \le \phi`
   * - 5
     - :math:`H`
     - linear hardening modulus of the cohesion
   * - 6
     - :math:`A_s`
     - ductility parameter of the damage
   * - 7
     - :math:`\varepsilon_f`
     - softening modulus of the damage
   * - 8
     - :math:`\omega_{\max}`
     - maximum damage, :math:`0 \le \omega_{\max} < 1`
   * - 9
     - :math:`l`
     - nonlocal radius
   * - 10
     - :math:`m`
     - weighting of the nonlocal measure (:math:`m > 1`: over-nonlocal)
   * - 11
     - :math:`\rho`
     - density (optional)

State variables: ``Fp`` (9), ``alphaP``, ``alphaD`` (the local damage variable), ``kappa``, ``omega``. The plastic
deformation gradient must be initialized to the identity by ``initializeYourself`` — in EdelweissFE with
``>>initializematerial`` in the step.

Implementation
--------------

.. doxygenclass:: Marmot::Materials::GradientEnhancedFiniteStrainDruckerPrager
   :allow-dot-graphs:
