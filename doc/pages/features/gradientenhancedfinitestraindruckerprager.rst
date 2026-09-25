Gradient-enhanced finite-strain Drucker-Prager
==============================================

A finite-strain Drucker-Prager damage-plasticity model with implicit-gradient (nonlocal) damage, implementing
``MarmotMaterialGradientEnhancedFiniteStrain``. It can be used with every consumer of that interface, e.g. the
gradient-enhanced finite-strain displacement element and the gradient-enhanced meshfree material point and particle.
It is registered as ``GRADIENTENHANCEDFINITESTRAINDRUCKERPRAGER``.

Theory
------

**Kinematics.** The deformation gradient is split multiplicatively, :math:`\boldsymbol{F} = \boldsymbol{F}^{\rm e}
\boldsymbol{F}^{\rm p}`, with the plastic deformation gradient as state variable. The trial elastic state
:math:`\boldsymbol{F}^{\rm e,trial} = \boldsymbol{F}\,(\boldsymbol{F}^{\rm p}_n)^{-1}` is decomposed spectrally
through its left Cauchy-Green tensor, and the principal elastic logarithmic strains
:math:`\varepsilon^{\rm e}_a = \ln\lambda^{\rm e}_a` are the variables of the return map. The plastic flow is
integrated with the exponential map and the elastic rotation is frozen over the step. For isotropic elasticity and
an isotropic yield function this makes the return map exactly the small-strain one in these variables, and the
plastic deformation gradient is updated as

.. math::

   \boldsymbol{F}^{\rm p}_{n+1} = (\boldsymbol{F}^{\rm e,trial})^{-1}\exp(\Delta\boldsymbol{\varepsilon}^{\rm p})\,\boldsymbol{F},
   \qquad\text{so that}\qquad \boldsymbol{F}(\boldsymbol{F}^{\rm p}_{n+1})^{-1} = \boldsymbol{V}^{\rm e}_{n+1}\boldsymbol{R}^{\rm e}.

**Elasticity** is that of :doc:`compressibleneohooke` (Pence-Gou, variant B) in the elastic stretch,

.. math::

   \Psi = \frac{K}{8}\left(J^2 + J^{-2} - 2\right) + \frac{G}{2}\left(I_1 J^{-2/3} - 3\right),

whose principal Mandel stresses, equal to the principal Kirchhoff stresses by isotropy, are available in closed form,

.. math::

   \Sigma_a = \frac{\partial\Psi}{\partial\varepsilon^{\rm e}_a}
            = \frac{K}{2}\sinh(2\theta) + G\left(e^{2e_a} - \tfrac13\sum_b e^{2e_b}\right),
   \qquad \theta = \sum_a\varepsilon^{\rm e}_a,\quad e_a = \varepsilon^{\rm e}_a - \theta/3 .

**Plasticity.** A Drucker-Prager yield function on the Mandel stress with linear hardening of the cohesion and a
non-associated plastic potential,

.. math::

   f = \sqrt{J_2} + \eta\,p - \xi\,(c_0 + H\alpha), \qquad g = \sqrt{J_2} + \bar\eta\,p,
   \qquad p = \tfrac13\operatorname{tr}\boldsymbol{\Sigma}\ \text{(tension positive)},

where the cone passes through the outer edges of the Mohr-Coulomb pyramid (compressive meridian),

.. math::

   \eta = \frac{6\sin\phi}{\sqrt3\,(3-\sin\phi)},\qquad \xi = \frac{6\cos\phi}{\sqrt3\,(3-\sin\phi)},\qquad
   \bar\eta = \frac{6\sin\psi}{\sqrt3\,(3-\sin\psi)} .

The return to the cone solves
:math:`\varepsilon^{\rm e}_a = \varepsilon^{\rm e,trial}_a - \Delta\lambda\,\partial g/\partial\Sigma_a`,
:math:`f = 0`, :math:`\alpha = \alpha_n + \xi\,\Delta\lambda` with Newton's method and an analytic Jacobian. When
the deviatoric stress would reverse, the state is returned to the apex of the cone instead, with the volumetric
plastic strain increment as unknown and :math:`\alpha = \alpha_n + (\xi/\bar\eta)\,\Delta\varepsilon^{\rm p}_v`
(de Souza Neto, Peric & Owen, *Computational Methods for Plasticity*, Sec. 8.3).

**Damage** is the implicit-gradient damage of the finite-strain damage-plasticity models (GMCDPFiniteStrain and its
gradient-enhanced siblings). The local variable grows with the volumetric plastic logarithmic strain, weighted by a
ductility measure of the compressive part of the plastic flow,

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
:math:`(1-\omega)\,\boldsymbol{\Sigma}:\Delta\boldsymbol{\varepsilon}^{\rm p} + \Psi_{\rm eff}\,\Delta\omega`.
The algorithmic tangents are computed by forward finite differences of the full state update.

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
