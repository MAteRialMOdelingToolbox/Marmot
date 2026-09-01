Bergström-Boyce model
=====================

Theory
------

This model implements the classical Bergström-Boyce finite-strain viscoelastic-viscoplastic
formulation for elastomers/polymers as two networks acting in parallel:

- Network A (equilibrium, hyperelastic): sees the total deformation directly.
- Network B ("Maxwell-like", viscous): a spring in series with a viscous flow element.

Both networks use the compressible neo-Hookean potential

.. math::

   \Psi(\boldsymbol{C}) = \frac{\mu}{2}\left(I_1 - 3 - \ln\det\boldsymbol{C}\right)
                        + \frac{\kappa}{8}\left(\ln\det\boldsymbol{C}\right)^2,

evaluated on the total right Cauchy-Green tensor :math:`\boldsymbol{C}` (network A, parameters
:math:`\mu_A,\kappa_A`) and on the elastic right Cauchy-Green tensor :math:`\boldsymbol{C}^{\rm e}`
of network B (parameters :math:`\mu_B,\kappa_B`).

Network B assumes the multiplicative split of the deformation gradient into elastic and
viscous parts,

.. math::

   \boldsymbol{F} = \boldsymbol{F}^{\rm e} \boldsymbol{F}^{\rm v}, \qquad
   \boldsymbol{C}^{\rm e} = (\boldsymbol{F}^{\rm e})^{\mathsf T}\boldsymbol{F}^{\rm e},

with a purely symmetric viscous velocity gradient (:math:`\boldsymbol{W}^{\rm v}=\boldsymbol 0`)
and a chain-stretch / deviatoric-Mandel-stress power-law flow rule:

.. math::

   \boldsymbol{D}^{\rm v} = \dot\gamma \boldsymbol{N}, \qquad
   \dot\gamma = c_1 \left(\lambda^{\rm v}_{\rm chain} - 1\right)^{c_2} \rho^{c_3}, \qquad
   \boldsymbol{N} = \frac{\boldsymbol{\mathcal S}}{\rho},

where :math:`\boldsymbol{\mathcal S}` is the deviatoric Mandel stress of network B's spring,
:math:`\rho = \sqrt{\boldsymbol{\mathcal S}:\boldsymbol{\mathcal S}}`, and
:math:`\lambda^{\rm v}_{\rm chain} = \sqrt{I_1(\boldsymbol{C}^{\rm e})/3}`. Unlike classical
plasticity, there is no yield surface -- the flow rate is a continuous power law that is always
active. Note that :math:`\boldsymbol{N}` is traceless by construction, so network B's flow is
isochoric: :math:`\det\boldsymbol{F}^{\rm e}(t) = \det\boldsymbol{F}(t)` identically, for every
time :math:`t`. This does not imply that network B's hydrostatic stress contribution stays frozen
during a hold, since this potential's :math:`I_1` term still couples the volumetric and isochoric
parts of :math:`\boldsymbol{C}^{\rm e}` -- only :math:`\det\boldsymbol{F}^{\rm e}` itself is exactly
invariant under the isochoric flow.

.. admonition:: Stress update algorithm at a quadrature point for current step :math:`n+1`

   **Trial state:**

   - :math:`\boldsymbol{F}^{\rm e,tr} = \boldsymbol{F}\,(\boldsymbol{F}^{\rm v,old})^{-1}`

   **Return mapping (always active, no yield check):**

   - Solve the nonlinear system for unknowns :math:`\boldsymbol{X}=\{\boldsymbol{F}^{\rm e}_{11},
     \ldots,\boldsymbol{F}^{\rm e}_{33},\Delta\gamma\}^{\mathsf T}` by Newton-Raphson iteration:

     - :math:`\boldsymbol{R}_1`: :math:`\boldsymbol{F}^{\rm e}\exp(\Delta\gamma\boldsymbol{N}) - \boldsymbol{F}^{\rm e,tr} = \boldsymbol{0}`
     - :math:`R_2`: :math:`\Delta\gamma/\Delta t - c_1(\lambda^{\rm v}_{\rm chain}-1)^{c_2}\rho^{c_3} = 0`

     with :math:`\boldsymbol{N}`, :math:`\rho`, :math:`\lambda^{\rm v}_{\rm chain}` evaluated at
     the current iterate (implicit/backward-Euler).

   - Update :math:`\boldsymbol{F}^{\rm v,new} = \left[(\boldsymbol{F}^{\rm e,new})^{-1}\boldsymbol{F}^{\rm e,tr}\right]\boldsymbol{F}^{\rm v,old}`.
   - Total stress: network A's PK2 stress (from :math:`\boldsymbol{C}`) plus network B's PK2
     stress (from :math:`\boldsymbol{C}^{\rm e,new}`) pushed forward through their respective
     deformation gradients, then superposed as Kirchhoff stress.

Reference: J. S. Bergström, M. C. Boyce (1998). Constitutive modeling of the large strain
time-dependent behavior of elastomers. *Journal of the Mechanics and Physics of Solids*
46(5):931-954.

.. note::

   Only the complex-step-differentiated (CSDA) tangent variant (``implementationType = 0``) is
   currently implemented. The fully analytic tangent (``implementationType = 1``) is reserved
   for future work.

Implementation
--------------

.. doxygenclass:: Marmot::Materials::BergstromBoyce
   :allow-dot-graphs:
