Bergström-Boyce model
=====================

Theory
------

This model implements the classical Bergström-Boyce finite-strain viscoelastic-viscoplastic
formulation for elastomers/polymers as two networks acting in parallel:

- Network A (equilibrium, hyperelastic): sees the total deformation directly.
- Network B ("Maxwell-like", viscous): a spring in series with a viscous flow element.

Both networks share the same selectable hyperelastic base potential
(``hyperelasticBase``), evaluated on the total right Cauchy-Green tensor :math:`\boldsymbol{C}`
for network A (coefficients :math:`A_1,A_2,A_3`, bulk modulus :math:`\kappa_A`) and on the
elastic right Cauchy-Green tensor :math:`\boldsymbol{C}^{\rm e}` for network B's spring
(coefficients :math:`B_1,B_2,B_3`, bulk modulus :math:`\kappa_B`). NeoHooke, Yeoh and
Mooney-Rivlin share the same volumetric term and are constructed to be stress-free at
:math:`\boldsymbol{C}=\boldsymbol{I}`:

.. math::

   \Psi_{\rm NeoHooke}(\boldsymbol{C}) &= \frac{\mu}{2}\left(I_1 - 3 - \ln\det\boldsymbol{C}\right)
                        + \frac{\kappa}{8}\left(\ln\det\boldsymbol{C}\right)^2
   \qquad (\text{coefficient 1} = \mu,\ 2,3\text{ unused}) \\[4pt]
   \Psi_{\rm Yeoh}(\boldsymbol{C}) &= C_{10}(I_1-3) + C_{20}(I_1-3)^2 + C_{30}(I_1-3)^3
                        - C_{10}\ln\det\boldsymbol{C} + \frac{\kappa}{8}\left(\ln\det\boldsymbol{C}\right)^2
   \qquad (\text{coefficients} = C_{10},C_{20},C_{30}) \\[4pt]
   \Psi_{\rm MooneyRivlin}(\boldsymbol{C}) &= C_{10}(I_1-3) + C_{01}(I_2-3)
                        - (C_{10}+2C_{01})\ln\det\boldsymbol{C} + \frac{\kappa}{8}\left(\ln\det\boldsymbol{C}\right)^2
   \qquad (\text{coefficients} = C_{10},C_{01},\ 3\text{ unused})

with :math:`I_2=\tfrac12\left(I_1^2-{\rm tr}(\boldsymbol C^2)\right)`. Yeoh reduces exactly to
NeoHooke when :math:`C_{20}=C_{30}=0,\ C_{10}=\mu/2`, and Mooney-Rivlin reduces exactly to
NeoHooke when :math:`C_{01}=0,\ C_{10}=\mu/2` (both verified as exact regression checks in the
unit test).

A fourth base, Arruda-Boyce (8-chain), is built on the isochoric invariant
:math:`\bar I_1 = I_1\det\boldsymbol{C}^{-1/3}` rather than the raw :math:`I_1` the other three
bases use, via the closed-form Cohen (1991) Padé approximation to the inverse Langevin function,
integrated in :math:`\bar I_1`:

.. math::

   \Psi_{\rm ArrudaBoyce}(\boldsymbol{C}) = \underbrace{\frac{\mu}{6}\left(\bar I_1-3\right)
      - \mu\lambda_L^2\ln\left(\frac{1-\bar I_1/(3\lambda_L^2)}{1-1/\lambda_L^2}\right)}_{\Psi_{\rm iso}(\bar I_1)}
      + \frac{\kappa}{8}\left(\ln\det\boldsymbol{C}\right)^2
   \qquad (\text{coefficients 1,2} = \mu,\lambda_L,\ 3\text{ unused})

where :math:`\mu` is the shear-modulus-like parameter and :math:`\lambda_L` is the locking
stretch. Because :math:`\bar I_1` is invariant under :math:`\boldsymbol C\to\lambda\boldsymbol C`
for any scalar :math:`\lambda`, :math:`\Psi_{\rm iso}` is automatically stress-free at
:math:`\boldsymbol{C}=\boldsymbol{I}` -- unlike Yeoh/Mooney-Rivlin, no linear-shift correction term
is needed to keep this base stress-free at the reference configuration. As
:math:`\lambda_L\to\infty`, :math:`\Psi_{\rm iso}(\bar I_1)\to\frac{\mu}{2}(\bar I_1-3)`, the
isochoric neo-Hookean energy expressed in :math:`\bar I_1` -- this is *not* the same tensor field
as :math:`\Psi_{\rm NeoHooke}` above (which is expressed in the raw :math:`I_1`): the two agree
only at :math:`\boldsymbol{C}=\boldsymbol{I}`, so ArrudaBoyce's full compressible stress response
does not converge to NeoHooke's as :math:`\lambda_L\to\infty` away from the reference
configuration. This isochoric potential is implemented once, in a shared core header
(``Marmot::ContinuumMechanics::EnergyDensityFunctions::ArrudaBoyce8ChainPotential`` and its
``FirstOrderDerived`` counterpart), and reused verbatim by
:doc:`compressiblefinitestrainlinearviscoelasticity`, which adds its own (different) volumetric
convention on top.

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

where :math:`\boldsymbol{\mathcal S}={\rm dev}\left(2\boldsymbol{C}^{\rm e}\,\partial\Psi_B/
\partial\boldsymbol{C}^{\rm e}\right)` is the deviatoric Mandel stress of network B's spring
(computed from the general potential above; for the NeoHooke base this reduces to the closed
form :math:`\boldsymbol{\mathcal S}=\mu_B\,{\rm dev}(\boldsymbol{C}^{\rm e})`, but no such
shortcut exists for Yeoh/Mooney-Rivlin since their potential derivative is not purely linear in
:math:`\boldsymbol{C}^{\rm e}`), :math:`\rho = \sqrt{\boldsymbol{\mathcal S}:\boldsymbol{\mathcal S}}`, and
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
