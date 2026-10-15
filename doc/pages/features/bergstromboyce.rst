Bergström-Boyce model
=====================

Theory
------

This model implements the classical Bergström-Boyce finite-strain viscoelastic-viscoplastic
formulation for elastomers/polymers with parallel networks:

- Network A (equilibrium, hyperelastic): sees the total deformation directly.
- Network B ("Maxwell-like", viscous): a spring in series with a viscous flow element.
- Network C ("Maxwell-like", viscous, optional): an additional parallel viscous network,
  allowing multi-mechanism relaxation with distinct elastic and flow characteristics.

All networks share the same selectable hyperelastic base potential
(``hyperelasticBase``), evaluated on the total right Cauchy-Green tensor :math:`\boldsymbol{C}`
for network A (coefficients :math:`A_1,A_2,A_3`, bulk modulus :math:`\kappa_A`) and on the
elastic right Cauchy-Green tensor :math:`\boldsymbol{C}^{\rm e}` for network B's spring
(coefficients :math:`B_1,B_2,B_3`, bulk modulus :math:`\kappa_B`) and network C's spring
(coefficients :math:`C_1,C_2,C_3`, bulk modulus :math:`\kappa_C`). NeoHooke, Yeoh and
Mooney-Rivlin are all built on the isochoric invariants :math:`\bar I_1=I_1\det\boldsymbol{C}^{-1/3}`,
:math:`\bar I_2=I_2\det\boldsymbol{C}^{-2/3}` and share the same volumetric penalty term, and are
therefore all exactly stress-free at :math:`\boldsymbol{C}=\boldsymbol{I}` with no linear-shift
correction needed:

.. math::

   \Psi_{\rm NeoHooke}(\boldsymbol{C}) &= \frac{\mu}{2}\left(\bar I_1 - 3\right)
                        + \frac{\kappa}{8}\left(\ln\det\boldsymbol{C}\right)^2
   \qquad (\text{coefficient 1} = \mu,\ 2,3\text{ unused}) \\[4pt]
   \Psi_{\rm Yeoh}(\boldsymbol{C}) &= C_{10}(\bar I_1-3) + C_{20}(\bar I_1-3)^2 + C_{30}(\bar I_1-3)^3
                        + \frac{\kappa}{8}\left(\ln\det\boldsymbol{C}\right)^2
   \qquad (\text{coefficients} = C_{10},C_{20},C_{30}) \\[4pt]
   \Psi_{\rm MooneyRivlin}(\boldsymbol{C}) &= C_{10}(\bar I_1-3) + C_{01}(\bar I_2-3)
                        + \frac{\kappa}{8}\left(\ln\det\boldsymbol{C}\right)^2
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

   \boldsymbol{F} = \boldsymbol{F}^{\rm e}_B \boldsymbol{F}^{\rm v}_B, \qquad
   \boldsymbol{C}^{\rm e}_B = (\boldsymbol{F}^{\rm e}_B)^{\mathsf T}\boldsymbol{F}^{\rm e}_B,

with a purely symmetric viscous velocity gradient (:math:`\boldsymbol{W}^{\rm v}_B=\boldsymbol 0`)
and a chain-stretch / deviatoric-Mandel-stress power-law flow rule:

.. math::

   \boldsymbol{D}^{\rm v}_B = \dot\gamma_B \boldsymbol{N}_B, \qquad
   \dot\gamma_B = c_1 \left(\lambda^{\rm v}_{{\rm chain},B} - 1\right)^{c_2} \rho_B^{c_3}, \qquad
   \boldsymbol{N}_B = \frac{\boldsymbol{\mathcal S}_B}{\rho_B},

where :math:`\boldsymbol{\mathcal S}_B={\rm dev}\left(2\boldsymbol{C}^{\rm e}_B\,\partial\Psi_B/
\partial\boldsymbol{C}^{\rm e}_B\right)` is the deviatoric Mandel stress of network B's spring,
:math:`\rho_B = \sqrt{\boldsymbol{\mathcal S}_B:\boldsymbol{\mathcal S}_B}`, and
:math:`\lambda^{\rm v}_{{\rm chain},B} = \sqrt{I_1(\boldsymbol{C}^{\rm e}_B)/3}`. Unlike classical
plasticity, there is no yield surface -- the flow rate is a continuous power law that is always
active. Note that :math:`\boldsymbol{N}_B` is traceless by construction, so network B's flow is
isochoric: :math:`\det\boldsymbol{F}^{\rm e}_B(t) = \det\boldsymbol{F}(t)` identically, for every
time :math:`t`.

Network C (optional, activated when :math:`\ge 20` material properties are supplied) acts in parallel
to Network B with its own spring and power-law flow rule:

.. math::

   \boldsymbol{F} = \boldsymbol{F}^{\rm e}_C \boldsymbol{F}^{\rm v}_C, \qquad
   \boldsymbol{C}^{\rm e}_C = (\boldsymbol{F}^{\rm e}_C)^{\mathsf T}\boldsymbol{F}^{\rm e}_C, \\[4pt]
   \boldsymbol{D}^{\rm v}_C = \dot\gamma_C \boldsymbol{N}_C, \qquad
   \dot\gamma_C = c_4 \left(\lambda^{\rm v}_{{\rm chain},C} - 1\right)^{c_5} \rho_C^{c_6}, \qquad
   \boldsymbol{N}_C = \frac{\boldsymbol{\mathcal S}_C}{\rho_C}.

.. admonition:: Stress update algorithm at a quadrature point for current step :math:`n+1`

   **Trial state:**

   - :math:`\boldsymbol{F}^{\rm e,tr}_B = \boldsymbol{F}\,(\boldsymbol{F}^{\rm v,old}_B)^{-1}`
   - If Network C active: :math:`\boldsymbol{F}^{\rm e,tr}_C = \boldsymbol{F}\,(\boldsymbol{F}^{\rm v,old}_C)^{-1}`

   **Return mapping (always active, no yield check):**

   - Solve the decoupled local 10x10 nonlinear system(s) for unknowns
     :math:`\boldsymbol{X}_K=\{\boldsymbol{F}^{\rm e}_{K,11},\ldots,\boldsymbol{F}^{\rm e}_{K,33},\Delta\gamma_K\}^{\mathsf T}`
     (:math:`K = B`, and conditionally :math:`K = C`) by Newton-Raphson iteration:

     - :math:`\boldsymbol{R}_{1,K}`: :math:`\boldsymbol{F}^{\rm e}_K\exp(\Delta\gamma_K\boldsymbol{N}_K) - \boldsymbol{F}^{\rm e,tr}_K = \boldsymbol{0}`
     - :math:`R_{2,K}`: :math:`\Delta\gamma_K/\Delta t - \dot\gamma_K = 0`

     with :math:`\boldsymbol{N}_K`, :math:`\rho_K`, :math:`\lambda^{\rm v}_{{\rm chain},K}` evaluated at
     the current iterate (implicit/backward-Euler).

   - Update viscous deformation gradient(s):
     :math:`\boldsymbol{F}^{\rm v,new}_K = \left[(\boldsymbol{F}^{\rm e,new}_K)^{-1}\boldsymbol{F}^{\rm e,tr}_K\right]\boldsymbol{F}^{\rm v,old}_K`.
   - Total stress: network A's PK2 stress (from :math:`\boldsymbol{C}`) plus network B's PK2
     stress (from :math:`\boldsymbol{C}^{\rm e,new}_B`) and network C's PK2 stress (from
     :math:`\boldsymbol{C}^{\rm e,new}_C`) pushed forward through their respective deformation
     gradients, then superposed as Kirchhoff stress:
     :math:`\boldsymbol\tau = \boldsymbol\tau_A + \boldsymbol\tau_B + \boldsymbol\tau_C`.

Reference: J. S. Bergström, M. C. Boyce (1998). Constitutive modeling of the large strain
time-dependent behavior of elastomers. *Journal of the Mechanics and Physics of Solids*
46(5):931-954.

.. note::

   Only the complex-step-differentiated (CSDA) tangent variant (``implementationType = 0``) is
   currently implemented. The fully analytic tangent (``implementationType = 1``) is reserved
   for future work.

Material parameters
-------------------

The material parameters are provided as a flat vector with the following layout:

.. list-table::
   :header-rows: 1
   :align: left

   * - **Index**
     - **Parameter**
     - **Description**
   * - 0
     - ``hyperelasticBase``
     - Hyperelastic base model selector: ``0`` = NeoHooke, ``1`` = Yeoh,
       ``2`` = MooneyRivlin, ``3`` = ArrudaBoyce
   * - 1
     - :math:`\kappa_A`
     - Bulk modulus of Network A
   * - 2
     - :math:`\kappa_B`
     - Bulk modulus of Network B spring
   * - 3, 4, 5
     - :math:`A_1, A_2, A_3`
     - Network A hyperelastic coefficients (interpretation depends on base model)
   * - 6, 7, 8
     - :math:`B_1, B_2, B_3`
     - Network B hyperelastic coefficients (interpretation depends on base model)
   * - 9
     - :math:`c_1`
     - Network B flow-rate prefactor
   * - 10
     - :math:`c_2`
     - Network B chain-stretch exponent
   * - 11
     - :math:`c_3`
     - Network B effective-stress exponent
   * - 12
     - ``implementationType``
     - Tangent algorithm selector: ``0`` = CSDA (complex-step derivative approximation)
   * - 13 *(2-net)* / 20 *(3-net)* *(optional)*
     - :math:`\rho`
     - Mass density

.. rubric:: Optional Network C parameters (indices 13–19, activated when :math:`\ge 20` parameters provided)

.. list-table::
   :header-rows: 1
   :align: left

   * - **Index**
     - **Parameter**
     - **Description**
   * - 13
     - :math:`\kappa_C`
     - Bulk modulus of Network C spring
   * - 14, 15, 16
     - :math:`C_1, C_2, C_3`
     - Network C hyperelastic coefficients (interpretation depends on base model)
   * - 17
     - :math:`c_4`
     - Network C flow-rate prefactor
   * - 18
     - :math:`c_5`
     - Network C chain-stretch exponent
   * - 19
     - :math:`c_6`
     - Network C effective-stress exponent
   * - 20 *(optional)*
     - :math:`\rho`
     - Mass density

.. rubric:: State variables

.. list-table::
   :header-rows: 1
   :align: left

   * - **State Variable**
     - **Size**
     - **Description**
   * - ``Fv``
     - 9
     - Viscous deformation gradient :math:`\boldsymbol{F}^{\rm v}_B` of Network B
   * - ``FvC`` *(optional)*
     - 9
     - Viscous deformation gradient :math:`\boldsymbol{F}^{\rm v}_C` of Network C (present when Network C is enabled)

Implementation
--------------

.. doxygenclass:: Marmot::Materials::BergstromBoyce
   :allow-dot-graphs:

