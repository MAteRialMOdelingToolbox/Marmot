Compressible Finite-Strain Fractional Schapery
================================================

A finite-strain nonlinear viscoelastic material model combining a selectable
hyperelastic base potential with a fractional-order power-law relaxation spectrum and
Schapery-type strain-dependent nonlinear modifier functions.

.. rubric:: Formulation

The model formulates the second Piola–Kirchhoff stress in terms of an equilibrium
volumetric response, an equilibrium isochoric response, and a non-equilibrium hereditary
integral:

.. math::

   \boldsymbol{S} = \boldsymbol{S}_{\text{vol}} + g_\infty \boldsymbol{S}_{\text{iso}} + g_1 \boldsymbol{S}_{\text{neq}}

where the non-equilibrium driving stress is modified by Schapery nonlinear strain functions:

.. math::

   \boldsymbol{S}_{\text{eff}} = g_2(\bar{I}_1) \boldsymbol{S}_{\text{iso}}

with the strain-dependent functions defined in terms of the first isochoric invariant :math:`\bar{I}_1 = J^{-2/3} \text{tr}(\boldsymbol{C})`:

.. math::

   g_1(\bar{I}_1) &= \max\left(1 + a_1 (\bar{I}_1 - 3),\, 0.01\right) \\
   g_2(\bar{I}_1) &= \max\left(1 + a_2 (\bar{I}_1 - 3),\, 0.01\right) \\
   a_S(\bar{I}_1) &= \exp\left(-b_S (\bar{I}_1 - 3)\right)

The continuous fractional power-law kernel :math:`E(t) \sim (t/\tau_0)^{-\alpha}` is
discretized into a logarithmic spectrum of 5 Maxwell modes spanning :math:`10^{-2}\tau_0`
to :math:`10^2\tau_0`, integrated numerically in reduced time :math:`\Delta \psi = \Delta t / a_S`.

.. rubric:: Material parameters

The material parameters are provided as a flat vector with the following layout:

.. list-table::
   :header-rows: 1
   :align: left

   * - **Index**
     - **Model Parameter**
     - **Description**
   * - 0
     - hyperelasticBase
     - Hyperelastic base model selector: ``0`` = NeoHooke, ``1`` = Yeoh,
       ``2`` = MooneyRivlin, ``3`` = PenceGouNeoHooke, ``4`` = ArrudaBoyce
   * - 1
     - onlyShearCreep
     - Flag: ``1`` restricts viscoelastic creep to the deviatoric/isochoric part only,
       ``0`` applies creep to the full stress
   * - 2 … 2+n\ :sub:`e`\ −1
     - elastic properties
     - Elastic properties of the selected hyperelastic base model (see below)
   * - 2+n\ :sub:`e`
     - :math:`\alpha`
     - Fractional power-law relaxation exponent (:math:`0 < \alpha < 1`)
   * - 2+n\ :sub:`e`\ +1
     - :math:`\tau_0`
     - Reference relaxation time scale in seconds
   * - 2+n\ :sub:`e`\ +2
     - :math:`g_\infty`
     - Long-term equilibrium stiffness ratio (:math:`0 < g_\infty < 1`)
   * - 2+n\ :sub:`e`\ +3
     - :math:`a_1`
     - Schapery :math:`g_1` strain nonlinearity parameter
   * - 2+n\ :sub:`e`\ +4
     - :math:`a_2`
     - Schapery :math:`g_2` strain nonlinearity parameter
   * - 2+n\ :sub:`e`\ +5
     - :math:`b_S`
     - Schapery :math:`a_S` strain-time shift parameter
   * - 2+n\ :sub:`e`\ +6 *(optional)*
     - :math:`\rho`
     - Density

.. rubric:: Elastic properties by base model

.. list-table::
   :header-rows: 1
   :align: left

   * - **Base model**
     - **selector value**
     - **n**\ :sub:`e`
     - **Properties (in order)**
   * - NeoHooke
     - 0
     - 2
     - :math:`\mu` (shear-modulus-like parameter), :math:`\kappa` (bulk modulus)
   * - Yeoh
     - 1
     - 4
     - :math:`C_{10}`, :math:`C_{20}`, :math:`C_{30}`, :math:`\kappa`
   * - MooneyRivlin
     - 2
     - 3
     - :math:`C_{10}`, :math:`C_{01}`, :math:`\kappa`
   * - PenceGouNeoHooke
     - 3
     - 2
     - :math:`K` (bulk modulus), :math:`G` (shear modulus)
   * - ArrudaBoyce
     - 4
     - 3
     - :math:`\mu`, :math:`\lambda_L`, :math:`\kappa`

.. rubric:: State variables

.. list-table::
   :header-rows: 1
   :align: left

   * - **State Variable**
     - **Description**
   * - ``S0_old`` (9 components)
     - Previous step effective driving stress :math:`\boldsymbol{S}_{\text{eff}}`
   * - ``creepStateVars`` (45 components)
     - Internal history tensors for the 5 discrete Maxwell modes
