Compressible Finite-Strain Linear Viscoelasticity
=================================================

A finite-strain viscoelastic material model based on a generalized Maxwell model.
The instantaneous hyperelastic response is provided by a selectable compressible
hyperelastic base potential, and the viscous contribution is formulated in terms of
the second Piola–Kirchhoff stress.

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
       ``2`` = MooneyRivlin, ``3`` = PenceGouNeoHooke (variant B)
   * - 1
     - onlyShearCreep
     - Flag: ``1`` restricts viscoelastic creep to the deviatoric (shear) part only,
       ``0`` applies creep to the full stress
   * - 2 … 2+n\ :sub:`e`\ −1
     - elastic properties
     - Elastic properties of the selected hyperelastic base model (see below)
   * - 2+n\ :sub:`e`
     - n\ :sub:`M`
     - Number of Maxwell elements
   * - 2+n\ :sub:`e`\ +1 … 2+n\ :sub:`e`\ +2n\ :sub:`M`
     - :math:`\gamma_i,\,\tau_i`
     - Pairs of relative stiffness weight :math:`\gamma_i` and relaxation time
       :math:`\tau_i` for each Maxwell element (:math:`i=1,\ldots,n_M`)
   * - 2+n\ :sub:`e`\ +1+2n\ :sub:`M` *(optional)*
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
     - :math:`K` (bulk modulus), :math:`G` (shear modulus)
   * - Yeoh
     - 1
     - 4
     - :math:`C_{10}`, :math:`C_{20}`, :math:`C_{30}`, :math:`K`
   * - MooneyRivlin
     - 2
     - 3
     - :math:`C_{10}`, :math:`C_{01}`, :math:`K`
   * - PenceGouNeoHooke
     - 3
     - 2
     - :math:`K` (bulk modulus), :math:`G` (shear modulus)

.. list-table::
   :header-rows: 1
   :align: left

   * - **State Variable**
     - **Description**
   * - ``S0_old`` (9 components)
     - Second Piola–Kirchhoff stress from the previous increment
   * - ``creepStateVars`` (9 × n\ :sub:`M` components)
     - Creep state tensors for each Maxwell element

Theory
------

The model is a finite-strain generalization of the linear viscoelastic model
described, e.g., in Liu et al. (2021). A generalized Maxwell model is employed
in which the long-term elastic part is provided by the chosen hyperelastic base
potential.

The total second Piola–Kirchhoff stress reads

.. math::

   \boldsymbol{S} = \boldsymbol{S}_{\infty} + \sum_{i=1}^{n_M} \boldsymbol{h}_i,

where :math:`\boldsymbol{S}_{\infty}` is the long-term (hyperelastic) stress and
:math:`\boldsymbol{h}_i` are the non-equilibrium (viscous) overstresses of each
Maxwell element.

The evolution of each overstress is governed by

.. math::

   \dot{\boldsymbol{h}}_i + \frac{\boldsymbol{h}_i}{\tau_i}
   = \gamma_i\,\mathbb{C}_0^{-1} : \dot{\boldsymbol{S}}_0,

where :math:`\tau_i` is the relaxation time, :math:`\gamma_i` the relative
stiffness weight, :math:`\mathbb{C}_0` the instantaneous (reference) stiffness
tensor, and :math:`\boldsymbol{S}_0` the instantaneous stress.

The algorithmic update uses a mid-point rule (see Liu et al. 2021):

.. math::

   \boldsymbol{h}_i^{n+1}
   = e^{-\Delta t/\tau_i}\,\boldsymbol{h}_i^n
     + \gamma_i\,\mathbb{C}_0^{-1} :
       \frac{1 - e^{-\Delta t/\tau_i}}{\Delta t/\tau_i}\,\Delta\boldsymbol{S}_0.

Primary reference: Liu et al. 2021, *A continuum and computational framework for
viscoelastodynamics: I. Finite deformation linear models*, Comput. Meth. Appl.
Mech. Engrg. 385, 114059.

Implementation
--------------

.. doxygenclass:: Marmot::Materials::CompressibleFiniteStrainLinearViscoelasticity
   :allow-dot-graphs:
