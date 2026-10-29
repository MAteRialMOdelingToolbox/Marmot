.. _linearviscoelasticwiechert:

Linear Viscoelastic Wiechert model
==================================

Theory
------

The model is an isotropic, linear viscoelastic material in the small-strain (hypoelastic) setting.
Its relaxation behaviour is a **power law** that is represented by a generalized Maxwell chain
(Wiechert model) with logarithmically spaced relaxation times.
The stress is obtained from the strain history by the hereditary integral

.. math::

   \sig(t) = \int_0^t \Psi(t - t')\, \DelNu : \dot{\eps}(t')\, \mathrm{d}t',
   \qquad
   \Psi(t) = E_\infty + m\, t^{-n},

with the equilibrium Young's modulus :math:`E_\infty = E`, the power-law relaxation parameter :math:`m`,
the power-law exponent :math:`n`, and the unit stiffness tensor

.. math::

   \DelNu = \Cel(E=1, \nu),

where :math:`\nu` is the Poisson's ratio. The Poisson's ratio is constant in time.

The decaying part :math:`m\,t^{-n}` is replaced by a finite sum of exponentials,

.. math::

   \Psi_N(t) = E + E_0 + \sum_{i=1}^{N} E_i\, e^{-t/\tau_i},

with :math:`N` Maxwell branches. The relaxation times are spaced logarithmically,

.. math::

   \tau_i = \tau_{\min}\, s^{\,i-1}, \qquad s = \sqrt{10}, \qquad i = 1, \dots, N,

with :math:`\tau_{\min}` the smallest relaxation time. The branch moduli :math:`E_i` follow from the
**Post-Widder inversion formula** applied to :math:`m\,t^{-n}` (order :math:`k = 2`), multiplied by the
logarithmic spacing :math:`\ln s`; compare the
:ref:`Linear Viscoelastic Power Law model <linearviscoelasticpowerlaw>`, which uses the same idea
for the creep compliance and a Kelvin chain.
The additional stiffness :math:`E_0` accounts for the part of the spectrum slower than the slowest
resolved branch, so that the tangent stays consistent with the power law beyond
:math:`\tau_N`.

Each branch carries a 6-component stress-like state variable that is updated exponentially over the
time increment :math:`\Delta t`,

.. math::

   \boldsymbol{q}_i^{n+1} = \lambda_i\, E_i\, \DelNu : \Delta\eps + \beta_i\, \boldsymbol{q}_i^{n},
   \qquad \beta_i = e^{-\Delta t/\tau_i},

with :math:`\lambda_i` the standard exponential-integrator coefficient of branch :math:`i`
(see ``computeLambdaAndBeta``). The stress increment and the consistent tangent are

.. math::

   \Delta\sig = \mathbb{C}_{\mathrm{eff}} : \Delta\eps - \sum_i (1-\beta_i)\, \boldsymbol{q}_i^{n},
   \qquad
   \mathbb{C}_{\mathrm{eff}} = \Big(E + E_0 + \sum_i \lambda_i E_i\Big)\, \DelNu .

Since the Poisson's ratio is constant, :math:`\DelNu` is computed once at construction and only
scaled by the effective modulus in every stress update.

Properties
----------

.. list-table::
   :header-rows: 1
   :widths: 10 20 70

   * - Index
     - Name
     - Description
   * - 0
     - ``E``
     - Equilibrium Young's modulus :math:`E`
   * - 1
     - ``nu``
     - Poisson's ratio :math:`\nu`
   * - 2
     - ``m``
     - Power-law relaxation parameter :math:`m \geq 0`
   * - 3
     - ``n``
     - Power-law exponent :math:`n > 0`
   * - 4
     - ``nMaxwell``
     - Number of Maxwell branches :math:`N \geq 1`
   * - 5
     - ``minTau``
     - Smallest relaxation time :math:`\tau_{\min} > 0`
   * - 6
     - ``timeToDays``
     - Factor converting the analysis time to days; the relaxation times are given in days
   * - 7
     - ``density``
     - Optional mass density

The material is registered as ``LINEARVISCOELASTICWIECHERT`` in the hypoelastic material factory.
Because it is a ``MarmotMaterialHypoElastic``, it can also serve as the bulk material of the
:ref:`interface material <interfacematerialhypoelastic>`.
The state variables are the ``6 * nMaxwell`` branch stresses.

Implementation
--------------

.. doxygenclass:: Marmot::Materials::LinearViscoElasticWiechert
   :allow-dot-graphs:
