.. _gradientEnhancedCorotationalWrapper:

Gradient-enhanced corotational wrapper
======================================

The wrapper makes a small-strain gradient-enhanced material, ``Marmot::Materials::MaterialGradientEnhancedHypoElastic<n>``,
usable as a finite-strain one, ``Marmot::Materials::MaterialGradientEnhancedFiniteStrain<n>``, in the frame of the polar
rotation. It is the wrapper for a material that **recomputes its stress from total strains**, such as :doc:`gm7`, for
which the :ref:`Hughes-Winget wrapper <gradientEnhancedHughesWingetWrapper>` is not objective: that one rotates the
stress it hands over, and GM7 overwrites it.

Theory
------

The deformation gradient is split into the polar rotation and the right stretch,

.. math::

   \boldsymbol F = \boldsymbol R\,\boldsymbol U, \qquad \boldsymbol U = \boldsymbol C^{1/2}, \qquad
   \boldsymbol C = \boldsymbol F^T\boldsymbol F .

The wrapped material lives in the frame that rotates with :math:`\boldsymbol R`. It is driven by the increment of the
logarithmic (Hencky) strain and returns the rotated Kirchhoff stress, which is rotated out:

.. math::

   \Delta\hat{\boldsymbol\varepsilon} = \ln\boldsymbol U^{(n+1)} - \ln\boldsymbol U^{(n)}, \qquad
   \ln\boldsymbol U = \tfrac12 \ln\boldsymbol C, \qquad
   \boldsymbol\tau = \boldsymbol R\,\hat{\boldsymbol\tau}\,\boldsymbol R^T .

- A superposed rigid rotation changes :math:`\boldsymbol R` only. Everything the wrapped material sees is unchanged, so
  the stress rotates exactly, whatever the material does with its state; its tensor-valued internal variables (the
  microplane state of GM7) live in the rotating frame and are objective as well.
- The logarithmic strain increments sum to :math:`\ln\boldsymbol U`, independent of the path. Under a pure stretch the
  material sees :math:`\ln\lambda = (\lambda-1) + O((\lambda-1)^2)`: the small-strain model is recovered up to second
  order.
- :math:`\hat{\boldsymbol\tau}` is work-conjugate to :math:`\ln\boldsymbol U` for an isotropic response, so the energy
  densities of the wrapped material are per unit reference volume and are passed through unchanged.

The nonlocal fields are handed over as they are, together with their increment since the last call. They are referred
to the rotating frame too (the nonlocal strain of GM7 is the nonlocal counterpart of :math:`\hat{\boldsymbol\varepsilon}`),
and their balance is written in the reference configuration by the consumers, so they are invariant under a rigid
rotation.

Tangent
^^^^^^^

The tangent is analytic and consistent,

.. math::

   \frac{\partial\tau_{ij}}{\partial F_{kL}} =
     \frac{\partial R_{ia}}{\partial F_{kL}}\,\hat\tau_{ab}\,R_{jb} + R_{ia}\,\hat\tau_{ab}\,
     \frac{\partial R_{jb}}{\partial F_{kL}} + R_{ia}R_{jb}\,\mathbb C_{abcd}\,
     \frac{\partial(\ln U)_{cd}}{\partial F_{kL}},
   \qquad
   \frac{\partial\boldsymbol R}{\partial F_{kL}} = \Bigl(\frac{\partial\boldsymbol F}{\partial F_{kL}}
     - \boldsymbol R\,\frac{\partial\boldsymbol U}{\partial F_{kL}}\Bigr)\boldsymbol U^{-1},

with :math:`\mathbb C` the algorithmic tangent of the wrapped material. The derivatives of the isotropic tensor
functions :math:`\tfrac12\ln\boldsymbol C` and :math:`\boldsymbol C^{1/2}` follow from the spectral decomposition
:math:`\boldsymbol C = \sum_a \lambda_a \boldsymbol N_a\otimes\boldsymbol N_a` (Daleckii-Krein),

.. math::

   \frac{\partial f(\boldsymbol C)_{ij}}{\partial C_{IJ}} = \sum_{a,b} \theta_{ab}\, N_{a,i} N_{b,j} N_{a,I} N_{b,J},
   \qquad \theta_{ab} = \frac{f(\lambda_a) - f(\lambda_b)}{\lambda_a - \lambda_b}, \quad \theta_{aa} = f'(\lambda_a),

which is continuous where eigenvalues coincide. The divided differences are evaluated without cancellation
(:math:`\ln`: by ``log1p``; the square root: :math:`1/(\sqrt{\lambda_a}+\sqrt{\lambda_b})`). The explicit update,
``computeStressExplicit``, skips the tangent. The tangent of the nonlocal variables is
:math:`\partial\boldsymbol\tau/\partial\bar N_K = \boldsymbol R\,(\partial\hat{\boldsymbol\tau}/\partial\bar N_K)\,
\boldsymbol R^T` and :math:`\partial L_K/\partial\boldsymbol F = \partial L_K/\partial\hat{\boldsymbol\varepsilon} :
\partial\ln\boldsymbol U/\partial\boldsymbol F`.

Usage
-----

The wrapper is registered as ``GM7/COROTATIONAL`` in the ``Marmot::Registry::MaterialGradientEnhancedFiniteStrainFactory<6>``.
Its material properties are those of :doc:`gm7`, unchanged; the wrapper consumes none. It drives the particles
``GradientEnhanced6FiniteStrain...`` of :doc:`gradientenhancedfinitestrainparticle`, e.g. in EdelweissMeshfree:

.. code-block:: python

  material = {"material": "GM7/COROTATIONAL", "properties": np.array([...])}
  particle = MarmotParticleWrapper("GradientEnhanced6FiniteStrainSQCNIxNSNI/3D/Hexa", ..., material)

State variables
^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1

   * - Name
     - Length
     - Meaning
   * - ``Corotational_logU_n``
     - 6
     - :math:`\ln\boldsymbol U` of the last call, Voigt notation with engineering shear
   * - ``Corotational_tau_n``
     - 6
     - rotated Kirchhoff stress :math:`\hat{\boldsymbol\tau}` of the last call, Voigt notation
   * - ``Corotational_N_n``
     - :math:`n`
     - nonlocal fields of the last call
   * - ``materialstate``
     - as the wrapped material
     - its state; its own names are resolved as well

A zero state is the undeformed, stress-free one.

What is verified
----------------

``testGM7Corotational.cpp`` (GM7 module) checks

- the spectral derivative of :math:`\tfrac12\ln\boldsymbol C` against central differences, for distinct and for
  repeated eigenvalues;
- a rigid rotation superposed on every increment of an anisotropic stretch path that runs into softening: the Kirchhoff
  stress rotates to :math:`10^{-10}` of its norm, the local driving forces stay unchanged;
- a pure stretch against small-strain GM7: the difference is of second order in the strain (halving the strain divides
  it by four);
- the analytic tangent (all four blocks) against central differences at a damaged, rotated state, and that the explicit
  update gives the same stress;
- that density, nonlocal viscosity and micro-inertia of GM7 are forwarded and its state names resolved.

Implementation
--------------

.. doxygenclass:: Marmot::Materials::GradientEnhancedCorotationalWrapper
   :allow-dot-graphs:

.. doxygenclass:: Marmot::Materials::GradientEnhancedHypoElasticWrapper
   :allow-dot-graphs:
