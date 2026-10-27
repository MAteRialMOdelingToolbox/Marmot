.. _gm7:

Gradient-enhanced microplane model M7 (GM7)
===========================================

Theory
------

GM7 is a nonlocal, gradient-enhanced extension of the :ref:`microplane model M7 <microplanem7>` (Caner and Bažant,
2013) for regularized, mesh-independent softening. It is a ``MarmotMaterialGeneralGradientEnhancedHypoElastic``
with six nonlocal fields: the nonlocal macroscopic strain :math:`\bar{\eps}`, which is the solution of the
implicit gradient equation driven by the local strain :math:`\eps` with the interaction parameter :math:`c`.

The elastic predictor on the microplanes is evaluated with the local strain, while the tensile, deviatoric,
volumetric and shear boundaries are evaluated with the *over-nonlocal* strain increment

.. math::

   \Delta \hat{\eps} = m \, \Delta \bar{\eps} + \left( 1 - m \right) \Delta \eps,

in which :math:`m` weights the nonlocal against the local strain. The consistent tangents with respect to the local
and the nonlocal strain follow accordingly.

The microplane discretization is taken from the shared ``Marmot::Microplane`` module (37-microplane scheme ``N37``,
fixed at compile time). The closed-form transition function for the elastic unloading stiffness and the gradient
enhancement of the boundaries are introduced in Basche et al. (2026). The constitutive laws of M7 are
implemented once in ``MicroplaneM7Formulation`` and shared with :ref:`MicroplaneM7 <microplanem7>`; GM7 therefore
requires the ``MicroplaneM7`` module.

Properties
----------

The fixed parameters of Caner and Bažant (2013) are calibrated in N and mm. Therefore :math:`E`, :math:`f_c` and
``transitionWidth`` have to be given in N/mm² (MPa); the M7 parameters :math:`\nu`, :math:`k_1, \ldots, k_5`,
:math:`c_3` and ``transitionExponent`` are dimensionless.

.. list-table::
   :header-rows: 1
   :widths: 10 20 70

   * - Index
     - Name
     - Description
   * - 0
     - ``E``
     - Young's modulus :math:`E` in N/mm²
   * - 1
     - ``nu``
     - Poisson's ratio :math:`\nu`
   * - 2
     - ``k1``
     - Radial scaling parameter :math:`k_1` (typical value 1.5e-4)
   * - 3
     - ``k2``
     - Horizontal asymptote of the frictional boundary :math:`k_2` (typical value 110)
   * - 4
     - ``k3``
     - Shape of the volumetric boundary :math:`k_3` (typical value 30)
   * - 5
     - ``k4``
     - Shape of the volumetric boundary :math:`k_4` (typical value 100)
   * - 6
     - ``k5``
     - Volumetric-deviatoric coupling at low pressures :math:`k_5` (typical value 1e-4)
   * - 7
     - ``c``
     - Nonlocal interaction parameter :math:`c` (length squared, e.g. mm²)
   * - 8
     - ``m``
     - Weighting of the nonlocal strain :math:`m` (0 = local, 1 = nonlocal, :math:`m > 1`: over-nonlocal)
   * - 9
     - ``c3``
     - Optional, slope of the postpeak in uniaxial tension :math:`c_3` (default 4.0)
   * - 10
     - ``fc``
     - Optional, compressive strength :math:`f_c` in N/mm², used by the deviatoric boundary (default 36.0)
   * - 11
     - ``transitionWidth``
     - Optional, stress offset in N/mm² from the virtual elastic loading path at which the unloading transition
       starts (default 4.0)
   * - 12
     - ``transitionExponent``
     - Optional, exponent of the unloading transition function (default 0.4)
   * - 13
     - ``density``
     - Optional mass density (requires all preceding properties)
   * - 14
     - ``nonlocalViscosity``
     - Optional viscosity of the nonlocal fields, applied to all six components (requires all preceding
       properties)

The material is registered as ``GM7`` in the general gradient-enhanced hypoelastic material factory (six nonlocal
fields), and is used with the six-field ``GeneralGradientEnhancedDisplacementFiniteElement`` elements (``G6GCPS8``,
``G6GC3D8``, ``G6GC3D20``, ...).

The state variables are the volumetric stress ``sig_v0``, ``zeta``, the microplane stresses ``sig_N0``, ``sig_M0``,
``sig_L0``, the extremal over-nonlocal normal strains ``eps_N0pos``, ``eps_N0neg``, the macroscopic over-nonlocal
strain ``strain_macro0OverNonLocal`` and the macroscopic local strain ``strain_macro0Local``, which is returned as
the local driving variable of the nonlocal fields. For postprocessing, ``activeTensileBoundary``,
``activeCompressiveBoundary`` and ``activeShearBoundary`` hold the fractions of the microplanes on the tensile,
compressive and shear boundary in the last increment.

Implementation
--------------

.. doxygenclass:: Marmot::Materials::GM7
   :allow-dot-graphs:
