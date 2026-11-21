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
and the nonlocal strain follow accordingly. For explicit dynamics, ``computeStressExplicit`` performs the same update
without homogenizing the tangents of the microplanes (one update, with the tangents an optional output): the stress
and the state variables are those of ``computeStress``, bitwise.

The microplane discretization is taken from the shared ``Marmot::ContinuumMechanics::Microplane`` module (37-microplane scheme ``N37``,
fixed at compile time). The closed-form transition function for the elastic unloading stiffness and the gradient
enhancement of the boundaries are introduced in Basche et al. (2026). The constitutive laws of M7 are
implemented once in ``MicroplaneM7Formulation`` and shared with :ref:`MicroplaneM7 <microplanem7>`; GM7 therefore
requires the ``MicroplaneM7`` module.

Properties
----------

The fixed parameters of Caner and Bažant (2013) are calibrated in N and mm, and most of them are dimensionless. The
dimensional ones (:math:`c_{18} = 62.5` MPa, :math:`c_{21} = 250` MPa and the reference values
:math:`f_{c0} = 15.08` MPa, :math:`E_0 = 20000` MPa), the default values of :math:`f_c` and ``transitionWidth`` and
an internal stiffness tolerance are converted to the user's stress unit with the optional property ``stressUnit``,
the value of 1 MPa in that unit (1 for MPa, the default; 1e6 for Pa; 1e-3 for GPa). :math:`E`, :math:`f_c` and
``transitionWidth`` have to be given in the same stress unit; the M7 parameters :math:`\nu`,
:math:`k_1, \ldots, k_5`, :math:`c_3` and ``transitionExponent`` are dimensionless. As ``stressUnit`` follows
:math:`f_c` and ``transitionWidth``, these have to be given explicitly for any unit other than MPa. Note that
``stressUnit`` was inserted at index 13, in front of the properties only required for explicit simulations, so
these moved up by one.

.. list-table::
   :header-rows: 1
   :widths: 10 20 70

   * - Index
     - Name
     - Description
   * - 0
     - ``E``
     - Young's modulus :math:`E`
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
     - Optional, compressive strength :math:`f_c`, used by the deviatoric boundary (default 36 MPa)
   * - 11
     - ``transitionWidth``
     - Optional, stress offset from the virtual elastic loading path at which the unloading transition starts
       (default 4 MPa)
   * - 12
     - ``transitionExponent``
     - Optional, exponent of the unloading transition function (default 0.4)
   * - 13
     - ``stressUnit``
     - Optional, value of 1 MPa in the user's stress unit, e.g. 1e6 for Pa (default 1.0, requires all preceding
       properties)
   * - 14
     - ``density``
     - Optional mass density (requires all preceding properties)
   * - 15
     - ``nonlocalViscosity``
     - Optional viscosity of the nonlocal fields, applied to all six components (requires all preceding
       properties)
   * - 16
     - ``nonlocalMicroInertia``
     - Optional micro-inertia of the nonlocal fields in time squared, applied to all six components (requires all
       preceding properties; default 0, i.e. first-order nonlocal fields). It makes the nonlocal fields second order
       in time, so that the stable increment of an explicit analysis falls off with the element size rather than
       with its square. Admissible up to :math:`\eta^2/4` (with :math:`\eta` the ``nonlocalViscosity``), where the
       reaction mode of the nonlocal fields is critically damped; a larger (underdamped), negative or non-finite value
       is rejected.

The material is registered as ``GM7`` in the general gradient-enhanced hypoelastic material factory (six nonlocal
fields), and is used with the six-field ``GeneralGradientEnhancedDisplacementFiniteElement`` elements (``G6GCPS8``,
``G6GC3D8``, ``G6GC3D20``, ...). At finite strains it is registered as ``GM7/COROTATIONAL`` (same properties) in the
gradient-enhanced finite-strain material factory with six nonlocal fields, for the particles
``GradientEnhanced6FiniteStrain...`` (:doc:`gradientenhancedfinitestrainparticle`): GM7 is then driven by the
logarithmic strain in the frame of the polar rotation, see :doc:`gradientenhancedcorotational`. It is not registered
with the Hughes-Winget wrapper, since it recomputes its stress from its stored strains.

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
