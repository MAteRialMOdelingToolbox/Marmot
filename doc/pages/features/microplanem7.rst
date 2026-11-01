.. _microplanem7:

Microplane model M7
===================

Theory
------

The model implements the microplane model M7 for plain concrete of Caner and Bažant (2013), *Microplane model M7
for plain concrete. I: Formulation* and *II: Calibration and verification*, Journal of Engineering Mechanics,
139(12). The macroscopic strain increment is projected kinematically onto a set of microplanes (normal component
:math:`\varepsilon_N` and the two tangential components :math:`\varepsilon_M`, :math:`\varepsilon_L`), on which an
incremental elastic predictor is limited by the tensile, deviatoric, volumetric and frictional (shear) stress-strain
boundaries of M7. The macroscopic stress follows from the principle of virtual work by homogenization over the
microplanes.

The microplane discretization (normal vectors, weights and the normal-tangential projection tensors) is taken from
the shared ``Marmot::Microplane`` module of ``MarmotMechanicsCore``; the quadrature scheme is fixed at compile time
to the 37-microplane scheme ``N37``.

In addition to the original formulation, a closed-form transition function for the elastic unloading stiffness on
the microplanes is used, which improves the convergence behavior in implicit simulations
(Basche et al., 2026). The consistent tangent is computed analytically.

The constitutive laws are implemented in ``MicroplaneM7Formulation``, which the
:ref:`gradient-enhanced extension GM7 <gm7>` reuses.

The fixed parameters :math:`c_1, \ldots, c_{21}` (except :math:`c_3`) take the values given by Caner and Bažant
(2013).

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
     - ``c3``
     - Optional, slope of the postpeak in uniaxial tension :math:`c_3` (default 4.0)
   * - 8
     - ``fc``
     - Optional, compressive strength :math:`f_c` in N/mm², used by the deviatoric boundary (default 36.0)
   * - 9
     - ``transitionWidth``
     - Optional, stress offset in N/mm² from the virtual elastic loading path at which the unloading transition
       starts (default 4.0)
   * - 10
     - ``transitionExponent``
     - Optional, exponent of the unloading transition function (default 0.4)
   * - 11
     - ``density``
     - Optional mass density (requires all preceding properties)

The material is registered as ``MICROPLANEM7`` in the hypoelastic material factory. The state variables are the
volumetric stress ``sig_v0``, ``zeta``, the microplane stresses ``sig_N0``, ``sig_M0``, ``sig_L0``, the extremal
normal strains ``eps_N0pos``, ``eps_N0neg`` and the macroscopic strain ``strain_macro0``. For postprocessing,
``activeTensileBoundary``, ``activeCompressiveBoundary`` and ``activeShearBoundary`` hold the fractions of the
microplanes on the tensile, compressive and shear boundary in the last increment.

Implementation
--------------

.. doxygenclass:: Marmot::Materials::MicroplaneM7
   :allow-dot-graphs:

.. doxygenclass:: Marmot::Materials::MicroplaneM7Formulation
   :allow-dot-graphs:
