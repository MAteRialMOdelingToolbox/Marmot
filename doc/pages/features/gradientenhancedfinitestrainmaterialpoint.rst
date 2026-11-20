Gradient-Enhanced Finite-Strain Material Point
==============================================

Theory
------
The gradient-enhanced finite-strain material point carries the displacement :math:`\boldsymbol u` and
:math:`n_\mathrm{nonlocal}` nonlocal fields :math:`\bar N_K` and drives a material of the interface
``Marmot::Materials::MaterialGradientEnhancedFiniteStrain<nNonlocal>``: with one nonlocal variable e.g.
:doc:`gradientenhancedcompressibleneohookedamage` or :doc:`gradientenhancedfinitestraindruckerprager`, with six (a
nonlocal strain tensor) ``GM7/COROTATIONAL`` (:doc:`gm7`, :doc:`gradientenhancedcorotational`). The material returns
the Kirchhoff stress :math:`\boldsymbol\tau`, the local driving forces :math:`L_K`, the nonlocal radii :math:`R_K`
(with the interactions :math:`c_K = R_K^2`) and the four algorithmic tangents of the implicit-gradient problem

.. math::

   \bar N_K - c_K\,\nabla_X^2 \bar N_K = L_K(\boldsymbol F, \bar{\boldsymbol N}).

Quantities of the nonlocal space are Eigen vectors of size :math:`n_\mathrm{nonlocal}`; a tangent that mixes a spatial
tensor with a nonlocal variable is one tensor per nonlocal variable (``tangents.dS_dN[K]``, ``tangents.dL_dDeltaF[K]``).

It is the material point of the MPM cell :doc:`gradientenhancedfinitestraincell` and is owned by each RKPM particle of
:doc:`gradientenhancedfinitestrainparticle`. The weak forms themselves are assembled there; the material point provides
the kinematics and the constitutive response. See :doc:`meshfreecells` for the material point interface and
:doc:`gradientenhancedfinitestraindisplacementelement` for the finite-element counterpart of the formulation.

Kinematics
^^^^^^^^^^
The material point follows the semi-Lagrangian convention of the meshfree layer. The state stores the deformation
gradient :math:`\boldsymbol F_n = \partial\boldsymbol Y/\partial\boldsymbol X` of the last accepted state (the
*intermediate reference configuration* :math:`\boldsymbol Y`) and the increment
:math:`\Delta\boldsymbol F = \partial\boldsymbol x/\partial\boldsymbol Y` of the current step. The host (a cell or a
particle) interpolates the increments of the current step and adds them with ``incrementDeformation``,

.. math::

   \Delta F_{ij} = \delta_{ij} + \frac{\partial \Delta u_i}{\partial Y_j}, \qquad
   \bar N_K = \bar N_{K,n} + \Delta\bar N_K, \qquad
   \nabla_X \bar N_K \mathrel{+}= \boldsymbol F_n^T \nabla_Y \Delta\bar N_K,

and the material is evaluated at

.. math::

   \boldsymbol F = \Delta\boldsymbol F\,\boldsymbol F_n .

``prepareYourself`` resets :math:`\Delta\boldsymbol u = \boldsymbol 0` and :math:`\Delta\boldsymbol F = \boldsymbol I`
(the nonlocal fields and their gradients are stored as totals and are not reset); ``acceptStateAndPosition`` sets
:math:`\boldsymbol u \leftarrow \boldsymbol u + \Delta\boldsymbol u` and
:math:`\boldsymbol F_n \leftarrow \Delta\boldsymbol F\,\boldsymbol F_n`.

Response and tangents
^^^^^^^^^^^^^^^^^^^^^
After ``computeYourself`` the material point exposes

- ``response.S``: the Kirchhoff stress :math:`\boldsymbol\tau`,
- ``response.dL``: :math:`\Delta L_K`, the driving forces of this evaluation minus the values stored in the state
  variable ``local damage`` (``local strain``), which is then overwritten by the new values. It is the increment of
  :math:`L_K` over the step as long as the host evaluates the point from the state of the last accepted increment,
- ``response.nonLocalRadius``: :math:`R_K`,

and the tangents with respect to the unknowns of the increment. Since :math:`\boldsymbol F = \Delta\boldsymbol F\,
\boldsymbol F_n`, the material tangents with respect to :math:`\boldsymbol F` are converted by the chain rule

.. math::

   \frac{\partial F_{iI}}{\partial \Delta F_{jJ}} = \delta_{ij}\,F_{n,JI}, \qquad
   \frac{\partial\tau_{ij}}{\partial\Delta F_{kL}} = \frac{\partial\tau_{ij}}{\partial F_{mn}}
   \frac{\partial F_{mn}}{\partial\Delta F_{kL}}, \qquad
   \frac{\partial L}{\partial\Delta F_{kL}} = \frac{\partial L}{\partial F_{mn}}
   \frac{\partial F_{mn}}{\partial\Delta F_{kL}},

while :math:`\partial\boldsymbol\tau/\partial\bar N` and :math:`\partial L/\partial\bar N` are passed on unchanged
(fields ``tangents.dS_dDeltaF``, ``tangents.dS_dN``, ``tangents.dL_dDeltaF``, ``tangents.dL_dN``).

In plane strain all quantities are stored in 3D, the material is evaluated through ``computePlaneStrain`` (the 3D
response of a deformation gradient whose out-of-plane entries are those of the identity), and the in-plane components
of the response and the tangents are handed out.

Initialization
^^^^^^^^^^^^^^
``initializeYourself`` sets :math:`\boldsymbol F_n = \boldsymbol I`, a unit eigen deformation and a zero increment,
initializes the material state and queries the reference density from the material, so that the inertia can be
assembled before the first computation. The density is updated after each ``computeYourself``. The undeformed volume
:math:`V_0` and coordinates :math:`\boldsymbol X` are set at construction.

The total gradient :math:`\nabla_X\bar N_K` and the rate :math:`\dot{\bar N}_K` are kept for the explicit scheme of
the particles, whose nonlocal residual is written in total form: the gradient is accumulated by ``incrementDeformation``
in every scheme, so that an explicit step may follow an implicit one.

The initial condition ``geostaticstress`` finds, through the material, the eigen deformation that produces a
hydrostatic stress :math:`\tau_{XX} = \tau_{YY} = \tau_{ZZ} =` ``value[0]`` and applies it in all subsequent material
evaluations.

Registered types
----------------

The material points are registered in the ``Marmot::Registry::MaterialPointFactory``:

.. list-table::
   :header-rows: 1

   * - Name
     - Class
     - Dimension
     - Material evaluation
   * - ``GradientEnhancedFiniteStrain/PlaneStrain``
     - ``GradientEnhancedFiniteStrainMaterialPoint<2, 1>``
     - 2
     - ``computePlaneStrain``
   * - ``GradientEnhancedFiniteStrain/3D``
     - ``GradientEnhancedFiniteStrainMaterialPoint<3, 1>``
     - 3
     - ``computeStress``
   * - ``GradientEnhanced6FiniteStrain/3D``
     - ``GradientEnhancedFiniteStrainMaterialPoint<3, 6>``
     - 3
     - ``computeStress``

The material is created by name from the ``Marmot::Registry::MaterialGradientEnhancedFiniteStrainFactory<nNonlocal>``
(``assignMaterial``), the registry of the materials with the same number of nonlocal variables.

State variables
---------------

The state variables of the material point come first, followed by those of the material (``getStateView`` looks up a
name in the material point first, then in the material). Vectors and tensors always have 3 and 9 entries. The names of
the nonlocal states depend on the number of nonlocal variables: ``nonlocal damage`` and ``local damage`` for one,
``nonlocal strain`` and ``local strain`` for six. Every entry starts at an even offset, so that the mapped tensors stay
16-byte aligned.

.. list-table::
   :header-rows: 1

   * - Name
     - Length
     - Meaning
   * - ``displacement``
     - 3
     - :math:`\boldsymbol u` of the last accepted state
   * - ``velocity``
     - 3
     - velocity (set by the particles' Newmark-beta update)
   * - ``acceleration``
     - 3
     - acceleration (set by the particles' Newmark-beta update)
   * - ``delta displacement``
     - 3
     - :math:`\Delta\boldsymbol u` of the current step
   * - ``delta deformation gradient``
     - 9
     - :math:`\Delta\boldsymbol F = \partial\boldsymbol x/\partial\boldsymbol Y`
   * - ``deformation gradient``
     - 9
     - :math:`\boldsymbol F_n = \partial\boldsymbol Y/\partial\boldsymbol X` of the last accepted state
   * - ``nonlocal damage`` / ``nonlocal strain``
     - :math:`n_\mathrm{nonlocal}`
     - nonlocal fields :math:`\bar N_K` (totals)
   * - ``local damage`` / ``local strain``
     - :math:`n_\mathrm{nonlocal}`
     - local driving forces :math:`L_K` of the last material evaluation
   * - ``stress``
     - 9
     - Kirchhoff stress :math:`\boldsymbol\tau`
   * - ``F0 XX``, ``F0 YY``, ``F0 ZZ``
     - 1 each
     - eigen deformation of a geostatic initial state
   * - ``nonlocal damage gradient`` / ``nonlocal strain gradient``
     - :math:`3\,n_\mathrm{nonlocal}`
     - :math:`\nabla_X\bar N_K` (totals), one column of 3 per nonlocal variable
   * - ``nonlocal damage rate`` / ``nonlocal strain rate``
     - :math:`n_\mathrm{nonlocal}`
     - :math:`\dot{\bar N}_K` of the last explicit increment
   * - ``begin of material state``
     - 0
     - start of the material state variables

Implementation
--------------

.. doxygenclass:: Marmot::Discretizations::MaterialPoints::GradientEnhancedFiniteStrainMaterialPoint
   :allow-dot-graphs:

.. doxygenstruct:: Marmot::Discretizations::MaterialPoints::GradientEnhancedFiniteStrainNames
