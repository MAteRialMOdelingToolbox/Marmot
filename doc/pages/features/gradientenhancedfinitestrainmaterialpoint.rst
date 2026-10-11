Gradient-Enhanced Finite-Strain Material Point
==============================================

Theory
------
The gradient-enhanced finite-strain material point carries the displacement :math:`\boldsymbol u` and a single scalar
nonlocal field :math:`\bar N` and drives a material of the interface ``MarmotMaterialGradientEnhancedFiniteStrain``,
e.g. :doc:`gradientenhancedcompressibleneohookedamage` or :doc:`gradientenhancedfinitestraindruckerprager`. The
material returns the Kirchhoff stress :math:`\boldsymbol\tau`, the local driving force :math:`L`, the nonlocal radius
:math:`R` (with the interaction :math:`c = R^2`) and the four algorithmic tangents of the implicit-gradient problem

.. math::

   \bar N - c\,\nabla_X^2 \bar N = L(\boldsymbol F, \bar N).

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
   \bar N = \bar N_n + \Delta\bar N,

and the material is evaluated at

.. math::

   \boldsymbol F = \Delta\boldsymbol F\,\boldsymbol F_n .

``prepareYourself`` resets :math:`\Delta\boldsymbol u = \boldsymbol 0` and :math:`\Delta\boldsymbol F = \boldsymbol I`
(the nonlocal field is stored as a total and is not reset); ``acceptStateAndPosition`` sets
:math:`\boldsymbol u \leftarrow \boldsymbol u + \Delta\boldsymbol u` and
:math:`\boldsymbol F_n \leftarrow \Delta\boldsymbol F\,\boldsymbol F_n`.

Response and tangents
^^^^^^^^^^^^^^^^^^^^^
After ``computeYourself`` the material point exposes

- ``response.S``: the Kirchhoff stress :math:`\boldsymbol\tau`,
- ``response.dL``: :math:`\Delta L`, the driving force of this evaluation minus the value stored in the state variable
  ``local damage``, which is then overwritten by the new value. It is the increment of :math:`L` over the step as long as
  the host evaluates the point from the state of the last accepted increment,
- ``response.nonLocalRadius``: :math:`R`,

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

The initial condition ``geostaticstress`` finds, through the material, the eigen deformation that produces a
hydrostatic stress :math:`\tau_{XX} = \tau_{YY} = \tau_{ZZ} =` ``value[0]`` and applies it in all subsequent material
evaluations.

Registered types
----------------

The material points are registered in the ``MarmotMaterialPointFactory``:

.. list-table::
   :header-rows: 1

   * - Name
     - Class
     - Dimension
     - Material evaluation
   * - ``GradientEnhancedFiniteStrain/PlaneStrain``
     - ``GradientEnhancedFiniteStrainMaterialPoint2D``
     - 2
     - ``computePlaneStrain``
   * - ``GradientEnhancedFiniteStrain/3D``
     - ``GradientEnhancedFiniteStrainMaterialPoint3D``
     - 3
     - ``computeStress``

The material is created by name from the ``MarmotMaterialGradientEnhancedFiniteStrainFactory`` (``assignMaterial``);
a material of another interface is rejected.

State variables
---------------

The state variables of the material point come first, followed by those of the material (``getStateView`` looks up a
name in the material point first, then in the material). Vectors and tensors always have 3 and 9 entries.

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
   * - ``nonlocal damage``
     - 1
     - nonlocal field :math:`\bar N` (total)
   * - ``local damage``
     - 1
     - local driving force :math:`L` of the last material evaluation
   * - ``stress``
     - 9
     - Kirchhoff stress :math:`\boldsymbol\tau`
   * - ``F0 XX``, ``F0 YY``, ``F0 ZZ``
     - 1 each
     - eigen deformation of a geostatic initial state
   * - ``begin of material state``
     - 0
     - start of the material state variables

Implementation
--------------

.. doxygenclass:: Marmot::MaterialPoints::GradientEnhancedFiniteStrainMaterialPoint
   :allow-dot-graphs:

.. doxygenclass:: Marmot::MaterialPoints::GradientEnhancedFiniteStrainMaterialPoint2D
   :allow-dot-graphs:

.. doxygenclass:: Marmot::MaterialPoints::GradientEnhancedFiniteStrainMaterialPoint3D
   :allow-dot-graphs:
