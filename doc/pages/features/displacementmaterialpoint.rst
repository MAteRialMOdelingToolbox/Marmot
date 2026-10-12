Displacement Material Point
===========================

Theory
------

The displacement material point is the carrier of the material state in the finite-strain displacement formulation of
the material point method (MPM). It consumes any material of the ``MarmotMaterialFiniteStrain`` interface and is
hosted by a :doc:`displacementcell` (MPM) or owned by a :doc:`displacementparticle` (RKPM). The general interfaces of
material points and cells are described in :doc:`meshfreecells`, those of the particles in :doc:`meshfreeparticles`.

Three configurations are distinguished: the undeformed configuration :math:`\boldsymbol{X}`, the intermediate
configuration :math:`\boldsymbol{Y}` (the last accepted state, i.e., the configuration at the beginning of the current
increment) and the current configuration :math:`\boldsymbol{x}`. The deformation gradient is updated
multiplicatively,

.. math::

   \boldsymbol{F} = \frac{\partial \boldsymbol{x}}{\partial \boldsymbol{X}} = \Delta\boldsymbol{F}\,\boldsymbol{F}_n,
   \qquad
   \Delta\boldsymbol{F} = \frac{\partial \boldsymbol{x}}{\partial \boldsymbol{Y}}
      = \boldsymbol{I} + \frac{\partial \Delta\boldsymbol{u}}{\partial \boldsymbol{Y}},
   \qquad
   \boldsymbol{F}_n = \frac{\partial \boldsymbol{Y}}{\partial \boldsymbol{X}},

where the host supplies the displacement increment :math:`\Delta\boldsymbol{u}` and its gradient with respect to
:math:`\boldsymbol{Y}` (``incrementDeformation``; the contributions are accumulated after ``prepareYourself`` has reset
:math:`\Delta\boldsymbol{u} = \boldsymbol{0}`, :math:`\Delta\boldsymbol{F} = \boldsymbol{I}`). ``computeYourself``
evaluates the material with :math:`\boldsymbol{F}` and provides

- the Kirchhoff stress :math:`\boldsymbol{\tau}` (``response.S``),
- its derivative with respect to the increment (``tangents.dS_dDeltaF``),

  .. math::

     \frac{\partial \tau_{ij}}{\partial \Delta F_{kL}} = \frac{\partial \tau_{ij}}{\partial F_{kN}}\,F_{n,LN} ,

  which the hosts contract with the gradients :math:`\partial N_B/\partial Y_L` of their shape functions.

On acceptance of the increment, :math:`\boldsymbol{u} \leftarrow \boldsymbol{u} + \Delta\boldsymbol{u}` and
:math:`\boldsymbol{F}_n \leftarrow \Delta\boldsymbol{F}\,\boldsymbol{F}_n`. The position reported to the host
(``getCoordinatesAtCenter``) is :math:`\boldsymbol{Y} = \boldsymbol{X} + \boldsymbol{u}`, the volume
(``getVolumeUndeformed``) and the density (``getDensityUndeformed``, from the material) refer to the undeformed
configuration. Velocity and acceleration are stored, but set by the host (e.g. by a Newmark-beta scheme in the
particles).

In plane strain, the in-plane increment is expanded to 3D (the out-of-plane components of :math:`\Delta\boldsymbol{F}`
and :math:`\boldsymbol{F}_n` remain those of the identity), the material is evaluated by ``computePlaneStrain`` with the
3D deformation gradient, and the in-plane parts of :math:`\boldsymbol{\tau}` and of its tangent are passed to the host.
The 3D material point evaluates the material with the 3D deformation gradient; the call also goes through
``computePlaneStrain``, whose default implementation in ``MarmotMaterialFiniteStrain`` forwards to ``computeStress``.

Registered material points
--------------------------

.. list-table::
   :header-rows: 1

   * - Name
     - Class
     - Dimension
   * - ``Displacement/PlaneStrain``
     - ``DisplacementMaterialPoint2D``
     - 2D, plane strain
   * - ``Displacement/3D``
     - ``DisplacementMaterialPoint3D``
     - 3D

State variables
---------------

All kinematic states are stored in 3D, regardless of the dimension. They are followed by the state variables of the
material; ``getStateView`` looks up the material point's states first and those of the material otherwise.

.. list-table::
   :header-rows: 1

   * - Name
     - Length
     - Content
   * - ``displacement``
     - 3
     - total displacement :math:`\boldsymbol{u}` of the last accepted state
   * - ``velocity``
     - 3
     - velocity
   * - ``acceleration``
     - 3
     - acceleration
   * - ``delta displacement``
     - 3
     - displacement increment :math:`\Delta\boldsymbol{u}` of the current increment
   * - ``delta deformation gradient``
     - 9
     - incremental deformation gradient :math:`\Delta\boldsymbol{F}`
   * - ``deformation gradient``
     - 9
     - deformation gradient :math:`\boldsymbol{F}_n` of the last accepted state
   * - ``stress``
     - 9
     - Kirchhoff stress :math:`\boldsymbol{\tau}` (written by the plane-strain material point only)
   * - ``begin of material state``
     - --
     - state variables of the material

The material point has no properties of its own; the material is assigned by name and properties
(``assignMaterial``). Initial conditions are not supported (``setInitialCondition`` has no effect).

Implementation
--------------

.. doxygenclass:: Marmot::MaterialPoints::DisplacementMaterialPoint
   :allow-dot-graphs:

.. doxygenclass:: Marmot::MaterialPoints::DisplacementMaterialPoint2D
   :allow-dot-graphs:

.. doxygenclass:: Marmot::MaterialPoints::DisplacementMaterialPoint3D
   :allow-dot-graphs:
