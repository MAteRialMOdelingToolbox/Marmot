Gradient-Enhanced Finite-Strain Displacement Element
====================================================

Theory
------
This element couples the displacement field :math:`\boldsymbol{u}` (``nDim`` dofs per node) to a single scalar
nonlocal field :math:`\bar N` (1 dof per node) at finite strains. It is the finite-element consumer of
``MarmotMaterialGradientEnhancedFiniteStrain`` (see :doc:`gradientenhancedhugheswinget` and the material models that
implement that interface).

The momentum balance is written with the Kirchhoff stress and the spatial gradients of the shape functions,
integrated over the reference volume,

.. math::

   \boldsymbol{r}_{U,A} = \int_{\Omega_0} \boldsymbol{\tau}\,\nabla_x N_A \,\mathrm{d}V,
   \qquad \nabla_x N_A = \boldsymbol{F}^{-\mathsf T}\nabla_X N_A,

and the nonlocal balance :math:`\bar N - c\,\nabla^2\bar N = L(\boldsymbol F,\bar N)` in the reference
configuration, with the interaction :math:`c = R^2` from the material's nonlocal radius,

.. math::

   r_{N,A} = \int_{\Omega_0} \left( N_A\,\bar N + c\,\nabla_X N_A\cdot\nabla_X\bar N - N_A\,L \right) \mathrm{d}V .

The consistent tangent contains the material blocks :math:`\partial\boldsymbol\tau/\partial\boldsymbol F`,
:math:`\partial\boldsymbol\tau/\partial\bar N`, :math:`\partial L/\partial\boldsymbol F`,
:math:`\partial L/\partial\bar N` and the geometric stiffness of the momentum balance. Plane strain is evaluated as the
3D response with :math:`F_{33} = 1`; plane stress is not supported. Distributed loads: follower pressure (with its
load stiffness) and surface traction; body forces; geostatic initial stresses through the material's eigen
deformation.

Registered elements:

.. list-table::
   :header-rows: 1

   * - Name
     - Shape
     - Integration
     - Section
   * - ``GCPE8UL`` / ``GCPE8RUL``
     - Quad8
     - full / reduced
     - plane strain
   * - ``GC3D8UL``
     - Hexa8
     - full
     - solid
   * - ``GC3D20UL`` / ``GC3D20RUL``
     - Hexa20
     - full / reduced
     - solid

The nodal fields are ``displacement`` and ``nonlocal damage``. A material that stores a plastic deformation gradient
must be initialized (``MarmotMaterialInitialization``, in EdelweissFE ``>>initializematerial``).

Implementation
--------------

.. doxygenclass:: Marmot::Elements::GradientEnhancedFiniteStrainDisplacementElement
   :allow-dot-graphs:
