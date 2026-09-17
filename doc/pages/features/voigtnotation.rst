.. _voigtnotation:

Voigt Notation
==============

Symmetric second order tensors, such as the stress tensor :math:`\sig` and the linearized strain
tensor :math:`\eps`, are stored in Marmot as column vectors, and the fourth order tensors relating
them as matrices. This page defines the component ordering and the scaling factors that are used
consistently throughout the library.

.. warning::
   Marmot orders the shear components as **12, 13, 23** (Abaqus convention). This differs from the
   ordering **23, 13, 12** found in much of the continuum mechanics literature and in several other
   finite element codes. When porting material parameters, stiffness matrices or state variable
   layouts from another code, check the shear ordering first — a mismatch swaps
   :math:`G_{12}` and :math:`G_{23}` and produces results that look plausible but are wrong.

Component ordering
------------------

The vector length follows from the spatial dimension as :math:`n_\mathrm{Voigt} = (n^2 + n)/2`,
provided by the macro ``VOIGTFROMDIM``. The orderings are

.. list-table::
   :header-rows: 1
   :align: left

   * - **Representation**
     - **Size**
     - **Components**
   * - 1D
     - 1
     - :math:`\sigma_{11}`
   * - 2D (plane strain / plane stress)
     - 3
     - :math:`\sigma_{11},\ \sigma_{22},\ \sigma_{12}`
   * - Axisymmetric
     - 4
     - :math:`\sigma_{rr},\ \sigma_{zz},\ \sigma_{\theta\theta},\ \sigma_{rz}`
   * - 3D
     - 6
     - :math:`\sigma_{11},\ \sigma_{22},\ \sigma_{33},\ \sigma_{12},\ \sigma_{13},\ \sigma_{23}`

For the 3D case the mapping between the tensor indices :math:`(i,j)` and the Voigt index is

.. math::
   \begin{bmatrix}
     0 & 3 & 4 \\
     3 & 1 & 5 \\
     4 & 5 & 2
   \end{bmatrix}

and is available in both directions as ``toVoigt<nDim>(i, j)`` and ``fromVoigt<nDim>(ij)`` in the
namespace ``Marmot::ContinuumMechanics::TensorUtility::IndexNotation``. Prefer these helpers over
hard-coded indices — code written against them stays correct if the ordering ever changes.

Stress and strain
-----------------

Stress components are stored as they appear in the tensor, whereas the shear components of the
strain are stored as **engineering shear strains**, i.e. twice the tensor components:

.. math::
   \sig^\mathrm{Voigt} &= \begin{bmatrix}
     \sigma_{11} & \sigma_{22} & \sigma_{33} & \sigma_{12} & \sigma_{13} & \sigma_{23}
   \end{bmatrix}^\mathsf{T} \\[1ex]
   \eps^\mathrm{Voigt} &= \begin{bmatrix}
     \varepsilon_{11} & \varepsilon_{22} & \varepsilon_{33} & \gamma_{12} & \gamma_{13} & \gamma_{23}
   \end{bmatrix}^\mathsf{T}
   = \begin{bmatrix}
     \varepsilon_{11} & \varepsilon_{22} & \varepsilon_{33} & 2\,\varepsilon_{12} & 2\,\varepsilon_{13} & 2\,\varepsilon_{23}
   \end{bmatrix}^\mathsf{T}

This asymmetry is what makes the inner product :math:`\sig^\mathrm{Voigt} \cdot \eps^\mathrm{Voigt}`
equal to the strain energy density :math:`\sig : \eps`. It also means stress and strain vectors must
never be converted with the same routine. The dedicated conversions are

.. list-table::
   :header-rows: 1
   :align: left

   * - **Function**
     - **Purpose**
   * - ``stressToVoigt``, ``voigtToStress``
     - Tensor :math:`\leftrightarrow` Voigt for stress-like quantities (shear copied directly)
   * - ``strainToVoigt``, ``voigtToStrain``
     - Tensor :math:`\leftrightarrow` Voigt for strain-like quantities (shear scaled by 2 resp. 1/2)
   * - ``stressMatrixFromVoigt<nDim>``, ``voigtFromStrainMatrix<nDim>``
     - Dimension-templated variants for the 1D, 2D and 3D representations

The scaling vectors :math:`\mathbf{P} = [1, 1, 1, 2, 2, 2]^\mathsf{T}` and
:math:`\mathbf{P}^{-1} = [1, 1, 1, \tfrac{1}{2}, \tfrac{1}{2}, \tfrac{1}{2}]^\mathsf{T}` are provided
for converting between the two scalings. The constants
:math:`\mathbf{I} = [1, 1, 1, 0, 0, 0]^\mathsf{T}`, :math:`\mathbf{I}_\mathrm{hyd}` and the deviatoric
projector :math:`\mathbf{I}_\mathrm{dev}` are defined in the same namespace.

Fourth order tensors
--------------------

A stiffness tensor in Voigt notation is a :math:`6\times6` matrix that maps the **engineering** strain
vector to the stress vector,

.. math::
   \sig^\mathrm{Voigt} = \Cel^\mathrm{Voigt}\, \eps^\mathrm{Voigt},

so that no additional factors of two appear when the tangent is used. As a consequence, the diagonal
shear entries of an orthotropic stiffness matrix are the shear moduli themselves,

.. math::
   \Cel^\mathrm{Voigt}_{44} = G_{12}, \qquad
   \Cel^\mathrm{Voigt}_{55} = G_{13}, \qquad
   \Cel^\mathrm{Voigt}_{66} = G_{23},

using one-based matrix indices, i.e. entries ``(3,3)``, ``(4,4)`` and ``(5,5)`` in the zero-based C++
code. The corresponding compliance matrix holds :math:`1/G_{12}`, :math:`1/G_{13}` and
:math:`1/G_{23}` in the same positions.

Conversion to and from the full fourth order tensor is available as ``voigtToStiffness`` (Eigen),
``voigtToStiffnessFastor`` (Fastor) and ``stiffnessToVoigt``. These populate

.. math::
   \Cel_{ijkl} = \Cel^\mathrm{Voigt}_{\,\mathrm{toVoigt}(i,j)\,,\,\mathrm{toVoigt}(k,l)},

without halving the shear entries: summing over both :math:`kl` and :math:`lk` reproduces the factor
of two that the engineering shear strain carries in the Voigt product.

Coordinate transformations
--------------------------

Because stress and strain use different shear scalings, they also require different transformation
matrices under a change of basis. Both are built from a local coordinate system in
``Marmot::ContinuumMechanics::VoigtNotation::Transformations``:

- ``transformationMatrixStressVoigt`` — for stress-like Voigt vectors
- ``transformationMatrixStrainVoigt`` — for strain-like Voigt vectors; identical to the stress
  variant except that the upper right :math:`3\times3` block is halved and the lower left block
  doubled

A stiffness matrix is rotated from a local material frame into the global frame as

.. math::
   \Cel^\mathrm{global} = \mathbf{T}_\sigma^{-1}\; \Cel^\mathrm{local}\; \mathbf{T}_\varepsilon,

which is how the transversely isotropic and orthotropic variants of the
:doc:`linearelastic` model are assembled.

Reduced representations
-----------------------

Plane and axisymmetric states are stored in shortened vectors, and the conversions to and from the
full 3D vector are position-based rather than value-based:

.. list-table::
   :header-rows: 1
   :align: left

   * - **Function**
     - **Mapping**
   * - ``voigtToPlaneVoigt``
     - takes entries 0, 1 and 3 of the 3D vector (:math:`11`, :math:`22`, :math:`12`)
   * - ``planeVoigtToVoigt``
     - writes them back; the remaining entries are set to zero
   * - ``voigtToAxisymmetricVoigt``
     - takes entries 0 to 3 (:math:`rr`, :math:`zz`, :math:`\theta\theta`, :math:`rz`)
   * - ``axisymmetricVoigtToVoigt``
     - writes them back; the two out-of-plane shear entries are set to zero

Note that the in-plane shear sits at index 3 of the 3D vector but at index 2 of the plane vector,
which is why the plane-strain and plane-stress tangent reductions in
``Marmot::ContinuumMechanics::PlaneStrain`` and ``Marmot::ContinuumMechanics::PlaneStress`` pick row
and column 3 of the :math:`6\times6` tangent.

.. note::
   ``planeVoigtToVoigt`` assumes that the out-of-plane components are zero. It must not be used to
   restore a state in which :math:`\sigma_{33}` is non-zero, such as a plane strain stress state.

Strain-displacement operators
-----------------------------

The B-operators in ``Marmot::FiniteElement`` produce strain vectors in exactly this ordering. For the
3D case,

.. math::
   \mathbf{B}_a = \begin{bmatrix}
     \partial N_a / \partial x_1 & 0 & 0 \\
     0 & \partial N_a / \partial x_2 & 0 \\
     0 & 0 & \partial N_a / \partial x_3 \\
     \partial N_a / \partial x_2 & \partial N_a / \partial x_1 & 0 \\
     \partial N_a / \partial x_3 & 0 & \partial N_a / \partial x_1 \\
     0 & \partial N_a / \partial x_3 & \partial N_a / \partial x_2
   \end{bmatrix}

where the last three rows again correspond to :math:`\gamma_{12}`, :math:`\gamma_{13}` and
:math:`\gamma_{23}`. The same ordering is used by the B-bar operator, the Green-Lagrange operator
and the enhanced assumed strain interpolations.

Implementation
--------------

The definitions above live in ``MarmotVoigt.h`` and ``MarmotTensor.h``; see
:doc:`../codedocumentation/mechanicscore` and :doc:`../codedocumentation/mathcore` for the full
interface.
