.. _energyGradients:

Energy gradients
================

Gradients of the elastic energy and the dissipation of an element with respect to its degrees of freedom, consistent
with the algorithmic update. They linearize constraints on the energies of a structure, e.g., the dissipation based
arc-length control of EdelweissFE, where an increment dissipates a prescribed amount of energy.

They are **optional capabilities**: elements and materials that provide them inherit from an additional interface
class. The standard interfaces (``Element``, ``MaterialHypoElastic``, ``MaterialGradientEnhancedHypoElastic``) are
unchanged, and elements and materials without the capability behave exactly as before.

Theory
------

The element's elastic energy and dissipation are the sums over its quadrature points,

.. math:: \Psi = \sum_{q} \psi_q\, J_0 w_q, \qquad D = \sum_{q} D_q\, J_0 w_q,

with the energy densities :math:`\psi_q` and :math:`D_q` after the update of the increment. Their gradients follow with
the kinematics of the element, e.g., :math:`\Delta\boldsymbol{\varepsilon} = \boldsymbol{B}\,\Delta\boldsymbol{q}` for
the displacement element,

.. math:: \frac{\partial \Psi}{\partial \boldsymbol{q}}
   = \sum_{q} \boldsymbol{B}^T \frac{\partial \psi_q}{\partial \Delta\boldsymbol{\varepsilon}}\, J_0 w_q, \qquad
   \frac{\partial D}{\partial \boldsymbol{q}}
   = \sum_{q} \boldsymbol{B}^T \frac{\partial D_q}{\partial \Delta\boldsymbol{\varepsilon}}\, J_0 w_q.

The gradient-enhanced element adds the derivatives with respect to the non-local variables, interpolated with the
shape functions of the non-local field. The material provides the derivatives of its energy densities, the
**energy tangents**, consistent with its stress update; they depend on the history, so they are not the derivatives of
a potential.

Plane stress
^^^^^^^^^^^^

In plane stress, the material's update is iterated until :math:`\sigma_f = 0` for the free components
:math:`f = \{33, 13, 23\}` (``Marmot::ContinuumMechanics::ReducedStress``), which determines the strain increments
:math:`\Delta\varepsilon_f` as functions of the in-plane strain increment and the non-local variables. The energy
tangents of the converged 3D update are condensed by the total derivative

.. math:: \frac{\mathrm{d} g}{\mathrm{d} \Delta\varepsilon_{p}}
   = \frac{\partial g}{\partial \Delta\varepsilon_{p}}
   + \frac{\partial g}{\partial \Delta\varepsilon_f} \frac{\partial \Delta\varepsilon_f}{\partial
   \Delta\varepsilon_{p}}, \qquad
   \frac{\mathrm{d} g}{\mathrm{d} \bar\kappa}
   = \frac{\partial g}{\partial \bar\kappa}
   + \frac{\partial g}{\partial \Delta\varepsilon_f} \frac{\partial \Delta\varepsilon_f}{\partial \bar\kappa},

for :math:`g = \psi, D` and the in-plane components :math:`p = \{11, 22, 12\}`, with
:math:`\partial \Delta\varepsilon_f / \partial (\cdot) = -\mathbf{C}_{ff}^{-1}\, \partial \sigma_f / \partial (\cdot)`
from :math:`\mathrm{d}\sigma_f = 0`, :math:`\mathbf{C}_{ff} = \partial \sigma_f / \partial \Delta\varepsilon_f`
(its pseudo-inverse if it is singular, e.g. for a fully softened crack without out-of-plane shear stiffness).

Usage
-----

A framework finds the capability of an element once, assigns two arrays of the size of the element's degrees of
freedom, and reads them after each evaluation of the element:

.. code-block:: cpp

   using namespace Marmot::Discretizations::Elements;

   // after the material section is assigned
   auto* provider = dynamic_cast< EnergyGradientProvider* >( element );
   if ( !provider )
     ; // the element provides no energy gradients

   std::vector< double > dElasticEnergy_dQ( element->getNDofPerElement() );
   std::vector< double > dDissipation_dQ( element->getNDofPerElement() );
   provider->assignEnergyGradientOutput( dElasticEnergy_dQ.data(), dDissipation_dQ.data() );

   element->computeKernels( QTotal, dQ, Pe, Ke, time, dT ); // also writes both gradients

   provider->assignEnergyGradientOutput( nullptr, nullptr ); // stops the output

- **When they are written:** in every call of ``computeKernels()``, for the state of that call, i.e., of the current
  iterate. The arrays are overwritten, not accumulated. ``computeKernelsExplicit()`` does not write them. Elements
  without assigned outputs evaluate their kernels as before.
- **Order:** the element's own order of the degrees of freedom, as the residual ``Pe`` of ``computeKernels()``; the
  framework maps them to its global vector as the residual
  (``Element::getDofIndicesPermutationPattern()``).
- **Checks:** ``assignEnergyGradientOutput()`` throws ``std::invalid_argument`` if the section or the material cannot
  provide the gradients, so an unsupported setup is rejected before the analysis. A material section assigned while
  the outputs are active is checked as well; if it is rejected, the element keeps its material.
- **Lifetime:** the element writes to the arrays until the outputs are reset with ``nullptr``, so they must outlive
  the assignment.

Supported elements and materials
--------------------------------

.. list-table::
   :header-rows: 1
   :align: left

   * - **Element**
     - **Sections**
     - **Material capability**
   * - ``DisplacementFiniteElement``
     - solid, plane strain, axisymmetric
     - ``HypoElasticEnergyTangentProvider``
   * - ``GeneralGradientEnhancedDisplacementFiniteElement``
     - solid, plane strain, plane stress, axisymmetric
     - ``GradientEnhancedHypoElasticEnergyTangentProvider``

Materials with energy tangents: :ref:`vonMisesModel`, Rankine softening plasticity and Rankine damage-plasticity
(``HypoElasticEnergyTangentProvider``), and their gradient-enhanced versions
(``GradientEnhancedHypoElasticEnergyTangentProvider``).

A material provides the energy tangents by inheriting from the capability class next to its base class and
implementing ``computeStressWithEnergyTangents()``, a variant of ``computeStress()`` that also returns the derivatives
of the updated energy densities. Both should share one implementation of the update, so that they give the same
stress. For plane stress, ``GradientEnhancedHypoElasticEnergyTangentProvider::computePlaneStressWithEnergyTangents()``
runs the plane stress iteration of ``MaterialGradientEnhancedHypoElastic::computePlaneStress()`` around the variant; a
material that overrides ``computePlaneStress()`` overrides it consistently.

Implementation
--------------

.. doxygenclass:: Marmot::Discretizations::Elements::EnergyGradientProvider
   :members:

.. doxygenclass:: Marmot::Materials::HypoElasticEnergyTangentProvider
   :members:

.. doxygenclass:: Marmot::Materials::GradientEnhancedHypoElasticEnergyTangentProvider
   :members:
