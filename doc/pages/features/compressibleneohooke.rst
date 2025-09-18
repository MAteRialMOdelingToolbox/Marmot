Compressible Neo Hooke model
============================

Theory
------
Background
..........
This is a standard hyperelastic baseline used to model moderately large deformations
with a simple, robust constitutive response. We adopt an isochoric–volumetric split
so that shear and volume changes are controlled independently by shear and bulk
moduli. No history variables are involved (purely elastic), which makes the model
fast and stable; it is also the elastic core used by several inelastic models.

Kinematics
..........
With deformation gradient :math:`\mathbf F`, define the right Cauchy–Green tensor
:math:`\mathbf C=\mathbf F^\mathrm T \mathbf F` and the Jacobian :math:`J=\det\mathbf F`.

Strain energy and stresses
..........................
We use a compressible Neo-Hookean energy of the form

.. math::

   \Psi(\mathbf C) \;=\; \tfrac{K}{2}\,\bigl(J-1\bigr)^2 \;+\; G\bigl(J^{-2/3}\,I_1(\mathbf C)-3\bigr),

with bulk modulus :math:`K`, shear modulus :math:`G`, and
:math:`I_1(\mathbf C)=\mathrm{tr}\,\mathbf C`.
The stresses follow from standard hyperelastic relations:

.. math::

   \mathbf S \,=\, 2\,\frac{\partial \Psi}{\partial \mathbf C}, \qquad
   \boldsymbol\tau \,=\, \mathbf F\,\mathbf S\,\mathbf F^\mathrm T,

where :math:`\mathbf S` is the second Piola–Kirchhoff stress and
:math:`\boldsymbol\tau` the Kirchhoff stress.

Consistent tangent
..................
The FE solver uses the algorithmic (consistent) tangent
:math:`\partial \boldsymbol\tau / \partial \mathbf F` derived from :math:`\Psi(\mathbf C)`.
Because the model is purely elastic, no state updates are performed.

Material parameters and output
..............................
- Parameters: :math:`K` (bulk), :math:`G` (shear).
- State: none.
- Output per step: :math:`\boldsymbol\tau`, elastic energy density, consistent tangent.

Notes
.....
A compact summary of this elastic baseline and its role inside finite-strain updates
is given in Dummer et al. (2024, CMA), which we follow for consistency with the
inelastic framework.


Implementation
--------------

.. doxygenclass:: Marmot::Materials::CompressibleNeoHooke
   :allow-dot-graphs:
