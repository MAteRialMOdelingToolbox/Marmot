Meshfree Approximation
======================

Overview
--------

The RKPM particles (:doc:`meshfreeparticles`) take their shape functions from a *meshfree approximation*: a set of
*kernel functions* :math:`\phi_A`, one per node :math:`A` with center :math:`\boldsymbol{x}_A`, is corrected such that
polynomials up to a given *completeness order* :math:`n` are reproduced exactly. ``MarmotMeshfreeCore`` provides

- the interface ``MarmotMeshfreeKernelFunction`` and the boxed tensor-product B-spline kernels of 2nd and 3rd order,
- the interface ``MarmotMeshfreeApproximation`` and the reproducing kernel approximation with explicit
  (``MarmotMeshfreeReproducingKernelApproximation``) and implicit (``MarmotMeshfreeReproducingKernelApproximationImplicit``)
  gradients,
- the complete monomial basis (namespace ``Marmot::Math``),
- B-spline basis functions of arbitrary degree (``MarmotBSpline.h``), used by the B-spline cell geometries of
  :doc:`meshfreecells`.

An approximation is evaluated at a point :math:`\boldsymbol{x}` for a list of :math:`n_\mathrm{c}` candidate kernels. It
returns the values :math:`\Psi_A(\boldsymbol{x})` (an array of length :math:`n_\mathrm{c}`) and the gradients as a
column-major :math:`d \times n_\mathrm{c}` array (entry :math:`i + d\,A` holds :math:`\partial\Psi_A/\partial x_i`),
both indexed by the position of the kernel in the candidate list. Candidates that do not cover :math:`\boldsymbol{x}`
get zero entries. Whether :math:`\boldsymbol{x}` and :math:`\boldsymbol{x}_A` are reference or current coordinates is
decided by the caller (the particle); the approximation only uses their difference.

Theory
------

Reproducing kernel approximation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The reproducing kernel shape function of node :math:`A` is its kernel multiplied by a correction function,

.. math::

   \Psi_A(\boldsymbol{x}) = \boldsymbol{H}^T(\boldsymbol{x} - \boldsymbol{x}_A)\,\boldsymbol{b}(\boldsymbol{x})\,
   \phi_A(\boldsymbol{x}),

where :math:`\boldsymbol{H}` is the vector of all monomials of total degree :math:`\le n` (see
`Monomial basis`_), evaluated at the shifted coordinates :math:`\boldsymbol{x} - \boldsymbol{x}_A`. The coefficients
:math:`\boldsymbol{b}` follow from the *moment matrix* :math:`\boldsymbol{M}`,

.. math::

   \boldsymbol{M}(\boldsymbol{x})\,\boldsymbol{b}(\boldsymbol{x}) = \boldsymbol{H}_0, \qquad
   \boldsymbol{M}(\boldsymbol{x}) = \sum_B \boldsymbol{H}(\boldsymbol{x} - \boldsymbol{x}_B)\,
   \boldsymbol{H}^T(\boldsymbol{x} - \boldsymbol{x}_B)\,\phi_B(\boldsymbol{x}), \qquad
   \boldsymbol{H}_0 = \boldsymbol{H}(\boldsymbol{0}) = [1, 0, \dots, 0]^T,

with the sum over the kernels covering :math:`\boldsymbol{x}` (those with
:math:`|\phi_B(\boldsymbol{x})| > 10^{-14}`). Since :math:`\boldsymbol{M}` is symmetric, this is the familiar form
:math:`\Psi_A = \boldsymbol{H}^T(\boldsymbol{0})\,\boldsymbol{M}^{-1}\boldsymbol{H}(\boldsymbol{x} -
\boldsymbol{x}_A)\,\phi_A`. By construction,
:math:`\sum_A \Psi_A(\boldsymbol{x})\,\boldsymbol{H}(\boldsymbol{x} - \boldsymbol{x}_A) = \boldsymbol{H}_0`, i.e., the
*reproducing conditions*

.. math::

   \sum_A \Psi_A(\boldsymbol{x})\,(\boldsymbol{x} - \boldsymbol{x}_A)^{\boldsymbol{\alpha}} =
   \delta_{\boldsymbol{\alpha}\boldsymbol{0}}
   \quad\Longleftrightarrow\quad
   \sum_A \Psi_A(\boldsymbol{x})\,\boldsymbol{x}_A^{\boldsymbol{\alpha}} = \boldsymbol{x}^{\boldsymbol{\alpha}},
   \qquad |\boldsymbol{\alpha}| \le n,

hold for all multi-indices :math:`\boldsymbol{\alpha}`. In particular, the shape functions form a partition of unity
(:math:`\boldsymbol{\alpha} = \boldsymbol{0}`) and reproduce linear fields for :math:`n \ge 1`. The *completeness
order* :math:`n` is thus the order of the reproduced polynomials, while the *continuity* of :math:`\Psi_A` is that of
the kernels (:math:`C^1` or :math:`C^2` for the kernels below), since :math:`\boldsymbol{H}` is smooth and
:math:`\boldsymbol{M}` inherits the continuity of the kernels. The shape functions are in general not interpolating,
:math:`\Psi_A(\boldsymbol{x}_B) \ne \delta_{AB}`. All linear systems are solved with a column-pivoting Householder QR
decomposition of :math:`\boldsymbol{M}`, computed once per evaluation point.

**Reduction of the completeness order.** :math:`\boldsymbol{M}` is singular if fewer kernels cover
:math:`\boldsymbol{x}` than :math:`\boldsymbol{H}` has entries. At every evaluation point, the requested order is
therefore reduced until the number of covering kernels is at least :math:`\binom{d+n}{d}` (the size of
:math:`\boldsymbol{H}`). This is a necessary condition only; degenerate node arrangements (e.g., all covering nodes on
a line in 2D) can still give a singular moment matrix, so the support radius must be chosen large enough for the
requested order.

Explicit gradients
^^^^^^^^^^^^^^^^^^

``MarmotMeshfreeReproducingKernelApproximation`` computes the exact derivatives of :math:`\Psi_A` by the product rule,
with :math:`(\bullet)_{,i} = \partial(\bullet)/\partial x_i`,

.. math::

   \Psi_{A,i} = \left(\boldsymbol{b}_{,i}^T\boldsymbol{H} + \boldsymbol{b}^T\boldsymbol{H}_{,i}\right)\phi_A
   + \boldsymbol{b}^T\boldsymbol{H}\,\phi_{A,i},
   \qquad
   \boldsymbol{b}_{,i} = -\boldsymbol{M}^{-1}\boldsymbol{M}_{,i}\,\boldsymbol{b},

.. math::

   \boldsymbol{M}_{,i} = \sum_B \left(\boldsymbol{H}_{,i}\boldsymbol{H}^T + \boldsymbol{H}\boldsymbol{H}_{,i}^T\right)
   \phi_B + \boldsymbol{H}\boldsymbol{H}^T\phi_{B,i},

where :math:`\boldsymbol{H}` and :math:`\boldsymbol{H}_{,i}` are evaluated at :math:`\boldsymbol{x} - \boldsymbol{x}_A`
(resp. :math:`\boldsymbol{x} - \boldsymbol{x}_B`). The :math:`d` vectors :math:`\boldsymbol{b}_{,i}` are computed once
per point and reused for all nodes. These gradients are the derivatives of the values, and they reproduce the
derivatives of the polynomials up to order :math:`n`.

Implicit gradients
^^^^^^^^^^^^^^^^^^

``MarmotMeshfreeReproducingKernelApproximationImplicit`` returns the same values :math:`\Psi_A`, but constructs the
gradients directly as corrected kernels (implicit, or synchronized, gradients),

.. math::

   \Psi^{(i)}_A(\boldsymbol{x}) = \boldsymbol{H}^T(\boldsymbol{x} - \boldsymbol{x}_A)\,
   \boldsymbol{b}^{(i)}(\boldsymbol{x})\,\phi_A(\boldsymbol{x}),
   \qquad
   \boldsymbol{M}(\boldsymbol{x})\,\boldsymbol{b}^{(i)}(\boldsymbol{x}) = \boldsymbol{H}^{(i)}_0,
   \qquad
   \boldsymbol{H}^{(i)}_0 = -\left.\frac{\partial\boldsymbol{H}(\boldsymbol{z})}{\partial z_i}\right|_{\boldsymbol{z}
   = \boldsymbol{0}}.

:math:`\boldsymbol{H}^{(i)}_0` is :math:`-1` at the position of the linear monomial :math:`z_i` in
:math:`\boldsymbol{H}` and zero elsewhere; the negative sign stems from the shifted argument
:math:`\boldsymbol{z} = \boldsymbol{x} - \boldsymbol{x}_A`. The position is located from the gradient of the basis at
the origin, so it is correct for every completeness order and dimension (for :math:`n \ge 2` the linear monomials are
not the entries :math:`1, \dots, d` of :math:`\boldsymbol{H}`, see `Monomial basis`_). This enforces the *gradient
reproducing conditions*

.. math::

   \sum_A \Psi^{(i)}_A(\boldsymbol{x})\,\boldsymbol{x}_A^{\boldsymbol{\alpha}} =
   \frac{\partial\boldsymbol{x}^{\boldsymbol{\alpha}}}{\partial x_i}, \qquad |\boldsymbol{\alpha}| \le n.

The implicit gradients reproduce the derivatives of all polynomials up to order :math:`n`, like the explicit ones, but
they are **not** the derivatives of :math:`\Psi_A`. They need neither kernel gradients nor
:math:`\boldsymbol{M}_{,i}`: the :math:`1 + d` right-hand sides :math:`[\boldsymbol{H}_0, \boldsymbol{H}^{(1)}_0,
\dots, \boldsymbol{H}^{(d)}_0]` are solved with a single decomposition of :math:`\boldsymbol{M}`.

In both classes, only ``computeShapeFunctions`` and ``computeShapeFunctionsAndGradients`` are implemented;
``computeShapeFunctionGradients`` throws ``std::runtime_error``.

Kernel functions
^^^^^^^^^^^^^^^^

Both kernels are *boxed* tensor products of a one-dimensional B-spline :math:`w` with support radius :math:`a`,

.. math::

   \phi_A(\boldsymbol{x}) = \prod_{i=1}^{d} w(x_i - x_{A,i}), \qquad z = \frac{|r|}{a},

.. list-table::
   :header-rows: 1

   * - Class
     - :math:`w(r)` for :math:`z \le \tfrac12`
     - :math:`w(r)` for :math:`\tfrac12 < z \le 1`
     - :math:`w(0)`
     - Continuity
   * - ``MarmotMeshfreeKernelFunctionBSpline2ndOrderBoxed``
     - :math:`1 - 2z^2`
     - :math:`2(1-z)^2`
     - :math:`1`
     - :math:`C^1`
   * - ``MarmotMeshfreeKernelFunctionBSpline3rdOrderBoxed``
     - :math:`\tfrac23 - 4z^2 + 4z^3`
     - :math:`\tfrac43(1-z)^3`
     - :math:`\tfrac23`
     - :math:`C^2`

and :math:`w = 0` for :math:`z > 1`. "2nd/3rd order" denotes the polynomial degree: the 2nd order kernel is the
quadratic B-spline on the knots :math:`\{-a, -a/2, a/2, a\}` (scaled to :math:`w(0) = 1`), the 3rd order kernel the
standard cubic B-spline on the uniform knots :math:`\{-a, -a/2, 0, a/2, a\}`. The gradients follow by the product rule,
:math:`\partial\phi_A/\partial x_i = w'(x_i - x_{A,i})\prod_{j\ne i} w(x_j - x_{A,j})`. The kernels are not
normalized; a constant scaling of the kernel cancels in :math:`\Psi_A`.

The support is the open box :math:`|x_i - x_{A,i}| < a` in every direction, i.e., the *support radius* :math:`a` is
half the edge length of the box (not the radius of a sphere) and equal in all directions; ``isInSupport`` tests
:math:`\phi_A > 0`, and ``getBoundingBox`` returns :math:`[\boldsymbol{x}_A - a, \boldsymbol{x}_A + a]`. With a node
spacing :math:`h`, a support radius :math:`a` covers about :math:`2a/h` nodes per direction; the unit tests use
:math:`a = 1.6\,h` for :math:`n = 1` and :math:`a = 2.4\,h` for :math:`n = 2`.

A kernel does **not** own its center: it keeps the pointer passed to its constructor, and ``moveTo`` overwrites the
coordinates behind that pointer. The storage must outlive the kernel.

Monomial basis
^^^^^^^^^^^^^^

``Marmot::Math::computeMonomialBasis`` evaluates the complete basis of all monomials
:math:`x_1^{\alpha_1}\cdots x_d^{\alpha_d}` with :math:`\alpha_1 + \dots + \alpha_d \le n`; it has
``computeSizeOfMonomialBasisVector(n, d)`` :math:`= \binom{d+n}{d}` entries. The exponent of the last coordinate varies
slowest and that of the first coordinate fastest; the first entry is always :math:`1`:

.. list-table::
   :header-rows: 1

   * - :math:`d`
     - :math:`n`
     - :math:`\boldsymbol{H}(\boldsymbol{x})`
   * - 1
     - 2
     - :math:`[1,\ x_1,\ x_1^2]`
   * - 2
     - 1
     - :math:`[1,\ x_1,\ x_2]`
   * - 2
     - 2
     - :math:`[1,\ x_1,\ x_1^2,\ x_2,\ x_1x_2,\ x_2^2]`
   * - 3
     - 1
     - :math:`[1,\ x_1,\ x_2,\ x_3]`

``computeMonomialBasisGradient`` returns the matrix :math:`\partial H_k/\partial x_i` (rows as in
:math:`\boldsymbol{H}`, one column per coordinate), computed by the product rule. Both functions expect the output
already sized. Besides the RK approximation, the basis is used by the variationally consistent integration (VCI) of the
particles (:doc:`meshfreeparticles`).

B-spline basis functions
^^^^^^^^^^^^^^^^^^^^^^^^

``MarmotBSpline.h`` provides the B-spline basis function :math:`N_{i,p}(u)` of degree :math:`p` over a non-decreasing
knot vector :math:`\{z_0, z_1, \dots\}` (``B<p>(u, knotVec, i)``) and its derivative (``dB_dU<p>(u, knotVec, i)``),
evaluated by the Cox--de Boor recursion,

.. math::

   N_{i,0}(u) = \begin{cases} 1, & z_i \le u < z_{i+1}, \\ 0, & \text{otherwise}, \end{cases}
   \qquad
   N_{i,p}(u) = \frac{u - z_i}{z_{i+p} - z_i}N_{i,p-1}(u) + \frac{z_{i+p+1} - u}{z_{i+p+1} - z_{i+1}}N_{i+1,p-1}(u),

where a term with a knot difference below :math:`10^{-14}` (repeated knots) is omitted. The degree is a template
parameter, so the recursion is resolved at compile time. These functions are declared in the global namespace. They are
the basis of the B-spline cell geometries of the material point method (:doc:`meshfreecells`); they are independent of
the kernels above.

Usage
-----

The approximation is not registered with a factory; the host constructs it directly. In EdelweissMeshfree the
approximation types ``ReproducingKernel`` and ``ReproducingKernelImplicitGradient`` (argument ``completenessOrder``)
and the kernel type ``BSplineBoxed`` (arguments ``supportRadius`` and ``continuityOrder`` = 2 or 3, selecting the 2nd
or 3rd order kernel) map to the classes of this page; see :doc:`meshfree` for the overall workflow.

.. list-table::
   :header-rows: 1

   * - Class
     - Constructor arguments
   * - ``MarmotMeshfreeKernelFunctionBSpline2ndOrderBoxed``, ``MarmotMeshfreeKernelFunctionBSpline3rdOrderBoxed``
     - pointer to the center coordinates, dimension :math:`d`, support radius :math:`a`
   * - ``MarmotMeshfreeReproducingKernelApproximation``, ``MarmotMeshfreeReproducingKernelApproximationImplicit``
     - dimension :math:`d`, completeness order :math:`n`

Implementation
--------------

.. doxygenclass:: Marmot::Meshfree::MarmotMeshfreeKernelFunction
   :allow-dot-graphs:

.. doxygenclass:: Marmot::Meshfree::MarmotMeshfreeKernelFunctionBSpline2ndOrderBoxed
   :allow-dot-graphs:

.. doxygenclass:: Marmot::Meshfree::MarmotMeshfreeKernelFunctionBSpline3rdOrderBoxed
   :allow-dot-graphs:

.. doxygenclass:: Marmot::Meshfree::MarmotMeshfreeApproximation
   :allow-dot-graphs:

.. doxygenclass:: Marmot::Meshfree::MarmotMeshfreeReproducingKernelApproximation
   :allow-dot-graphs:

.. doxygenclass:: Marmot::Meshfree::MarmotMeshfreeReproducingKernelApproximationImplicit
   :allow-dot-graphs:

Monomial basis (``MarmotMonomialBasisFunctions.h``, namespace ``Marmot::Math``):

.. doxygenfunction:: Marmot::Math::computeSizeOfMonomialBasisVector

.. doxygenfunction:: Marmot::Math::computeMonomialBasis

.. doxygenfunction:: Marmot::Math::computeMonomialBasisGradient

B-spline basis functions (global namespace):

.. doxygenfile:: MarmotBSpline.h
