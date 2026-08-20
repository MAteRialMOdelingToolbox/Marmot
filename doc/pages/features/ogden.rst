.. _ogden:

Ogden model
===========

Theory
------

Hyperelastic materials postulate a strain-energy density function :math:`\Psi` from which stresses
are obtained by differentiation with respect to an appropriate strain measure. For an isotropic
material, :math:`\Psi` may equivalently be expressed in terms of the principal stretches
:math:`\lambda_1,\lambda_2,\lambda_3` of the deformation gradient :math:`\boldsymbol F`, with
:math:`J=\det\boldsymbol F=\lambda_1\lambda_2\lambda_3`.

The classical (multi-term) Ogden model implemented here is a **purely isochoric** hyperelastic
potential, i.e. it contains no volumetric energy term :math:`U(J)`. It is expressed in terms of the
isochoric principal stretches :math:`\bar\lambda_i = J^{-1/3}\lambda_i` (so that
:math:`\bar\lambda_1\bar\lambda_2\bar\lambda_3=1`) as

.. math::

   \Psi_{\rm iso}(\bar\lambda_1,\bar\lambda_2,\bar\lambda_3)
   =
   \sum_{p=1}^{N} \frac{\mu_p}{\alpha_p}
   \left( \bar\lambda_1^{\alpha_p} + \bar\lambda_2^{\alpha_p} + \bar\lambda_3^{\alpha_p} - 3 \right),

where :math:`\mu_p` and :math:`\alpha_p` are the :math:`N` pairs of Ogden moduli and exponents. The
classical incompressible Neo-Hookean (:ref:`incompressibleneohooke`) and Mooney-Rivlin
(:ref:`incompressiblemooneyrivlin`) models are mathematically special one- and two-term cases of
this potential, but - unlike the general (non-integer) Ogden exponents, which require the spectral
decomposition below - are quadratic in the stretches and are implemented directly in terms of the
invariants of :math:`\boldsymbol C`, without spectral decomposition; see those pages for details.

Because the model carries no volumetric stiffness, it cannot be driven by a plain displacement
element: the discrete problem would be under-constrained in volumetric response. It must instead be
used together with a mixed displacement-pressure element that enforces incompressibility as an
independent constraint, e.g. :ref:`displacementpressurefinitestrainelement`, which supplies the
pressure (the Lagrange multiplier for :math:`J=1`) as an additional field.

For an energy expressed in isochoric principal stretches, the principal Kirchhoff stress follows
from the standard result (see e.g. Holzapfel, *Nonlinear Solid Mechanics*, Eq. 6.98)

.. math::

   \tau_i
   =
   \sum_{p=1}^{N} \mu_p\,\bar\lambda_i^{\alpha_p}
   -
   \frac{1}{3}\sum_{j=1}^{3}\sum_{p=1}^{N} \mu_p\,\bar\lambda_j^{\alpha_p}
   ,\qquad i=1,2,3,

which is trace-free (deviatoric) by construction, consistent with the absence of a volumetric term.

Implementation
--------------

The implementation decomposes the right Cauchy-Green tensor :math:`\boldsymbol C=\boldsymbol
F^{\mathsf T}\boldsymbol F` spectrally, :math:`\boldsymbol C = \boldsymbol Q\,\mathrm{diag}(\lambda_i^2)\,
\boldsymbol Q^{\mathsf T}`, evaluates the principal Kirchhoff stress :math:`\tau_i` above, converts
it to a principal second Piola-Kirchhoff stress :math:`S_i=\tau_i/\lambda_i^2`, rotates
:math:`\mathrm{diag}(S_i)` back with :math:`\boldsymbol Q` to obtain :math:`\boldsymbol S`, and
pushes it forward to the Kirchhoff stress :math:`\boldsymbol\tau=\boldsymbol F\,\boldsymbol
S\,\boldsymbol F^{\mathsf T}`. The consistent algorithmic tangent
:math:`\partial\boldsymbol\tau/\partial\boldsymbol F` is obtained automatically via automatic
differentiation (the model derives from ``MarmotMaterialFiniteStrainAD``), so no analytic tangent is
implemented by hand.

.. list-table::
   :header-rows: 1
   :align: left

   * - **Index**
     - **Model Parameter**
     - **Description**
   * - 0
     - :math:`N`
     - Number of Ogden terms
   * - :math:`1 \dots N`
     - :math:`\mu_1 \dots \mu_N`
     - Ogden moduli
   * - :math:`N{+}1 \dots 2N`
     - :math:`\alpha_1 \dots \alpha_N`
     - Ogden exponents
   * - :math:`2N{+}1` (optional)
     - :math:`\rho`
     - Density

.. list-table::
   :header-rows: 1
   :align: left

   * - **State Variable**
     - **Name**
     - **Description**

   * - ---
     - ---
     - No state variables required.

.. doxygenclass:: Marmot::Materials::Ogden
   :allow-dot-graphs:
