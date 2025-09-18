Finite Strain J2 Plasticity model
=================================

Theory
------

Background
..........
This model combines a Neo-Hookean elastic response with a J2 plasticity mechanism
formulated at finite strains. It targets metals and metal-like behaviors where
pressure sensitivity is small and isotropic hardening is a reasonable first
approximation. The implementation uses standard hyperelasto-plastic ingredients
and a return-mapping scheme in Mandel-stress space.

   * The elastic response uses a compressible Neo-Hookean hyperelastic material  Pence–Gou potential (variant B).
   * Plasticity is formulated in Mandel stress with an associated flow rule. The yield function is
   * \f[
   *   f = \frac{\rho}{f_y} - \sqrt{\frac{2}{3}}\,\frac{\beta_p}{f_y},
   * \f]
   * where \f$\rho=\|\mathrm{dev}(\mathbf{M})\|\f$ and \f$\beta_p\f$ is the current flow stress (isotropic hardening).

Kinematics and elastic response
...............................

We use the multiplicative split

.. math:: \mathbf F = \mathbf F_e \mathbf F_p, \qquad \mathbf C_e = \mathbf F_e^\mathrm T \mathbf F_e,

and the Neo-Hookean elastic energy :math:`\Psi(\mathbf C_e)` as in the elastic model,
with bulk :math:`K` and shear :math:`G`. The Mandel stress is

.. math:: \mathbf M \,=\, \mathbf C_e \,\mathbf S(\mathbf C_e), \quad \mathbf S = 2\,\partial \Psi / \partial \mathbf C_e.

Yield function and hardening
............................
A J2 yield surface in Mandel stress with isotropic hardening:

.. math::

   f(\mathbf M,\beta(\alpha)) \;=\; \bigl\|\mathrm{dev}\,\mathbf M\bigr\|
                                   \;-\; \sqrt{\tfrac{2}{3}}\,\beta(\alpha) \;\le 0,

with accumulated plastic strain :math:`\alpha` and a simple saturation + linear law

.. math::

   \beta(\alpha) \;=\; f_{y\infty} + (f_y - f_{y\infty})\,e^{-\eta\,\alpha} + H\,\alpha.

Flow rule and update
....................
Associated flow is used in Mandel-stress space. The plastic update employs the
exponential map

.. math:: \mathbf F_p^{n+1} \;=\; \exp\!\bigl(\Delta\lambda\,\mathbf N\bigr)\,\mathbf F_p^{n},

with :math:`\mathbf N \propto \partial f/\partial \mathbf M` (normalized flow direction) and
plastic multiplier :math:`\Delta\lambda \ge 0`.

Return mapping (full system)
............................
The unknowns of the local problem are collected in

.. math:: \mathbf X \;=\; \bigl[\mathrm{vec}(\mathbf F_e),\; \alpha,\; \Delta\lambda \bigr]^\mathrm T,

and solved by Newton’s method using the residual blocks

- **Kinematics:** compatibility between :math:`\mathbf F_e`, :math:`\mathbf F_p` and the trial state,
- **Hardening:** :math:`\alpha^{n+1} = \alpha^{n} + \Delta\lambda`,
- **Yield:** :math:`f(\mathbf M(\mathbf F_e), \beta(\alpha)) = 0`.

The consistent (algorithmic) tangent :math:`\partial \boldsymbol\tau / \partial \mathbf F`
is assembled consistently with the same linearization used in the inner Newton solve.

Derivative strategies
.....................
The code provides several interchangeable ways to compute the derivatives needed by
the Newton step and the consistent tangent:

- **Analytic:** closed-form derivatives (baseline).
- **FD (forward/central):** finite-difference perturbations of the residual.
- **CSDA:** complex-step differentiation (machine-precision first derivatives for
  analytic functions).

Material parameters, state, output
..................................
- Parameters: :math:`K, G, f_y, f_{y\infty}, \eta, H`.
- State: :math:`\mathbf F_p` (3×3), :math:`\alpha` (scalar).
- Output per step: Kirchhoff stress :math:`\boldsymbol\tau`, updated :math:`\mathbf F_p`,
  updated :math:`\alpha`, consistent tangent.

Notes
.....

; see, e.g., Hashiguchi & Yamakawa
(2012) for an overview and Dummer et al. (2024) for an implementation-focused summary.

Implementation
--------------

.. doxygenclass:: Marmot::Materials::FiniteStrainJ2Plasticity
   :allow-dot-graphs:
