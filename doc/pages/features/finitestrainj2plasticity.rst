Finite-strain J2 plasticity
============================

Theory
------

This model is a finite-strain hyperelastic-plastic J2 plasticity formulation with isotropic hardening.
The strain energy potential is based on the compressible Neo-Hookean model (Pence–Gou, variant B).

We assume the multiplicative split of the deformation gradient into elastic and plastic parts,

.. math::

   \mathbf{F} = \mathbf{F}_e \mathbf{F}_p,\qquad
   \mathbf{C}_e = \mathbf{F}_e^{\mathsf{T}}\mathbf{F}_e,\qquad
   J_e = \det \mathbf{F}_e,\qquad
   I_1^e = \operatorname{tr}\mathbf{C}_e.

The Pence–Gou potential (variant B) applied to :math:`\mathbf{C}_e`:

.. math::

   \Psi(\mathbf{C}_e) = \Psi_d(I_1^e, J_e) + \Psi_h(J_e) = \frac{G}{2} \left(I_1^e J_e^{-2/3} - 3\right) + \frac{K}{8} \left(J_e - J_e^{-1}\right)^2.

From the potential we obtain the second Piola–Kirchhoff stress (elastic configuration)

.. math::

   \mathbf{S} = 2 \frac{\partial \Psi}{\partial \mathbf{C}_e},

and the Mandel stress in the intermediate configuration:

.. math::

   \mathbf{M} = \mathbf{C}_e \mathbf{S}.

Plastic admissibility is enforced by a J2 yield function in Mandel stress space with a stress-like hardening variable :math:`\beta_p(\alpha_p)`:

.. math::

   f(\mathbf{M}, \beta_p) = \frac{1}{f_y^0} \left(\sqrt{\mathbf{M}^{\text{dev}}:\mathbf{M}^{\text{dev}}} - \sqrt{\tfrac{2}{3}} \beta_p\right) \le 0,

where :math:`\mathbf{M}^{\text{dev}} = \mathbf{M} - \frac{1}{3}\text{tr}(\mathbf{M})\mathbf{I}` is the deviatoric part of the Mandel stress. The hardening function is:

.. math::

   \beta_p(\alpha_p) = f_y^\infty + (f_y^0 - f_y^\infty) e^{-\eta \alpha_p} + H \alpha_p,

with strain-like hardening internal variable :math:`\alpha_p`. Associated flow is used:

.. math::

   \mathbf{L}_p = \dot{\lambda} \frac{\partial f}{\partial \mathbf{M}},
   \qquad
   \dot{\alpha}_p = -\dot{\lambda} \frac{\partial f}{\partial \beta_p}.


The Kuhn-Tucker conditions for plastic loading are:

.. math::

   \Delta\lambda \geq 0, \quad f \leq 0, \quad \Delta\lambda f = 0.

The Kirchhoff stress follows by push-forward from :math:`\mathbf{S}` with :math:`\mathbf{F}_e`:

.. math::

   \boldsymbol{\tau} = \mathbf{F}_e \, \mathbf{S} \, \mathbf{F}_e^{\mathsf{T}}

Stress update (material point) uses a trial elastic split :math:`\mathbf{F}_e^{\rm tr}=\mathbf{F}\,\mathbf{F}_p^{\rm old\,-1}`, evaluation of :math:`f`, and, if plastic, a fully implicit return-mapping in the unknown vector

.. math::

   \mathbf{X} = \begin{bmatrix} \mathbf{F}_{e,11} \\ \mathbf{F}_{e,12} \\ \vdots \\ \mathbf{F}_{e,33} \\ \alpha_p \\ \Delta\lambda \end{bmatrix},

A residual system :math:`\mathbf{R}(\mathbf{X})=\mathbf{0}` enforces (i) elastic–plastic kinematics, (ii) hardening update, and (iii) consistency :math:`f=0`. The plastic update of :math:`\mathbf{F}_p` uses the matrix exponential of the flow direction in Mandel-space:

.. math::

   \Delta\mathbf{F}_p = \exp\left(\Delta\lambda\,\frac{\partial f}{\partial \mathbf{M}}\right)

.. math::

   \Delta\alpha_p = -\Delta\lambda\,\frac{\partial f}{\partial \beta_p} = \frac{\Delta\lambda}{f_y^0} \sqrt{\frac{2}{3}}.

After convergence, :math:`\mathbf{F}_p^{\rm new} = \Delta\mathbf{F}_p\,\mathbf{F}_p^{\rm old}`, and stresses/tangents are evaluated from the hyperelastic law at :math:`\mathbf{F}_e` consistent with the update.

.. admonition:: Algorithm — Stress update

   **1. Trial state:**

   - Compute trial elastic deformation gradient: :math:`\mathbf{F}_e^{\mathrm{tr}} = \mathbf{F} \mathbf{F}_p^{\mathrm{old},-1}`
   - Evaluate Mandel stress: :math:`\mathbf{M}^{\mathrm{tr}} = \mathbf{C}_e^{\mathrm{tr}} \mathbf{S}^{\mathrm{tr}}`
   - Compute hardening stress: :math:`\beta_p = f_y^\infty + (f_y^0 - f_y^\infty) e^{-\eta \alpha_p^{\mathrm{old}}} + H \alpha_p^{\mathrm{old}}`
   - Evaluate yield function: :math:`f = \frac{1}{f_y^0}\left(||\mathbf{M}^{\mathrm{tr},\text{dev}}|| - \sqrt{\frac{2}{3}} \beta_p\right)`

   **2. Elastic predictor:**

   - If :math:`f \leq 0` (elastic step):

     - Set :math:`\mathbf{F}_p^{\mathrm{new}} = \mathbf{F}_p^{\mathrm{old}}` and :math:`\alpha_p^{\mathrm{new}} = \alpha_p^{\mathrm{old}}`
     - Compute stress and consistent tangent from hyperelastic law
     - **Exit**

   **3. Plastic corrector (return mapping):**

   - If :math:`f > 0` (plastic step):

     - Initialize unknowns: :math:`\mathbf{X} = \{\mathbf{F}_{e,11}, \mathbf{F}_{e,12}, \ldots, \mathbf{F}_{e,33}, \alpha_p, \Delta\lambda\}^{\mathsf{T}}`
     - Solve nonlinear system :math:`\mathbf{R}(\mathbf{X}) = \mathbf{0}` by Newton-Raphson iteration:

       - **Residual equations:**

         - :math:`\mathbf{R}_1`: Elastic-plastic kinematics: :math:`\mathbf{F}_e - \mathbf{F} \mathbf{F}_p^{\mathrm{old}-1} (\Delta\mathbf{F}_p)^{-1} = \mathbf{0}`
         - :math:`\mathbf{R}_2`: Hardening evolution: :math:`\alpha_p - \alpha_p^{\mathrm{old}} - \frac{\Delta\lambda}{f_y^0} \sqrt{\frac{2}{3}} = 0`
         - :math:`\mathbf{R}_3`: Consistency condition: :math:`f(\mathbf{M}, \beta_p) = 0`

       - **Newton iteration:** :math:`\mathbf{X}^{(k+1)} = \mathbf{X}^{(k)} - \left[\frac{\partial \mathbf{R}}{\partial \mathbf{X}}\right]^{-1} \mathbf{R}(\mathbf{X}^{(k)})`

       - **Convergence criteria:** :math:`||\mathbf{R}|| < 10^{-12}` and :math:`||\Delta\mathbf{X}|| < 10^{-12}`

     - Update plastic variables:

       - :math:`\Delta\mathbf{F}_p = \exp\left(\Delta\lambda \frac{\partial f}{\partial \mathbf{M}}\right)`
       - :math:`\mathbf{F}_p^{\mathrm{new}} = \Delta\mathbf{F}_p \mathbf{F}_p^{\mathrm{old}}`
       - :math:`\alpha_p^{\mathrm{new}} = \alpha_p^{\mathrm{old}} + \frac{\Delta\lambda}{f_y^0} \sqrt{\frac{2}{3}}`

   **4. Final stress and tangent computation:**

   - Compute second Piola-Kirchhoff stress: :math:`\mathbf{S} = 2 \frac{\partial \Psi}{\partial \mathbf{C}_e}`
   - Push-forward to Kirchhoff stress: :math:`\boldsymbol{\tau} = \mathbf{F}_e \mathbf{S} \mathbf{F}_e^{\mathsf{T}}`
   - Evaluate consistent algorithmic tangent: :math:`\frac{\partial \boldsymbol{\tau}}{\partial \mathbf{F}}`

Primary reference: A. Dummer, M. Neuner, P. Gamnitzer, G. Hofstetter (2024). Robust and efficient implementation of finite strain generalized continuum models for material failure: Analytical, numerical, and automatic differentiation with hyper-dual numbers. *Computer Methods in Applied Mechanics and Engineering* 426:116987.


Implementation
--------------

.. doxygenclass:: Marmot::Materials::FiniteStrainJ2Plasticity
   :allow-dot-graphs:
