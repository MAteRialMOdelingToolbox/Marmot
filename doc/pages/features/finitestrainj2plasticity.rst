Finite-strain J2 plasticity
============================

Theory
------

This model is a finite-strain hyperelastic-plastic J2 plasticity formulation with isotropic hardening.
The strain energy potential is based on the compressible Neo-Hookean model (Pence–Gou, variant B).

We assume the multiplicative split of the deformation gradient into elastic and plastic parts,

.. math::

   \boldsymbol{F} = \boldsymbol{F}^{\rm e} \boldsymbol{F}^{\rm p},\qquad
   \boldsymbol{C}^{\rm e} = \big(\boldsymbol{F}^{\rm e})^{\mathsf{T}}\boldsymbol{F}^{\rm e},\qquad
   J^{\rm e} = \det \boldsymbol{F}^{\rm e},\qquad
   I_1^{\rm e} = \operatorname{tr}\boldsymbol{C}^{\rm e}.

The Pence–Gou potential (variant B) applied to :math:`\boldsymbol{C}^{\rm e}` reads:

.. math::

   \Psi(\boldsymbol{C}^{\rm e}) = \Psi_{\rm d}(I_1^{\rm e}, J^{\rm e}) + \Psi_{\rm h}(J^{\rm e}) = \frac{G}{2} \left(I_1^{\rm e} (J^{\rm e})^{-2/3} - 3\right) + \frac{K}{8} \left(J^{\rm e} - (J^{\rm e})^{-1}\right)^2.

From the potential we obtain the second Piola–Kirchhoff stress (elastic configuration)

.. math::

   \boldsymbol{S} = 2 \frac{\partial \Psi}{\partial \boldsymbol{C}^{\rm e}},

and the Mandel stress in the intermediate configuration:

.. math::

   \boldsymbol{M} = \boldsymbol{C}^{\rm e} \boldsymbol{S}.

Plastic admissibility is enforced by a J2 yield function in Mandel stress space with a stress-like hardening variable :math:`\beta_{\rm p}(\alpha_{\rm p})`:

.. math::

   f(\boldsymbol{M}, \beta_{\rm p}) = \frac{1}{f_{\rm y}^0} \left(\sqrt{\boldsymbol{M}^{\rm dev}:\boldsymbol{M}^{\rm dev}} - \sqrt{\tfrac{2}{3}} \beta_{\rm p}\right) \le 0,

where :math:`\boldsymbol{M}^{\rm dev} = \boldsymbol{M} - \frac{1}{3}\operatorname{tr}(\boldsymbol{M})\boldsymbol{I}` is the deviatoric part of the Mandel stress.


The hardening function is:

.. math::

   \beta_{\rm p}(\alpha_{\rm p}) = f_{\rm y}^\infty + (f_{\rm y}^0 - f_{\rm y}^\infty) e^{-\eta \alpha_{\rm p}} + H \alpha_{\rm p},

with a strain-like hardening internal variable :math:`\alpha_{\rm p}`.


An associated flow is used:

.. math::

   \boldsymbol{L}^{\rm p} = \dot{\lambda} \frac{\partial f}{\partial \boldsymbol{M}},
   \qquad
   \dot{\alpha}_{\rm p} = -\dot{\lambda} \frac{\partial f}{\partial \beta_{\rm p}}.


The Karush-Kuhn-Tucker (KKT) conditions for plastic loading are:

.. math::

   \dot{\lambda} \geq 0, \quad f \leq 0, \quad \dot{\lambda} f = 0.

The Kirchhoff stress follows by push-forward from :math:`\boldsymbol{S}` with :math:`\boldsymbol{F}^{\rm e}`:

.. math::

   \boldsymbol{\tau} = \boldsymbol{F}^{\rm e} \, \boldsymbol{S} \, \big(\boldsymbol{F}^{\rm e})^{\mathsf{T}}.

Stress update at a material point proceeds by computing the trial elastic split :math:`\boldsymbol{F}^{\rm e,tr}=\boldsymbol{F}\,\big(\boldsymbol{F}^{\rm p,old}\big)^{-1}` and evaluating the yield function :math:`f`. If plastic yielding occurs (:math:`f > 0`) a fully implicit return-mapping is solved for the unknown vector

.. math::

   \boldsymbol{X} = \begin{bmatrix} \boldsymbol{F}^{\rm e}_{11} \\ \boldsymbol{F}^{\rm e}_{12} \\ \vdots \\ \boldsymbol{F}^{\rm e}_{33} \\ \alpha_{\rm p} \\ \Delta\lambda \end{bmatrix}.

A residual system :math:`\boldsymbol{R}(\boldsymbol{X})=\boldsymbol{0}` enforces (i) elastic–plastic kinematics, (ii) hardening update, and (iii) consistency :math:`f=0`. The plastic update of :math:`\boldsymbol{F}^{\rm p}` uses the matrix exponential of the flow direction in Mandel-space:

.. math::

   \Delta\boldsymbol{F}^{\rm p} = \exp\left(\Delta\lambda\,\frac{\partial f}{\partial \boldsymbol{M}}\right),

and the strain-like hardening variable is updated as:

.. math::

   \Delta\alpha_{\rm p} = -\Delta\lambda\,\frac{\partial f}{\partial \beta_{\rm p}} = \frac{\Delta\lambda}{f_{\rm y}^0} \sqrt{\frac{2}{3}}.

Upon convergence of the return-mapping algorithm, the plastic variables are updated according to:

.. math::

   \boldsymbol{F}^{\rm p,new} = \Delta\boldsymbol{F}^{\rm p}\,\boldsymbol{F}^{\rm p,old},
   \qquad
   \alpha_{\rm p}^{\rm new} = \alpha_{\rm p}^{\rm old} - \Delta\lambda\,\frac{\partial f}{\partial \beta_{\rm p}}
   \;=\; \alpha_{\rm p}^{\rm old} + \frac{\Delta\lambda}{f_{\rm y}^0}\,\sqrt{\tfrac{2}{3}}.

The elastic deformation gradient follows from the multiplicative split as :math:`\boldsymbol{F}^{\rm e,new}=\boldsymbol{F}\,(\boldsymbol{F}^{\rm p,new})^{-1}`, and the Kirchhoff stress and consistent algorithmic tangent are then computed from the hyperelastic constitutive law using :math:`\boldsymbol{F}^{\rm e,new}`.

The consistent tangent :math:`\frac{\partial\boldsymbol{\tau}}{\partial\boldsymbol{F}}` is computed using the chain rule as:

.. math::

   \frac{\partial\boldsymbol{\tau}}{\partial\boldsymbol{F}}
   =
   \frac{\partial\boldsymbol{\tau}}{\partial\boldsymbol{S}} : 2\,\frac{\partial^{2}\Psi}{\partial\boldsymbol{C}^{\rm e}\,\partial\boldsymbol{C}^{\rm e}} : \frac{\partial\boldsymbol{C}^{\rm e}}{\partial\boldsymbol{F}^{\rm e}} : \frac{\partial\boldsymbol{F}^{\rm e}}{\partial\boldsymbol{F}}
   +
   \frac{\partial\boldsymbol{\tau}}{\partial\boldsymbol{F}^{\rm e}} : \frac{\partial\boldsymbol{F}^{\rm e}}{\partial\boldsymbol{F}},

where the individual terms are:

.. math::

   \frac{\partial\boldsymbol{\tau}}{\partial\boldsymbol{S}} = \boldsymbol{F}^{\rm e} \otimes \boldsymbol{F}^{\rm e},

.. math::

   \frac{\partial^{2}\Psi}{\partial\boldsymbol{C}^{\rm e}\,\partial\boldsymbol{C}^{\rm e}}
   =
   \frac{\partial^{2}\Psi}{\partial (J^{\rm e})^{2}}\,
   \frac{\partial J^{\rm e}}{\partial \boldsymbol{C}^{\rm e}}\otimes\frac{\partial J^{\rm e}}{\partial \boldsymbol{C}^{\rm e}}
   +
   \frac{\partial \Psi}{\partial J^{\rm e}}\,
   \frac{\partial^{2} J^{\rm e}}{\partial \boldsymbol{C}^{\rm e}\,\partial \boldsymbol{C}^{\rm e}}
   +
   \frac{\partial^{2}\Psi}{\partial J^{\rm e}\,\partial I_1^{\rm e}}
   \left(
     \frac{\partial J^{\rm e}}{\partial \boldsymbol{C}^{\rm e}}\otimes\frac{\partial I_1^{\rm e}}{\partial \boldsymbol{C}^{\rm e}}
     +
     \frac{\partial I_1^{\rm e}}{\partial \boldsymbol{C}^{\rm e}}\otimes\frac{\partial J^{\rm e}}{\partial \boldsymbol{C}^{\rm e}}
   \right),

.. math::

   \frac{\partial\boldsymbol{C}^{\rm e}}{\partial\boldsymbol{F}^{\rm e}} = (\boldsymbol{F}^{\rm e})^{\mathsf{T}}\otimes\boldsymbol{I} + \boldsymbol{I}\otimes(\boldsymbol{F}^{\rm e})^{\mathsf{T}},

.. math::

   \frac{\partial\boldsymbol{\tau}}{\partial\boldsymbol{F}^{\rm e}} = \frac{\partial}{\partial\boldsymbol{F}^{\rm e}}\left(\boldsymbol{F}^{\rm e}\,\boldsymbol{S}\,(\boldsymbol{F}^{\rm e})^{\mathsf{T}}\right) = \boldsymbol{I}\otimes(\boldsymbol{S}\boldsymbol{F}^{\rm e})^{\mathsf{T}} + (\boldsymbol{S}\boldsymbol{F}^{\rm e})\otimes\boldsymbol{I},

.. math::

   \frac{\partial\boldsymbol{F}^{\rm e}}{\partial\boldsymbol{F}} =
   \begin{cases}
   \boldsymbol{I} \otimes (\boldsymbol{F}^{\rm p,old})^{-\mathsf{T}} & \text{ (for elastic step)} \\[0.5em]
   \text{extracted from } \frac{\partial\boldsymbol{X}}{\partial\boldsymbol{F}} & \text{ (for plastic step)}
   \end{cases}.

For plastic step, :math:`\frac{\partial\boldsymbol{X}}{\partial\boldsymbol{F}}` is computed by solving the linear system :math:`\frac{\partial\boldsymbol{R}}{\partial\boldsymbol{X}} \frac{\partial\boldsymbol{X}}{\partial\boldsymbol{F}} = -\frac{\partial\boldsymbol{R}}{\partial\boldsymbol{F}}` using the converged Jacobian from the return mapping. Then :math:`\frac{\partial\boldsymbol{F}^{\rm e}}{\partial\boldsymbol{F}}` is extracted from the first 9 rows of :math:`\frac{\partial\boldsymbol{X}}{\partial\boldsymbol{F}}`.


.. admonition:: Stress update algorithm at a quadrature point for current step :math:`n+1`

   Notation: :math:`(\cdot)^{\rm old} := (\cdot)_n`, :math:`(\cdot)^{\rm new} := (\cdot)_{n+1}`.

   **Input (known quantities):**

   - Current deformation gradient: :math:`\boldsymbol{F}`
   - Plastic deformation gradient from previous step: :math:`\boldsymbol{F}^{\rm p,old}`
   - Strain-like hardening variable from previous step: :math:`\alpha_{\rm p}^{\rm old}`


   **Trial elastic state:**

   - Compute trial elastic deformation gradient: :math:`\boldsymbol{F}^{\rm e,tr}=\boldsymbol{F}\,\big(\boldsymbol{F}^{\rm p,old}\big)^{-1}`
   - Compute trial Mandel stress: :math:`\boldsymbol{M}^{\rm tr} = \boldsymbol{C}^{\rm e,tr} \boldsymbol{S}^{\rm tr}`
   - Compute stress-like hardening variable: :math:`\beta_{\rm p} = f_{\rm y}^\infty + (f_{\rm y}^0 - f_{\rm y}^\infty) e^{-\eta \alpha_{\rm p}^{\rm old}} + H \alpha_{\rm p}^{\rm old}`
   - Evaluate yield function: :math:`f^{\rm tr} = \frac{1}{f_{\rm y}^0}\left(||\boldsymbol{M}^{\rm tr,dev}|| - \sqrt{\frac{2}{3}} \beta_{\rm p}\right)`

   - **If** :math:`f^{\rm tr} \leq 0` **(elastic step):**

      - Accept the trial state:

        - :math:`\boldsymbol{F}^{\rm e,new} = \boldsymbol{F}^{\rm e,tr}`
        - :math:`\boldsymbol{F}^{\rm p,new} = \boldsymbol{F}^{\rm p,old}`
        - :math:`\alpha_{\rm p}^{\rm new} = \alpha_{\rm p}^{\rm old}`

      - Compute Kirchhoff stress :math:`\boldsymbol{\tau}` and consistent tangent :math:`\frac{\partial\boldsymbol{\tau}}{\partial\boldsymbol{F}}` using :math:`\boldsymbol{F}^{\rm e,new}`
      - **Return**

   **Return mapping:**

   - **If** :math:`f^{\rm tr} > 0` **(plastic step):**

     - Initialize unknowns: :math:`\boldsymbol{X} = \{\boldsymbol{F}^{\rm e}_{11}, \boldsymbol{F}^{\rm e}_{12}, \ldots, \boldsymbol{F}^{\rm e}_{33}, \alpha_{\rm p}, \Delta\lambda\}^{\mathsf{T}}`
     - Solve nonlinear system :math:`\boldsymbol{R}(\boldsymbol{X}) = \boldsymbol{0}` by Newton-Raphson iteration:

       - **Residual equations:**

         - :math:`\boldsymbol{R}_1`: Elastic-plastic kinematics: :math:`\boldsymbol{F}^{\rm e}\,\Delta\boldsymbol{F}^{\rm p} - \boldsymbol{F}\,\big(\boldsymbol{F}^{\rm p,old}\big)^{-1} = \boldsymbol{0}`
         - :math:`\boldsymbol{R}_2`: Hardening evolution: :math:`\alpha_{\rm p} - \alpha_{\rm p}^{\rm old} - \frac{\Delta\lambda}{f_{\rm y}^0} \sqrt{\frac{2}{3}} = 0`
         - :math:`\boldsymbol{R}_3`: Consistency condition: :math:`f(\boldsymbol{M}, \beta_{\rm p}) = 0`

       - **Newton-Raphson iteration:** :math:`\boldsymbol{X}^{(k+1)} = \boldsymbol{X}^{(k)} - \left[\frac{\partial \boldsymbol{R}}{\partial \boldsymbol{X}}\right]^{-1} \boldsymbol{R}(\boldsymbol{X}^{(k)})`

       - **Convergence criteria:** :math:`||\boldsymbol{R}|| < \text{TOL}` and :math:`||\Delta\boldsymbol{X}|| < \text{TOL}`

     - **Upon convergence**, update plastic variables:

       - :math:`\Delta\boldsymbol{F}^{\rm p} = \exp\left(\Delta\lambda \frac{\partial f}{\partial \boldsymbol{M}}\right)`
       - :math:`\boldsymbol{F}^{\rm p,new} = \Delta\boldsymbol{F}^{\rm p} \boldsymbol{F}^{\rm p,old}`
       - :math:`\alpha_{\rm p}^{\rm new} = \alpha_{\rm p}^{\rm old} + \frac{\Delta\lambda}{f_{\rm y}^0} \sqrt{\frac{2}{3}}`

     - Compute updated elastic deformation gradient: :math:`\boldsymbol{F}^{\rm e,new} = \boldsymbol{F}\,(\boldsymbol{F}^{\rm p,new})^{-1}`
     - Compute Kirchhoff stress :math:`\boldsymbol{\tau}` and consistent tangent :math:`\frac{\partial\boldsymbol{\tau}}{\partial\boldsymbol{F}}` using :math:`\boldsymbol{F}^{\rm e,new}`
     - **Return**

Primary reference: A. Dummer, M. Neuner, P. Gamnitzer, G. Hofstetter (2024). Robust and efficient implementation of finite strain generalized continuum models for material failure: Analytical, numerical, and automatic differentiation with hyper-dual numbers. *Computer Methods in Applied Mechanics and Engineering* 426:116987.


Implementation
--------------

.. doxygenclass:: Marmot::Materials::FiniteStrainJ2Plasticity
   :allow-dot-graphs:
