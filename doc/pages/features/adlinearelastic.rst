Linear Elastic model (Automatic Differentiation)
================================================

Theory
------

The constitutive law is given in total form as

.. math::
   \sig = \Cel : \eps


relating the nominal stress tensor :math:`\sig`
to the linearized strain tensor :math:`\eps`
with the fourth order stiffness tensor :math:`\Cel`.
The latter can be specified for isotropic, transversely isotropic
or orthotropic material behavior as follows:

Isotropic Behavior
..................

Number of independent material parameters:   2

.. math::
  \Cel^{-1} = \begin{bmatrix}
                  \frac{1}{E} & \frac{-\nu}{E} & \frac{-\nu}{E} & 0 & 0 & 0 \\
                  \frac{-\nu}{E} & \frac{1}{E} & \frac{-\nu}{E} & 0 & 0 & 0 \\
                  \frac{-\nu}{E} & \frac{-\nu}{E} & \frac{1}{E} & 0 & 0 & 0 \\
               0 & 0 & 0 & \frac{1}{G} & 0 & 0 \\
               0 & 0 & 0 & 0 & \frac{1}{G} & 0 \\
               0 & 0 & 0 & 0 & 0 & \frac{1}{G}
                \end{bmatrix}

with

.. math::
   \displaystyle G = \frac{E}{2\,(1 + \nu)}

Implementation
--------------

.. doxygenclass:: Marmot::Materials::ADLinearElastic
   :allow-dot-graphs:
