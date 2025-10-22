Linear Elastic model (Automatic Differentiation)
================================================

This material is intended for demonstration purposes, to illustrate how to use the capabilities of the AD Material.

Theory
------

The constitutive law is given in total form as

.. math::
   \sig = \Cel : \eps


relating the nominal stress tensor :math:`\sig`
to the linearized strain tensor :math:`\eps`
with the fourth order stiffness tensor :math:`\Cel`.
For the given definitions it is important to note, that the Voigt notation is used.
The stiffness tensor can be specified for isotropic material behavior as follows:

Isotropic Behavior
..................

Number of independent material parameters:	2

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
