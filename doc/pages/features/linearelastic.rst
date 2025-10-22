Linear Elastic model
====================

Theory
------

The constitutive law is given in total form as

.. math::
   \sig = \Cel : \eps


relating the nominal stress tensor :math:`\sig`
to the linearized strain tensor :math:`\eps`
with the fourth order stiffness tensor :math:`\Cel`.
For the given definitions it is important to note, that the Voigt notation is used.
The stiffness tensor can be specified for isotropic, transversely isotropic
or orthotropic material behavior as follows:

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

Transversely Isotropic Behavior
...............................
In case of transversely isotropic behavior, the user defined normal vector specifies the :math:`x_1` - axis
of a local coordinate system, representing the principal material directions in which the material
stiffness tensor is formulated. The isotropic plane is implementend with respect to the local :math:`x_2` and :math:`x_3` axes.
The two required vectors must be orthogonal to each other.

Number of independent material parameters:	5

.. math::
  \Cel^{-1} = \begin{bmatrix}
				    	\frac{1}{E_1} & \frac{-\nu_{12}}{E_2} & \frac{-\nu_{12}}{E_2} & 0 & 0 & 0 \\
				    	\frac{-\nu_{12}}{E_2} & \frac{1}{E_2} & \frac{-\nu_{23}}{E_2} & 0 & 0 & 0 \\
				    	\frac{-\nu_{12}}{E_2} & \frac{-\nu_{12}}{E_2} & \frac{1}{E_2} & 0 & 0 & 0 \\
					0 & 0 & 0 & \frac{1}{G_{12}} & 0 & 0 \\
					0 & 0 & 0 & 0 & \frac{1}{G_{12}} & 0 \\
					0 & 0 & 0 & 0 & 0 & \frac{1}{G_{23}}
				    \end{bmatrix}

with

.. math::
   G_{23} = \frac{E_2}{2\,(1 + \nu_{23})}

Orthotropic Behavior
....................

In case of orthotropic behavior, the user defined normal vector defines the :math:`x_1` - axis
of a local coordinate system, representing the principal material directions in which the material stiffness
tensor is formulated.
The two required vectors must be orthogonal to each other.

Number of independent material parameters:	9

.. math::
  \mathbb{ C }^{-1} = \begin{bmatrix}
				    	\frac{1}{E_1} & \frac{-\nu_{12}}{E_2} & \frac{-\nu_{13}}{E_3} & 0 & 0 & 0 \\
				    	\frac{-\nu_{12}}{E_2} & \frac{1}{E_2} & \frac{-\nu_{23}}{E_3} & 0 & 0 & 0 \\
				    	\frac{-\nu_{13}}{E_3} & \frac{-\nu_{23}}{E_3} & \frac{1}{E_3} & 0 & 0 & 0 \\
					0 & 0 & 0 & \frac{1}{G_{12}} & 0 & 0 \\
					0 & 0 & 0 & 0 & \frac{1}{G_{13}} & 0 \\
					0 & 0 & 0 & 0 & 0 & \frac{1}{G_{23}}
				    \end{bmatrix}


Implementation
--------------

.. doxygenclass:: Marmot::Materials::LinearElastic
   :allow-dot-graphs:
