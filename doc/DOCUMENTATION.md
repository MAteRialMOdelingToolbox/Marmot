## Linear Elastic Model based on Hooke's law

\page linearelastic Linear Elastic

## Implementation

The implementation can be found in \ref Marmot::Materials::LinearElastic

## Theory 

The constitutive law is given in total form as

\f[ 
  \displaystyle \boldsymbol{ \sigma } = \mathbb{ C } : \boldsymbol{ \varepsilon },  
\f] 

relating the nominal stress tensor \f$ \boldsymbol{ \sigma } \f$ 
to the linearized strain tensor\f$ \boldsymbol{ \varepsilon }  \f$ 
with the fourth order stiffness tensor \f$ \mathbb{ C } \f$. The latter is implemented
in \ref MarmotElasticity.h and can be specified for isotropic, transversely isotropic
or orthotropic material behavior as follows:

### Isotropic Behavior
Number of independent material parameters:	2
\f[ 
  \displaystyle \mathbb{ C }^{-1} = \begin{bmatrix}  
				    	\frac{1}{E} & \frac{-\nu}{E} & \frac{-\nu}{E} & 0 & 0 & 0 \\
				    	\frac{-\nu}{E} & \frac{1}{E} & \frac{-\nu}{E} & 0 & 0 & 0 \\
				    	\frac{-\nu}{E} & \frac{-\nu}{E} & \frac{1}{E} & 0 & 0 & 0 \\
					0 & 0 & 0 & \frac{1}{G} & 0 & 0 \\
					0 & 0 & 0 & 0 & \frac{1}{G} & 0 \\
					0 & 0 & 0 & 0 & 0 & \frac{1}{G}
				    \end{bmatrix}
\f] 

with

\f[ 
   \displaystyle G = \frac{E}{2\,(1 + \nu)}
\f] 

### Transversely Isotropic Behavior

In case of transversely isotropic behavior, the user defined normal vector specifies the \f$ x_1 \f$ - axis 
of a local coordinate system, representing the principal material directions in which the material
stiffness tensor is formulated. The isotropic plane is implementend with respect to the local \f$ x_2 \f$ and \f$ x_3 \f$ axes.


Number of independent material parameters:	5

\f[ 
  \displaystyle \mathbb{ C }^{-1} = \begin{bmatrix}  
				    	\frac{1}{E_1} & \frac{-\nu_{12}}{E_2} & \frac{-\nu_{12}}{E_2} & 0 & 0 & 0 \\
				    	\frac{-\nu_{12}}{E_2} & \frac{1}{E_2} & \frac{-\nu_{23}}{E_2} & 0 & 0 & 0 \\
				    	\frac{-\nu_{12}}{E_2} & \frac{-\nu_{12}}{E_2} & \frac{1}{E_2} & 0 & 0 & 0 \\
					0 & 0 & 0 & \frac{1}{G_{12}} & 0 & 0 \\
					0 & 0 & 0 & 0 & \frac{1}{G_{12}} & 0 \\
					0 & 0 & 0 & 0 & 0 & \frac{1}{G_{23}}
				    \end{bmatrix}
\f] 

with

\f[ 
   \displaystyle G_{23} = \frac{E_2}{2\,(1 + \nu_{23})}
\f] 

### Orthotropic Behavior

In case of orthotropic behavior, the user defined normal vector defines the \f$ x_1 \f$ - axis 
of a local coordinate system, representing the principal material directions in which the material stiffness 
tensor is formulated

Number of independent material parameters:	9

\f[ 
  \displaystyle \mathbb{ C }^{-1} = \begin{bmatrix}  
				    	\frac{1}{E_1} & \frac{-\nu_{12}}{E_2} & \frac{-\nu_{13}}{E_3} & 0 & 0 & 0 \\
				    	\frac{-\nu_{12}}{E_2} & \frac{1}{E_2} & \frac{-\nu_{23}}{E_3} & 0 & 0 & 0 \\
				    	\frac{-\nu_{13}}{E_3} & \frac{-\nu_{23}}{E_3} & \frac{1}{E_3} & 0 & 0 & 0 \\
					0 & 0 & 0 & \frac{1}{G_{12}} & 0 & 0 \\
					0 & 0 & 0 & 0 & \frac{1}{G_{13}} & 0 \\
					0 & 0 & 0 & 0 & 0 & \frac{1}{G_{23}}
				    \end{bmatrix}
\f]  
