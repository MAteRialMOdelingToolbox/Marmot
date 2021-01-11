# MechanicsCore

\page mechanicscore MechanicsCore

## Substepper

### Adaptive Substepper

#### Substepper for Implicit Return Mapping Algorithms

Adaptive Substepper, employing an error estimation and Richardson extrapolation for an implicit return mapping algorithm.

**Implementation:** \ref AdaptiveSubstepper.h

#### Substepper for Semi - Explicit Return Mapping Algorithms
 
Adaptive Substepper, employing an error estimation and Richardson extrapolation for an semi - explicit return mapping algorithm.

**Implementation:** \ref AdaptiveSubstepperExplicit.h

### Perez Fouget Substepper
#### Substepper for Semi - Explicit Return Mapping Algorithms

Perez Fouget Substepper for semi - explicit elastoplastic materials.

**Implementation:** \ref PerezFougetSubstepperExplicitMarkII.h

#### Substepper for Implicit Return Mapping Algorithms

Perez Fouget Substepper for elastoplastic materials with an implicit return mapping algorithm.

**Implementation:** \ref PerezFougetSubstepperMarkII.h

#### Modified Version to Consider Time-Variant Elastic Stiffness Tensor

**Implementation:** \ref PerezFougetSubstepperTime.h

---

## Mechanical Materials
Abstract base class for mechanical materials with scalar nonlocal interaction. 

**Implementation:** \ref MarmotMaterialMechanical

### Hypoelastic Materials

Derived abstract base class for elastic materials expressed purely in rate form. In general, the nominal stress rate tensor \f$ \dot{\boldsymbol{\sigma}} \f$ can be written as a function of the nominal stress tensor \f$ \boldsymbol{\sigma} \f$, the stretching rate tensor \f$ \dot{\boldsymbol{\varepsilon}} \f$ and the time \f$ t \f$.

\f[
  \displaystyle \dot{\boldsymbol{\sigma}} = f( \boldsymbol{\sigma}, \dot{\boldsymbol{\varepsilon}}, t, ...)
\f]

In course of numerical time integration, this relation will be formulated incrementally as 

\f[
  \displaystyle \Delta\boldsymbol{\sigma} = f ( \boldsymbol{\sigma}_n, \Delta\boldsymbol{\varepsilon}, \Delta t, t_n, ...)
\f]

with 

\f[
  \displaystyle \Delta\boldsymbol{\varepsilon} =  \dot{\boldsymbol{\varepsilon}}\cdot \Delta t
\f]

and the algorithmic tangent 

\f[
  \displaystyle \frac{d \boldsymbol{\sigma}}{d \boldsymbol{\varepsilon}} =  \frac{d \Delta \boldsymbol{\sigma}}{d \Delta \boldsymbol{\varepsilon}}
\f]

This formulation is compatible with an Abaqus interface.

**Implementation:** \ref MarmotMaterialHypoElastic

### Hyperelastic Materials

Derived abstract base class for _simple_, purely hyperelastic materials to be used for finite elements based on the total lagrangian kinematic description (TL elements). The second Piola - Kirchhoff stress tensor \f$ S \f$ will be derived by

\f[
  \displaystyle S = \frac{\partial f(\boldsymbol{E},t )}{\partial \boldsymbol{E}}
\f]

with the Green - Lagrange strain tensor \f$ \boldsymbol{E} \f$

\f[
  \displaystyle E  = \frac{1}{2}\,\left(\boldsymbol{F}^T\cdot \boldsymbol{F} - \boldsymbol{I} \right)
\f]

as work conjugated measure and the variable \f$ \boldsymbol{F} \f$ denoting the deformation gradient. 
The algorithmic tangent will be calculated by 

\f[
  \displaystyle \frac{d \boldsymbol{S}}{d \boldsymbol{E}}
\f]



**Implementation:** \ref MarmotMaterialHyperElastic

---

## Gradient Enhancend Mechanical Materials 

Base class for mechanical materials with gradient enhanced regularization to assure mesh indepency in finite element simulations. 

**Implementation:** \ref MarmotMaterialGradientEnhancedMechanical

### Gradient Enhanced Hypoelastic Materials
**Implementation:** \ref MarmotMaterialGradientEnhancedHypoElastic

---
## Mathematical Constants, Functions and Algorithms

**Implementation:** \ref MarmotConstants.h 
		\ref MarmotMath.h

---

## Data Types and Common Tensors

**Implementation:** \ref MarmotTypedefs.h
		\ref MarmotTensor.h

---

## Kinematic Definitions 

**Implementation:** \ref MarmotKinematics.h

---

## Additional
### Duvaut Lions Viscosity

**Implementation:** \ref DuvautLionsViscosity.h

### Convergency Check for Newton Iterations of the Stress Update Algorithm

**Implementation:**  \ref InnerNewtonIterationChecker.h
	 	 \ref InnerNewtonIterationCheckerMarkII.h


### Hughes Winget

**Implementation:** \ref HughesWinget

### Constants for Linear Elastic Material Models 

**Implementation:** \ref MarmotElasticity.h

### Utilities

**Implementation:** \ref MarmotUtility.h

### Voigt Notation

**Implementation:** \ref MarmotVoigt.h

### Menetrey Willam

**Implementation:** \ref MenetreyWillam.h

### Yield Surface Combination Manager for Multisurface Plasticity

**Implementation:** \ref YieldSurfaceCombinationManager.h
