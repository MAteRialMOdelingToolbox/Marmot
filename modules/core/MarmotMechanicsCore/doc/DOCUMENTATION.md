## MechanicsCore

\page continuummechanicsothers Others

### Additional implementations for Material Models based on Plasticity Theory

#### Duvaut Lions Viscosity

**Implementation:** \ref DuvautLionsViscosity.h

#### Menetrey Willam Yield Surfaces

**Implementation:** \ref MenetreyWillam.h


#### Yield Surface Combination Manager

Necessary for multisurface plasticity material models.

**Implementation:** \ref YieldSurfaceCombinationManager.h

\page mechanicalmaterials Mechanical Material Models

Abstract base class for mechanical materials with scalar nonlocal interaction.

**Implementation:** \ref MarmotMaterialMechanical

\page hypoelastic Hypoelastic Material Models

**Implementation:** \ref MarmotMaterialHypoElastic

Derived abstract base class for elastic materials expressed purely in rate form.

## Basic Theory

In general, the nominal stress rate tensor \f$ \sigRate \f$ can be written as a function of the nominal stress tensor \f$ \sig \f$, the stretching rate tensor \f$ \epsRate \f$ and the time \f$ t \f$.

\f[  \displaystyle \sigRate = f( \sig, \epsRate, t, ...) \f]

In course of numerical time integration, this relation will be formulated incrementally as

\f[  \displaystyle \Delta \sig = f ( \sig_n, \Delta\eps, \Delta t, t_n, ...) \f]

with

\f[  \displaystyle \Delta\eps =  \epsRate\, \Delta t \f]

and the algorithmic tangent

\f[ \displaystyle \frac{d \sig }{d \eps } =  \frac{d \Delta \sig }{d \Delta \eps } \f]

This formulation is compatible with an Abaqus interface.



\page hyperelastic Hyperelastic Material Models

## Basic Theory

**Implementation:** \ref MarmotMaterialHyperElastic

Derived abstract base class for _simple_, purely hyperelastic materials to be used for finite elements based on the total lagrangian kinematic description (TL elements). The second Piola - Kirchhoff stress tensor \f$ S \f$ will be derived by

\f[ \displaystyle S = \frac{\partial f(\boldsymbol{E},t )}{\partial \boldsymbol{E}} \f]

with the Green - Lagrange strain tensor \f$ \boldsymbol{E} \f$

\f[
  \displaystyle E  = \frac{1}{2}\,\left(\boldsymbol{F}^T\cdot \boldsymbol{F} - \boldsymbol{I} \right)
\f]

as work conjugated measure and the variable \f$ \boldsymbol{F} \f$ denoting the deformation gradient.
The algorithmic tangent will be calculated by

\f[
  \displaystyle \frac{d \boldsymbol{S}}{d \boldsymbol{E}}
\f]

\page gradmechanicalmaterials Gradient Enhanced Mechanical Material Models


**Implementation:** \ref MarmotMaterialGradientEnhancedMechanical

Base class for mechanical materials with gradient enhanced regularization to assure mesh indepency in finite element simulations.

\page gradhypoelastic Gradient Enhanced Hypoelastic Material Models

## Basic Theory
**Implementation:** \ref MarmotMaterialGradientEnhancedHypoElastic

\page substepper Substepping Algorithms

## Adaptive Substepper

### Substepper for Implicit Return Mapping Algorithms

**Implementation:** \ref AdaptiveSubstepper.h

Adaptive Substepper, employing an error estimation and Richardson extrapolation for an implicit return mapping algorithm.

### Substepper for Semi - Explicit Return Mapping Algorithms

**Implementation:** \ref AdaptiveSubstepperExplicit.h

Adaptive Substepper, employing an error estimation and Richardson extrapolation for an semi - explicit return mapping algorithm.

## Perez Fouget Substepper


### Substepper for Implicit Return Mapping Algorithms

**Implementation:** \ref PerezFougetSubstepperMarkII.h

Perez Fouget Substepper for elastoplastic materials with an implicit return mapping algorithm.

### Substepper for Semi - Explicit Return Mapping Algorithms

**Implementation:** \ref PerezFougetSubstepperExplicitMarkII.h

Perez Fouget Substepper for semi - explicit elastoplastic materials.

### Modified Version to Consider Time-Variant Elastic Stiffness Tensor

**Implementation:** \ref PerezFougetSubstepperTime.h


\page hugheswinget Hughes Winget

**Implementation:** \ref HughesWinget.h



\page voigtnotation Voigt Notation

In %Marmot, Voigt notation is used to simplify tensorial operations including second order tensors like stress and strain and fourth order tensors like stiffness and compliance.
For instance, tensor products are written as vector matrix products.

**Namespace:** \ref Marmot::ContinuumMechanics::VoigtNotation

**Implementation:** \ref MarmotVoigt.h

### Cauchy stress tensor

\f[\displaystyle
    \sig =
    \begin{bmatrix}
        \sigma_{11} & \sigma_{12} & \sigma_{13}\\
        \sigma_{21} & \sigma_{22} & \sigma_{23}\\
        \sigma_{31} & \sigma_{32} & \sigma_{33}
    \end{bmatrix}
    \rightarrow
    \begin{bmatrix}
        \sigma_{11}\\
        \sigma_{22}\\
        \sigma_{33}\\
        \sigma_{12}\\
        \sigma_{13}\\
        \sigma_{23}\\
    \end{bmatrix}
\f]

### Linearized strain tensor

\f[\displaystyle
    \eps =
    \begin{bmatrix}
        \varepsilon_{11} & \varepsilon_{12} & \varepsilon_{13}\\
        \varepsilon_{21} & \varepsilon_{22} & \varepsilon_{23}\\
        \varepsilon_{31} & \varepsilon_{32} & \varepsilon_{33}
    \end{bmatrix}
    \rightarrow
    \begin{bmatrix}
        \varepsilon_{11}\\
        \varepsilon_{22}\\
        \varepsilon_{33}\\
        2\,\varepsilon_{12}\\
        2\,\varepsilon_{13}\\
        2\,\varepsilon_{23}\\
    \end{bmatrix}
\f]


### Generalized Hooke's law

\f[
    \sig = \Cel : \eps = \mathbb{D}^{-1} : \eps \rightarrow
    \begin{bmatrix}
        \sigma_{11}\\
        \sigma_{22}\\
        \sigma_{33}\\
        \sigma_{12}\\
        \sigma_{13}\\
        \sigma_{23}\\
    \end{bmatrix}
    =
    \begin{bmatrix}
    D_{1111} & D_{1122} & D_{1133} & 2\,D_{1112} & xD_{1113} & D_{1123} \\
    D_{2211} & D_{2222} & D_{2233} & 2\,D_{2212} & D_{2213} & D_{2223} \\
    D_{3311} & D_{3322} & D_{3333} & 2\,D_{3312} & D_{3313} & D_{3323} \\
    2\,D_{1211} & 2\,D_{1222} & 2\,D_{1233} & 4\,D_{1212} & 2\,D_{1213} & 2\,D_{1223} \\
    2\,D_{1311} & 2\,D_{1322} & 2\,D_{1333} & 2\,D_{1312} & 4\,D_{1313} & 2\,D_{1323} \\
    2\,D_{2311} & 2\,D_{2322} & 2\,D_{2333} & 2\,D_{2312} & 2\,D_{2313} & 4\,D_{2323} \\
    \end{bmatrix}^{-1}
    \begin{bmatrix}
        \varepsilon_{11}\\
        \varepsilon_{22}\\
        \varepsilon_{33}\\
        2\,\varepsilon_{12}\\
        2\,\varepsilon_{13}\\
        2\,\varepsilon_{23}\\
    \end{bmatrix}
\f]



