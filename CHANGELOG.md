# Changelog

All notable changes to Marmot are documented in this file. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Changes that break existing input or code are marked
**Breaking**; changes that alter numerical results are marked **Results**.

## [Unreleased]

### Added
- The meshfree layer in-tree, imported from the BOKU-CMP-MEC-MAT repositories with their full history
  (`git subtree`): `core/MarmotMeshfreeCore` (reproducing kernel approximations, kernel functions, cells, particle
  domains, factories), `materialpoints/DisplacementMaterialPoint`,
  `materialpoints/GradientEnhancedFiniteStrainMaterialPoint`, `particles/DisplacementParticle`,
  `particles/GradientEnhancedFiniteStrainParticle` and `cells/GradientEnhancedFiniteStrainCell`. The separate
  module repositories are superseded.
- `cells/DisplacementCell`, the MPM cell for `DisplacementMaterialPoint` (Lagrangian and B-spline geometries).
- `BODYFORCE` for all displacement and gradient-enhanced finite-strain particles (it was advertised but empty, or
  threw for the SDI particles).
- `CWFCORRECTION`: `load[0]` selects the corrected traction components as a bit mask (1 = x, 2 = y, 4 = z;
  0 or no load = all), for faces where only the normal displacement is prescribed.
- Tests for all registered meshfree particles, cells and material points, and for the meshfree core.

### Changed
- **Breaking:** a rank-deficient moment matrix of the reproducing kernel approximations now throws
  (`factorizeMomentMatrix`) instead of silently returning shape functions without partition of unity. This
  includes every 3D RKPM model with a single layer of particles through the thickness (all kernel centres
  coplanar): use at least two layers of kernels, larger supports, or a lower completeness order.
- **Breaking:** `DisplacementMaterialPoint::setInitialCondition` throws for every condition, including
  `geostaticstress`, which was silently ignored.
- **Results:** the gradients of reproducing kernel shape functions of completeness order ≥ 2
  (`computeMonomialBasisGradient`), the VCI basis in 3D, the weak-form correction `CWFCORRECTION` (boundary term in
  the intermediate configuration, `1/J_Y`), and the nonlocal tangent `dL/dN` of the gradient-enhanced
  finite-strain particles.
- **Results:** the subdomains of the SDI particles are deformed about the centroid of the parent particle, so that
  they tile the deformed particle; the pressure load of `DisplacementParticleSQCNIxSDI` evaluates its test
  function at the face centre of the geometry, like the load vector and the VCI boundary term.
- The displacement NSNI particles are named like all other particles: `DisplacementSQCNIxNSNI`,
  `DisplacementSNNIxNSNI`, `DisplacementSQCNI_RxNSNI` and `DisplacementSQCNI_RUxNSNI` (2D and 3D). The former names
  `Displacement/SQCNIxNSNI`, `Displacement/SNNIxNSNI`, `Displacement/R-SNNIxNSNI` and `Displacement/RS-SNNIxNSNI`
  remain registered as aliases.
- Point particles have no faces and no longer advertise `PRESSURE` / `CWFCORRECTION`.
- The cell, particle and material point registries are function-local statics; the public `register*` functions
  are unchanged.
- The reference wave speed of the artificial bulk viscosity is cached also when it is zero (it was asked for again on
  every increment), in `DisplacementFiniteElement` and `GeneralGradientEnhancedDisplacementFiniteElement`. A material
  meant to run with the damage degradation should override `getMaximumWaveSpeed()` with its closed form; the new
  `Marmot::Testing::checkMaximumWaveSpeedOverrideMatchesDefault()` (`MarmotWaveSpeedContract.h`) checks an override
  against the default.
- `GeneralGradientEnhancedDisplacementFiniteElement` no longer stores the strain-displacement operator per quadrature
  point: `QuadraturePoint::B` and the unused `B_K` are removed (code reading them must be rebuilt and use
  `localGeometryElement.B( qp.dNdX )`), and `B` is built from `dNdX` where the kernels use it. 16.7 of 20.7 GB on a
  107,508-element `GC3D20` mesh; results are bitwise identical.
- **Results:** the kernels of `GeneralGradientEnhancedDisplacementFiniteElement` (`GC3D8`, `GC3D20`, `GCPS4`, ... and
  the mixed-order `M` variants) compute the strain increment and the internal force from `dNdX` through
  `FiniteElement::SmallStrain::strainIncrement` and `addInternalForce` (new, in `MarmotSmallStrainKinematics.h`)
  instead of the dense B-matrix, and the gradient term of the non-local residual and tangent in O(nNodes) per
  quadrature point. Residuals and tangents change in the last digits (about 1e-16 relative to the residual scale);
  the explicit edge-breakout run is 25 % faster. The quadrature points are smaller (about 5 KB to under 1 KB per
  point for `GC3D20R`); equal-order elements store the non-local shape functions once.

### Fixed
- **Results:** sign of the hoop term in the axisymmetric finite-strain UL elements (`CX8UL`, `CX8RUL`, residual and
  stiffness): it was subtracted, so any state with nonzero hoop stress was wrong (e.g. homogeneous radial expansion
  off by a factor of roughly 5–60). Uniaxial loading with zero hoop stress was unaffected.
- Uninitialized intermediate second moments of the SQCNIxNSNI particles in the first increment.
- The inverse isoparametric map of distorted Lagrangian cells (it threw) and their point location test.
- Stack overflows of the 64-node B-spline cells on Windows (1 MB default stack).
