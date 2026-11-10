# Changelog

All notable changes to Marmot are documented in this file. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Changes that break existing input or code are marked
**Breaking**; changes that alter numerical results are marked **Results**.

## [Unreleased]

### Added
- Linear simplices in `elements/DisplacementFiniteElement`: the triangles `CPS3` (plane stress) and `CPE3` (plane
  strain), and the tetrahedron `C3D4`, each with one Gauss point. The strain is constant in a linear simplex, so one
  point integrates the stiffness and the internal force exactly. Their center, characteristic element length and
  regularisation length come from `GeometryElement` (below).
- `Tri3` shape: shape functions, the 1-point (reduced) and 3-point (full, exact for quadratics) Gauss rules, the
  `GeometryElement<2, 3>` specialisations and the EnSight shape `tria3`.
- Distributed loads (pressure, surface traction) on `Tetra4` and `Tri3` elements (`C3D4`, `CPS3`, `CPE3`). This adds
  the face definitions `Tetra4::getBoundaryElementIndices` (Abaqus faces S1–S4; previously a stub that always threw)
  and `Tri3::getBoundaryElementIndices`, and `BoundaryElement` support for both parent shapes.
- `getReferenceDomain( shape )`, `getInterpolationOrder( shape )` (the order p of the shape functions) and
  `getNumberOfCornerNodes( domain )`, all `constexpr` like `getElementShapeByMetric`, which moved into the header
  and now also knows `Bar3` (1D, 3 nodes). Further facts of a `ReferenceDomain`: `isSimplex( domain )`,
  `getReferenceMeasure( domain )` (2^nDim for the cube [-1, 1]^nDim, 1 / nDim! for the unit simplex) and
  `getReferenceCentroidCoordinate( domain )`.
- Element geometry that depends on the reference domain, in `GeometryElement`, so the elements do not have to
  distinguish the cube from the simplex: the compile-time `shape` (previously a member initialised at construction),
  `referenceDomain` and `interpolationOrder`,
  - `center()`: the image of the centroid of the reference domain; the origin of the unit simplex is its node 0;
  - `characteristicLengthForCriticalTimeStepAt( xi )`: the length L_e of the stable time increment L_e / c, twice the
    smallest singular value of the Jacobian for quadrilaterals and hexahedra (their smallest extent), 1 / sqrt( sum_i
    |grad N_i|^2 ) for linear simplices. The latter does not depend on the node numbering and stays below the exact
    stable increment of the element; the singular values assume the reference cube and would exceed it by up to 4
    times;
  - `volumeEquivalentLength( detJ )`: the regularisation length for softening materials, sqrt( A ) or cbrt( V ).

  `characteristicLengthForCriticalTimeStepAt` is not derived for quadratic simplices (`Tri6`, `Tetra10`); asking for
  it is a compile error.
- Nodal integration with closed Newton-Cotes rules (`Quadrature::newtonCotesExactToDegree( degree )`) on `Line`, `Quadrilateral`
  and `Triangle`: the rules of degree 1 have the corners of the domain as points, i.e. one point at each node of the
  linear element (`Bar2`, `Quad4`, `Tri3`) in node order, with the integral of its shape function as weight.
- Gauss rules of degree 5 on simplices: the 7-point triangle rule (Radon 1948) and the 14-point tetrahedron rule
  (Walkington 2000), both with positive weights and interior points. The 14-point rule makes the consistent mass of
  a quadratic tetrahedron (`gaussExactToDegree( 4 )`) available.

### Changed
- **Breaking:** `PRESSURE` on the faceted particles (`DisplacementParticleSQCNI`, implicit and explicit,
  `DisplacementParticleSQCNIxSDI`, `GradientEnhancedFiniteStrainParticleSQCNI` and the NSNI particles derived from
  them) has the sign of the finite elements: a positive pressure pushes on the face, against its outward normal. It
  pulled before. Inputs keep their results by flipping the sign of the pressure.
- `DisplacementFiniteElement`, `GeneralGradientEnhancedDisplacementFiniteElement`,
  `DisplacementFiniteStrainULElement` and `GradientEnhancedFiniteStrainDisplacementElement` take their center,
  characteristic length for the critical time step and (the first) regularisation length from `GeometryElement`. The
  results of the quadrilaterals, hexahedra and bars are unchanged; instantiated for linear simplices (not registered
  yet, except `CPS3`, `CPE3` and `C3D4`), the gradient-enhanced and the finite-strain elements now also have the right
  center and length. The element methods `characteristicElementLengthAt( xi )` of `DisplacementFiniteElement` and
  `GeneralGradientEnhancedDisplacementFiniteElement` are replaced by
  `GeometryElement::characteristicLengthForCriticalTimeStepAt( xi )` (for the latter through its member
  `localGeometryElement`), so the name no longer collides with the material's `characteristicElementLength`, the
  regularisation length.
- **Breaking:** quadrature rules are requested by family and degree of exactness instead of `FullIntegration` /
  `ReducedIntegration`, whose meaning differed from shape to shape. All rules are kept in one catalogue,
  `Quadrature::quadratureRules()`, each entry with its reference domain (`Line`, `Triangle`, `Quadrilateral`,
  `Tetrahedron`, `Hexahedron`), its family (`Gauss`, `NewtonCotes`) and the degree it integrates exactly.
  `Quadrature::getQuadraturePoints( domain, rule )` returns the rule with the fewest points on the reference domain
  that reaches the requested degree, and throws if there is none; `getQuadraturePoints( shape, rule )` does the same
  on the reference domain of an element shape. Elements take a
  `Quadrature::QuadratureRule`, built by `gaussExactToDegree( degree )` or `newtonCotesExactToDegree( degree )`, instead of `IntegrationTypes`;
  the registrations state the degree their stiffness needs (e.g. `CPS4` `gaussExactToDegree( 2 )`, `CPS8R` `gaussExactToDegree( 3 )`). The
  consistent mass asks for `gaussExactToDegree( 2p )` (N^T N) and the faces of distributed loads for `gaussExactToDegree( 3p - 1 )` of the face
  shape (N times the area vector of the face). Every registered element keeps its points, in the same order, so
  results are unchanged. Migration:

  | Old | New |
  |---|---|
  | `IntegrationTypes::FullIntegration` | `gaussExactToDegree( 2p )`; for `Tetra4` `gaussExactToDegree( 1 )`, for `Tetra10` `gaussExactToDegree( 2 )` |
  | `IntegrationTypes::ReducedIntegration` | `gaussExactToDegree( 2p - 1 )`; for `Tetra4` `gaussExactToDegree( 1 )`, for `Tetra10` `gaussExactToDegree( 2 )` |
  | `getGaussPointInfo( shape, type )` | `getQuadraturePoints( shape, rule )` |
  | `getNumGaussPoints( shape, type )` | `getQuadraturePoints( shape, rule ).size()` |
  | `Spatial1D::gaussPointList1/2/3` | `getQuadraturePoints( ReferenceDomain::Line, gaussExactToDegree( 1/3/5 ) )` |
  | `Spatial2D::gaussPointList1x1/2x2/3x3` | `getQuadraturePoints( ReferenceDomain::Quadrilateral, gaussExactToDegree( 1/3/5 ) )` |
  | `Spatial3D::gaussPointList1x1x1/2x2x2/3x3x3` | `getQuadraturePoints( ReferenceDomain::Hexahedron, gaussExactToDegree( 1/3/5 ) )` |
  | `Spatial3D::gaussPointListTetra4/Tetra10` | `getQuadraturePoints( ReferenceDomain::Tetrahedron, gaussExactToDegree( 1/2 ) )` |

### Fixed
- **Results:** elastic strain energy of `LinearElastic`, `ADLinearElastic` and every material on
  `MaterialHypoElasticAD`. The fixes:
  - `LinearElastic` did not compute the elastic energy density at all, so the internal energy of the elements was
    zero.
  - `ADLinearElastic` accumulated ½ Δε : Δσ without the stress at the start of the increment, so it was too small
    by about the number of increments.
  - `MaterialHypoElasticAD::computeStress` returned only the stress of `computeStressAD` and dropped its energy and
    dissipation.
  - Both linear elastic laws now integrate W by the trapezoidal rule, which is exact for linear elasticity.
  - Tests compare W with ½ ε : C : ε after a non-proportional strain path in 3D and in plane stress.
- Result-file keyword of 2D interface elements (`ILINE2`): `GeometryInterfaceElement` named the four-node
  element `bar2`, a two-node cell, so the EnSight geometry was invalid; it is now `quad4`, the cell spanning both
  faces, as `hexa8` for the 3D interfaces.
- Lumped inertia and critical time step of the linear simplices `CPS3`, `CPE3` and `C3D4`. The manifold-based mass
  lumping blended in the shape functions of a quadrilateral or hexahedron and wrote past the end of the nodal
  vectors, which corrupted the heap of an explicit analysis. A linear simplex is now its own linear counterpart
  (`MassLumping::nLinearNodes()`, also used by the finite-strain and the gradient-enhanced elements), with the
  mass / nNodes at each node. `MassLumping::nLinearNodes()` is the number of corner nodes of the reference domain,
  which also corrects `Tetra10` (4 instead of 8; not registered).
- Consistent mass of `Tetra4` elements (`C3D4`). Its Gauss rule has one point for every integration type, which gave
  a mass of rank 3 out of 12. `ConsistentMass::integrationRule()` now takes the 4-point rule (`gaussExactToDegree( 2 )`), which is
  exact for the quadratic integrand.
- **Results:** `newmarkBetaIntegration` (`MarmotMechanicsCore`, used by all particles with the Newmark properties):
  a zero-length increment adds no inertia; velocity and acceleration are left as they are and the inertia tangent
  is zero. The time increment used to be clamped to 1e-16, which made a zero-length increment with a nonzero
  velocity produce an acceleration of order -v/(beta 1e-16). Increments of at least 1e-16 are unchanged.
- **Results:** `newmarkBetaIntegration` with `beta = 0` (no time integration, a quasi-static analysis) sets the velocity
  and the acceleration to zero: a quasi-static step after a dynamic one leaves the particles at rest. It used to keep
  the velocity, advanced once by `dT (1 - gamma) a_n`. Analyses that are quasi-static from the start are unchanged.
- **Results:** precision of the 3-point Gauss rules. The abscissa sqrt(3/5) was tabulated 4 ulp off and the 3x3x3
  weights up to 12 ulp off (their sum missed 8 by 7e-15); both are now correctly rounded, the weights as exact fractions
  over 729. The abscissas are no longer public (`Quadrature::gp2`, `Quadrature::gp3`); like all rule data they live in
  `quadratureRules()`. This changes the last digits (about 1e-15 relative) of all elements and faces that use the
  3-point, 3x3 or 3x3x3 rule: the full-integration Quad8 and Hexa20 elements, pressure on Quad8 and Hexa20 faces, and
  the consistent mass of Quad8 and Hexa20.

## [v26.11] - 2026-11-02

This changelog was started late in the v26.11 cycle. The entries below cover October 2026 onwards; the main
changes from earlier in the cycle, since v25.11, were:

- Refactored material interfaces (#25) and a new element internal-force sign convention without `pNewDT` and
  `stepTime` (#50). **Breaking**
- Finite-strain materials: a material point solver (#30), automatic-differentiation materials (#36), a substepper
  (#32), linear viscoelasticity (#35), the Hughes-Winget co-rotational wrapper (#64) and the gradient-enhanced
  finite-strain interface (#98).
- Explicit dynamics: artificial bulk viscosity, nonlocal micro-inertia, mass lumping and the critical time step (#45,
  #81, #84).
- The general gradient-enhanced displacement element with multiple nonlocal fields, and the AT2 phase-field model.
- Interface elements and materials (#40).
- A Python interface to the material point solvers (#53), an element named-property interface (#54) and custom
  material and solver exceptions (#43).
- Builds: one CMake target per module with package export (#47, #113, #114), Windows DLLs with MSVC (#110),
  macOS/arm64 (#57), Eigen 5 (#107), free-threading-safe statics (#55), a conda package (#121), and a coverage
  pipeline (#83).
- Fixes that change results: `projectVoigtStressToPlane` shear column (#90), the Hughes-Winget rotation (#89), the
  consistent mass integration (#101) and Hexa20 lumped corner masses (#74). **Results**

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
- `Marmot::Microplane` in `MarmotMechanicsCore`: the microplane quadrature schemes (N21 symmetric/asymmetric, N28,
  N37, N61) and the volumetric-deviatoric (`VDSplit`) and normal-tangential (`NTSplit`) projection tensors with
  their homogenization helpers, shared by all microplane models.
- `materials/MicroplaneM7` (`MICROPLANEM7`), the microplane model M7 for concrete (Caner and Bažant, 2013), and
  `materials/GM7` (`GM7`), its gradient-enhanced nonlocal extension. The former standalone module repositories are
  superseded.

### Changed
- **Breaking:** complete namespace rework. All code now lives in the global `Marmot` namespace with a closed set
  of second-level namespaces (`Utilities`, `Math`, `ContinuumMechanics`, `Materials`, `Discretizations`,
  `Registry`, `Testing`); classes dropped the redundant `Marmot` name prefix; the top-level `MarmotLibrary`
  namespace is gone. There are **no compatibility aliases** — consumers migrate with this table.
  Header file names are unchanged, with three exceptions: `Marmot/MarmotUtility.h` (empty, deleted — include
  `Marmot/MarmotTypedefs.h` instead) and `Marmot/VonMisesConstants.h` / `Marmot/ADVonMisesConstants.h`
  (dissolved into the model headers, see table):

  | Old | New |
  |---|---|
  | `MarmotLibrary::Marmot<X>Factory` | `Marmot::Registry::<X>Factory` (e.g. `ElementFactory`, `MaterialHypoElasticFactory`, `MaterialFiniteStrainFactory`, `MaterialGradientEnhancedFiniteStrainFactory`, `MaterialGradientEnhancedHypoElasticFactory`, `ParticleFactory`, `MaterialPointFactory`, `CellFactory`, `CellElementFactory`) |
  | `::MarmotMaterialHypoElastic(/AD)`, `::MarmotInterfaceMaterialHypoElastic` | `Marmot::Materials::MaterialHypoElastic(/AD)`, `::InterfaceMaterialHypoElastic` |
  | `::MarmotMaterialFiniteStrain(/AD)`, `::MarmotMaterialGradientEnhancedFiniteStrain` | `Marmot::Materials::MaterialFiniteStrain(/AD)`, `::MaterialGradientEnhancedFiniteStrain` |
  | `::MarmotMaterialGeneralGradientEnhancedHypoElastic<>` | `Marmot::Materials::MaterialGradientEnhancedHypoElastic<>` |
  | `::MarmotElement`, `::MarmotGeometryElement<>`, `::MarmotElementSpatialWrapper` | `Marmot::Discretizations::Elements::{Element, GeometryElement<>, ElementSpatialWrapper}` |
  | `::MarmotMaterialSection`, `::ElementProperties` | `Marmot::Discretizations::{MaterialSection, ElementProperties}` |
  | `::MarmotMaterialPoint`, `::MarmotCell`, `::MarmotCellElement` | `Marmot::Discretizations::{MaterialPoints::MaterialPoint, Cells::Cell, Cells::CellElement}` |
  | `::MarmotJournal`, `::MakeString`, `::StateView`, `::MarmotStateVarVectorManager` | `Marmot::Utilities::{Journal, MakeString, StateView, StateVarVectorManager}` |
  | `Marmot::Elements` / `Marmot::Meshfree` / `Marmot::MaterialPoints` / `Marmot::Cells` | `Marmot::Discretizations::{Elements, Meshfree, MaterialPoints, Cells}` (meshfree classes dropped prefixes: `Particle`, `Approximation`, `KernelFunction*`) |
  | `Marmot::FiniteElement::*` (shape functions, quadrature, EAS, …) | `Marmot::Discretizations::Elements::*` (`Spatial2D::axisymmetric` → `Axisymmetric`) |
  | `Marmot::Constants`, `Marmot::EigenTensors` + `ContinuumMechanics::{CommonTensors, TensorUtility}` | `Marmot::Math::{Constants, Tensors}` |
  | `Marmot::FastorStandardTensors`, `Marmot::FastorIndices`, `ContinuumMechanics::TensorUtility::TensorExponential`, plain-`Marmot` Fastor helpers (`mapEigenToFastor`, `makeDual`/`makeReal`, `deviatoric`, `expandTo3D`/`reduceTo2D`, …) | `Marmot::Math::Tensors::{FastorTensors, FastorIndices, TensorExponential}` and functions directly in `Marmot::Math::Tensors` |
  | `Marmot::AutomaticDifferentiation`, `Marmot::NumericalAlgorithms::Differentiation` | `Marmot::Math::Differentiation::{Automatic, Numerical}` |
  | `Marmot::NumericalAlgorithms::Integration` | `Marmot::Math::Integration` |
  | `Marmot::NumericalAlgorithms` (substeppers, `HughesWinget`, `YieldSurfaceCombinationManager`) | `Marmot::ContinuumMechanics::Algorithms` |
  | `Marmot::TimeIntegration`, `Marmot::GeostaticStress` | under `Marmot::ContinuumMechanics` |
  | `Marmot::GradientDamage`, `Marmot::PhaseField` | `Marmot::ContinuumMechanics::Degradation::{GradientDamage, PhaseField}` |
  | `Marmot::Solvers::MarmotMaterialPointSolver{HypoElastic, FiniteStrain}` | `Marmot::Materials::MaterialPointSolver{HypoElastic, FiniteStrain}` |
  | `ContinuumMechanics::FiniteStrain::*` (namespace dissolved) | `ContinuumMechanics::Kinematics::DeformationMeasures`, `::Elasticity::EnergyDensityFunctions`, `::StressMeasures`, `::Plasticity`, `::Viscoelasticity` (merged with the small-strain one) |
  | `ContinuumMechanics::CommonConstitutiveModels` (namespace dissolved) | `ContinuumMechanics::Plasticity::MenetreyWillam`, `ContinuumMechanics::Viscoplasticity::DuvautLionsViscosity` |
  | `Marmot::Materials::{KelvinChain, Wiechert, SolidificationTheory}` | `Marmot::ContinuumMechanics::Viscoelasticity::{KelvinChain, Wiechert, SolidificationTheory}` |
  | `ContinuumMechanics::VoigtNotation` (conversions, `Transformations`) | `ContinuumMechanics::Voigt` |
  | `ContinuumMechanics::VoigtNotation::{Invariants, Derivatives}`, `ContinuumMechanics::HaighWestergaard` | `ContinuumMechanics::Invariants`, `::Invariants::Derivatives`, `::Invariants::HaighWestergaard` |
  | derivative namespaces `*OrderDerived` spellings | unified `FirstOrder` / `SecondOrder` / `ThirdOrder` |
  | `Marmot::Materials::Shrinkage::B4` | `Marmot::Materials::B4Shrinkage` |
  | `Marmot::Materials::VonMisesConstants` / `::ADVonMisesConstants` (headers `VonMisesConstants.h`, `ADVonMisesConstants.h`) | class members `VonMisesModel::Constants` / `ADVonMises::Constants` (in `VonMises.h`, `ADVonMises.h`) |
  | `Marmot::Microplane` (quadrature schemes, `VDSplit`, `NTSplit`) | `Marmot::ContinuumMechanics::Microplane` |
  | plain-`Marmot` typedefs (`Vector6d`, `Matrix6d`, `mVector6d`, …), general math functions and B-spline/monomial bases | `Marmot::Math::*` (root) |

  Factory registration **string keys are unchanged**, so name-based material/element creation through the
  factories behaves identically. Registration translation units now use unnamed namespaces. The convention and
  its rules are documented in `AGENTS.md`; `scripts/checkNamespaceConventions.sh` checks them.
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
  `Discretizations::Elements::SmallStrain::strainIncrement` and `addInternalForce` (new, in
  `MarmotSmallStrainKinematics.h`)
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

[Unreleased]: https://github.com/MAteRialMOdelingToolbox/Marmot/compare/v26.11...HEAD
[v26.11]: https://github.com/MAteRialMOdelingToolbox/Marmot/compare/v25.11...v26.11
