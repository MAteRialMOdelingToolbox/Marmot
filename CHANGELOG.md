# Changelog

All notable changes to Marmot are documented in this file. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Changes that break existing input or code are marked
**Breaking**; changes that alter numerical results are marked **Results**.

## [Unreleased]

### Added
- `elements/TrussElement`: truss (bar) elements `T2D2`, `T2D3`, `T3D2`, `T3D3` (geometrically linear, hypoelastic
  materials in uniaxial stress) and `T2D2CR`, `T2D3CR`, `T3D2CR`, `T3D3CR` (total Lagrangian, finite strain
  materials evaluated in the frame co-rotated with the bar's tangent: rigid rotations are stress free, consistent
  material and geometric tangent; for 2 nodes identical to the co-rotational truss). Property: the cross section
  area. E.g., for reinforcement bars embedded in continuum elements.
- `Materials::MaterialFiniteStrain::computeUniaxialStress` (virtual): uniaxial stress reduction of finite strain
  materials whose response to a diagonal deformation gradient is diagonal, e.g. isotropic ones (Newton iteration on
  the lateral stretches with state reset, condensed tangent dτ11/dF11).
- Bond-slip of embedded reinforcement: the pluggable bond-slip law interface `Materials::BondSlipLaw` with its
  `Registry::BondSlipLawFactory` (bond laws are modules in the materials category with their own factory), the laws
  `materials/LinearElasticBondSlip` (`LINEARELASTICBONDSLIP`), `materials/ModelCode2010BondSlip`
  (`MODELCODE2010BONDSLIP`: fib Model Code 2010 envelope, elastic un-/reloading, optional Coulomb friction and a
  normal crushing cap), and `elements/EmbeddedBondElement` (`EB<nDim>D<nBarNodes><host>`, e.g. `EB2D2Q4`,
  `EB3D2H20`, `EB2D2T3`: the bond of a bar element to the host element it lies in, small slip, geometrically
  linear).
- `GeometryElement::parentCoordinates` (inverse isoparametric mapping on coordinates relative to the element, so it
  converges far from the origin) and `GeometryElement::isInsideParentDomain`; `LagrangianCellGeometry` and the
  embedded bond element use it.
- The Gauss-Lobatto family of the quadrature catalogue on Line (2 to 5 points), `gaussLobattoExactToDegree()`.
- Large slip of embedded reinforcement: `EmbeddedLargeSlipBondElement` (`EBLS<nDim>D<nBarNodes><host>W<n>`, e.g.
  `EBLS2D2Q4W3`): the bond is integrated over the channel fixed in the host element and each channel point searches
  its current partner (closest point projection) in a window of `n` bar elements, so the bar can slide out of its
  host and the bonded length shrinks as it is pulled out. Objective (current frame at the partner), with the
  linearization of the closest point; `EmbeddedBondKinematics.h` holds the frame and the integration rule shared with
  `EmbeddedBondElement`; `Bar2::d2NdXi2`, `Bar3::d2NdXi2`.
- Energy gradients of the elements for dissipation based arc-length control, as optional capabilities that leave the
  standard interfaces unchanged (see the feature page "Energy gradients"):
  - `EnergyGradientProvider` (elements): `assignEnergyGradientOutput()` assigns two arrays to which every
    `computeKernels()` call writes the gradients of the element's elastic energy and dissipation with respect to its
    degrees of freedom; unsupported sections or materials are rejected when the outputs are assigned.
    `DisplacementFiniteElement` (solid, plane strain) and `GeneralGradientEnhancedDisplacementFiniteElement` (solid,
    plane strain, plane stress, with the non-local degrees of freedom) provide it.
  - `HypoElasticEnergyTangentProvider` and `GradientEnhancedHypoElasticEnergyTangentProvider` (materials):
    `computeStressWithEnergyTangents()`, the stress update with the derivatives of the updated elastic energy density
    and dissipation, and `computePlaneStressWithEnergyTangents()` for the gradient-enhanced plane stress iteration.
  - Energy tangents of `VonMises`, `RankineSofteningPlasticity`, `RankineDamagePlasticity`,
    `RankineGradientSofteningPlasticity` and `RankineGradientDamagePlasticity`.
  - `MaterialGradientEnhancedHypoElastic::iteratePlaneStress()`: the plane stress iteration of `computePlaneStress()`
    around a given 3D stress update.
- `elements/BeamElement`: Euler-Bernoulli beams `B2D2`, `B2D3`, `B3D2`, `B3D3` (2 or 3 nodes, geometrically
  linear; cubic / quintic Hermite transverse interpolation), node fields displacement and rotation, sections given
  as section points (any shape, `[n, y_i, A_i ...]` in 2D, `[v, J, n, y_i, z_i, A_i ...]` in 3D) with one material
  instance per point (any hypoelastic material), Saint-Venant torsion scaled to G J. `BeamFiberSection` integrates
  the section independently of the kinematics; its 2D fibers are in uniaxial stress (`computeUniaxialStress`), its 3D
  fibers in the beam stress state σ22 = σ33 = σ23 = 0 (`ReducedStress::computeReducedStress` with its index sets
  `BeamFiberSection::beamStressPrescribed = {0, 3, 4}`, `beamStressFree = {1, 2, 5}`). A body force acts at the
  section points: line load and, for a section whose centroid is off the axis, its line moment.
- `B2D2CR` and `B3D2CR`: co-rotational (finite rotation) variants of `B2D2` / `B3D2` with the same properties
  and section (2D: Crisfield; 3D: Battini and Pacoste, exact geometric stiffness by forward automatic
  differentiation); arbitrary rigid rotations, small local strains.
- Quadrature catalogue: the 4-point Gauss-Legendre rule on `Line` (degree 7).
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
- `materials/RankineSofteningPlasticity` (`RANKINESOFTENINGPLASTICITY`), multi-surface Rankine plasticity with the
  isotropic hardening/softening law of VonMises, the yield stress floored at zero;
  `materials/RankineDamagePlasticity` (`RANKINEDAMAGEPLASTICITY`), perfect Rankine plasticity in the effective stress
  with exponential damage driven by the plastic strain (Grassl and Jirásek, 2006), bounded by a damage threshold; and `materials/RankineGradientDamagePlasticity` (`RANKINEGRADIENTDAMAGEPLASTICITY`), its
  implicit-gradient (over-)nonlocal extension. The hypoelastic models are also registered with the Hughes-Winget
  wrapper (`.../HUGHES-WINGET`). Imported with their history (`git subtree`); the module repositories are superseded.
- `materials/RankineGradientSofteningPlasticity` (`RANKINEGRADIENTSOFTENINGPLASTICITY`), the implicit-gradient
  over-nonlocal extension of `RankineSofteningPlasticity` (Vermeer and Brinkgreve, 1994): the yield stress depends on
  a combination of the local and the nonlocal softening variable.
- Small-strain axisymmetric elements in `elements/DisplacementFiniteElement`: `CAX4` (4 nodes, full integration),
  `CAX8` (8 nodes, full integration) and `CAX8R` (8 nodes, reduced integration), with r = x_1 and z = x_2. They are
  `DisplacementFiniteElement< 2, nNodes >` with the new section `SectionType::Axisymmetric`, which adds the hoop
  strain u_r / r and integrates over the full ring (2 pi r), so internal forces, energies, masses, body forces and
  distributed loads are those of the full revolution. The artificial bulk viscosity and the energy gradients are
  available. `initializeYourself()` rejects a node at r < -1e-12 h (h the size of the element) and a
  quadrature point at r <= 0.
- `BoundaryElement( ..., axisymmetric )`: with `axisymmetric = true`, its load vectors are integrated over the ring
  (2 pi r ds).
- `Spatial2D::Axisymmetric::BFromPlanarB( BPlanar, N, r )`: the axisymmetric B-operator from the planar one, the
  single definition of the hoop row; `Spatial2D::Axisymmetric::B` (`GeometryElement::B_axisymmetric`) and the
  axisymmetric section of `DisplacementFiniteElement` use it.
- `DisplacementFiniteElement::sectionFactorAt( xi )`: the thickness, cross section, 2 pi r or 1 between the volume
  and the reference measure, used for the quadrature weights and the consistent mass.
- Small-strain axisymmetric gradient-enhanced elements in `elements/GeneralGradientEnhancedDisplacementFiniteElement`:
  `GCAX4` (4 nodes, full integration), `GCAX8` (8 nodes, full integration) and `GCAX8R` (8 nodes, reduced
  integration), the gradient-enhanced counterparts of `CAX4`, `CAX8` and `CAX8R`. They are
  `GeneralGradientEnhancedDisplacementFiniteElement< 2, nNodes >` with the new section `SectionType::Axisymmetric`:
  the hoop strain u_r / r, and the full ring (2 pi r) for the internal forces, the Helmholtz equation of the
  non-local fields, energies, masses, body forces and pressure. Implicit and explicit, with the artificial bulk
  viscosity and the energy gradients; mixed-order variants (e.g. `< 2, 8, 1, 4 >`) work as well. Nodes and quadrature
  points are checked as for the CAX elements. Based on the axisymmetric section by Thomas Mader.
- `GeneralGradientEnhancedDisplacementFiniteElement::sectionFactorAt( xi )`, as for `DisplacementFiniteElement`.
- `Spatial2D::Axisymmetric::checkNodesAreInTheMeridianHalfPlane()`, `checkRadiusOfQuadraturePoint()` and
  `sectionTangent()`: the checks of the meridian section and the 4 × 4 block of a 3D tangent, shared by the
  axisymmetric sections of both elements.
- `ContinuumMechanics::Voigt::reduce3DVoigt< Axial >`, the counterpart of `make3DVoigt< Axial >`.
- `Marmot::ContinuumMechanics::ReducedStress` (`MarmotLowerDimensionalStress.h`, header only): one static
  condensation for reduced stress states with prescribed strain components `p` and stress-free components `f`.
  `computeReducedStress( computeStress3D, dStrain, p, f )` iterates a 3D stress update (a callable) on the strains at
  `f` (Newton, tolerance relative to the largest stress of the call, 1e-12, 1e-9 from the 9th evaluation, at most 15)
  and returns the condensed tangent, the sensitivity `dStrainF_dStrainP = -Cff^+ Cfp` and the complete orthogonal
  decomposition of `Cff` (pseudo-inverse `Cff^+`, the inverse if it is regular), so that callers condense further
  derivatives the same way (`freeStrainSensitivity( condensation, dStressF_dX )`, which checks that a singular `Cff`
  can compensate them); `condense( C, p, f )` is its algebra for a given tangent. The index sets of plane stress (`planeStressPrescribed`, `planeStressFree`) and uniaxial stress
  (`uniaxialStressPrescribed`, `uniaxialStressFree`) are defined there.
- Transport (diffusion) problems: the core module `MarmotTransportCore` with the material interface
  `Materials::TransportMaterial` (content, flux, source and their tangents for one or more scalar fields; backward
  Euler in content form, so the balance is conservative) and its `Registry::TransportMaterialFactory`; the element
  family `elements/TransportFiniteElement` (`DC2D3T`, `DC2D4T`, `DC2D8T`, `DCAX4T`, `DCAX8T`, `DC3D4T`, `DC3D8T`,
  `DC3D20T`, "diffusive continuum" with the field `temperature`); the material `materials/HeatConduction`
  (`HEATCONDUCTION`: anisotropic, temperature dependent conductivity and capacity; 3, 6 or 11 properties); the
  distributed load types `HeatFlux`, `HeatFilm` (convection) and `HeatRadiation`. A zero time increment lets
  nothing flow.
- `materials/MicroplaneM7`, `materials/GM7`: optional property `stressUnit`, the value of 1 MPa in the user's stress
  unit (default 1.0, e.g. 1e6 for Pa). It converts the dimensional fixed parameters of Caner and Bažant (2013)
  (`c18`, `c21`, the reference strength and modulus), the default values of `fc` and `transitionWidth` and the
  stiffness tolerance of the unloading transition, so that the models can be used in any consistent unit system.
  Results with the default are unchanged.

### Changed
- **Breaking:** Marmot requires Eigen 5; Eigen 3.4 is no longer supported (CMake rejects Eigen < 5.0, the conda package
  requires `eigen >=5`). CI builds and tests with Eigen 5.0.1 and autodiff 1.1.2 with the Eigen 5 patch
  (autodiff#397).
- **Breaking:** `T2D2` is now the 2D 2-node element of the `TrussElement` module instead of the spatial wrapper of
  the 1D `DisplacementFiniteElement` (uniaxial stress). Same property (the cross section area), same small strain
  response for hypoelastic materials; but 1 quadrature point instead of 2 (a deck selecting `quadraturePoint=1`
  must use 0), and the quadrature point states are those of the truss element (`stress`, `strain`,
  `normal force`, ...). `ElementSpatialWrapper` and the uniaxial-stress section of `DisplacementFiniteElement` stay as
  generic facilities (tested), now without a registered element.
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

- `DisplacementFiniteElement`: `computeKernels()` and `computeKernelsExplicit()` share one implementation, and the
  solid, plane strain and axisymmetric sections one material update through the 3D interface. Results are unchanged.
  The consistent mass takes the thickness or cross section from the element property instead of recovering it from
  the nearest quadrature point (identical up to rounding).
- `GeneralGradientEnhancedDisplacementFiniteElement`: likewise, `computeKernels()` and `computeKernelsExplicit()`
  share one material update and integration per quadrature point for the plane stress, plane strain, axisymmetric
  and solid sections, which removes the copies of the non-local coupling per section. Results are unchanged. The
  consistent mass takes the thickness from the element property (identical up to rounding), so the element
  properties must stay valid after `initializeYourself()`, as for `DisplacementFiniteElement`. A one-dimensional
  instance (none is registered), for which no section is implemented, throws instead of integrating nothing.
- `MarmotGradientMechanicsCore` requires `MarmotMechanicsCore` (for `ReducedStress`); a build with a `CORE_MODULES`
  filter that selects the former must select the latter, too.
- **Breaking:** `materials/MicroplaneM7`, `materials/GM7`: the new property `stressUnit` is inserted after
  `transitionExponent` (index 11 for `MICROPLANEM7`, 13 for `GM7`), in front of the properties only needed for explicit
  simulations. The density of `MICROPLANEM7` moves from index 11 to 12; the density of `GM7` from 13 to 14 and its
  nonlocal viscosity from 14 to 15. In an old input file giving these properties, the density would be read as
  `stressUnit`.

### Deprecated
- `BoundaryElement::computeSurfaceNormalVectorialLoadVectorForAxisymmetricElements()`: construct the boundary element
  with `axisymmetric = true` and call `computeSurfaceNormalVectorialLoadVector()`.

### Removed
- **Breaking:** `ContinuumMechanics::UniaxialStress::getUniaxialStressTangent()` (its only user was
  `MaterialHypoElastic::computeUniaxialStress()`) and `ContinuumMechanics::PlaneStress::getPlaneStressTangent()` (used
  by `MaterialHypoElastic::computePlaneStress()` and `GeneralGradientEnhancedDisplacementFiniteElement`): the condensed
  tangent comes from `ReducedStress`.
- **Breaking:** `ContinuumMechanics::PlaneStress::dStrainDStrainPlaneStress( C )` and `dStressPlaneStressDStress()`,
  without a caller since `getPlaneStressTangent()` is gone; they defined plane stress by σ33 = 0 only. Plane stress is
  `ReducedStress::planeStressPrescribed` / `planeStressFree` with `ReducedStress::condense()`.
- **Breaking:** `ContinuumMechanics::PlaneStress::planeStressCompensationStrain( strain, nu )` and
  `planeStressTangentTransformationMatrix( C )`, without a caller; they were a third definition of plane stress (σ33 = 0
  only, isotropic). Plane stress is `ReducedStress::planeStressPrescribed` / `planeStressFree`.

### Fixed
- **Results, Breaking:** reduced stress states of the small strain materials. `MaterialHypoElastic::computePlaneStress()`
  and `computeUniaxialStress()` (still virtual) and `MaterialGradientEnhancedHypoElastic::iteratePlaneStress()` (the
  default `computePlaneStress()` and `computePlaneStressWithEnergyTangents()`) are thin wrappers around
  `ReducedStress::computeReducedStress()` with the full definitions of the stress states: plane stress
  σ33 = σ13 = σ23 = 0, uniaxial stress all stresses but σ11 zero. They held γ13 = γ23 = 0 (plane stress) and
  γ12 = γ13 = γ23 = 0 (uniaxial stress) instead, which is exact only for materials without coupling between normal and
  out-of-plane shear components (isotropic, orthotropic in the axes of the plane); for other materials, e.g. an
  orthotropic `LinearElastic` with axes rotated out of the plane, stresses and tangents change. Further:
  - The tolerances are relative to the largest stress component of the call (1e-12, 1e-9 after 8 iterations) instead
    of absolute (1e-10 / 1e-8 for plane stress, 1e-13 / 1e-10 for uniaxial stress, in the units of the material) or,
    in the gradient-enhanced iteration, relative to max(1, |σ|) with 1e-10 / 1e-5 and a strain criterion; at most 15
    evaluations instead of 13.
  - A singular tangent of the stress-free components `Cff` throws `StressUpdateFailed` (cutback) only if the reduced
    problem is inconsistent (the stress at `p` depends on a strain along the null space of `Cff`, or `σ_f` or `Cfp`
    has a part outside its range; thresholds relative, `ReducedStress::consistencyTolerance`), instead of using a
    compliance clamped to 1e10; the state variables are then reset to the start of the increment. A consistent
    singular `Cff`, e.g. a fully softened Rankine crack without out-of-plane shear stiffness, is condensed with its
    pseudo-inverse (minimum norm corrections: the undetermined strains keep their initial guess).
  - The uniaxial tangent is the condensation of the converged iteration; it condensed all five other components while
    the iteration freed only the lateral normal strains.
  - No `Journal` warning before the throw.
  - `GeneralGradientEnhancedDisplacementFiniteElement` condenses the coupling tangents and the energy gradients in plane
    stress with all three free strains (`ReducedStress::condense()`), not only with ε33.
- `VonMises` and `ADVonMises`: the return mapping converges on a residual relative to the norm of the deviatoric
  trial stress, `|g| <= 1e-12 * ||s_trial||`, instead of an absolute one (`VonMises` 1e-12 MPa, `ADVonMises`
  1e-10 MPa). `g` is a difference of numbers of the size of `||s_trial||`, so its round-off is about
  `eps * ||s_trial||`: the absolute tolerance could not be met for large trial stresses, and large strain increments
  failed to converge depending on the platform's floating point evaluation order (e.g., AVX2 builds, macOS). Both
  models now use the same tolerance; for `ADVonMises` the effective criterion is 1e-12 relative instead of 1e-10
  absolute, i.e. looser for trial stresses above 100 MPa (Newton converges quadratically, so the result is at
  round-off either way). Test `TestVonMisesLargeIncrement`: one strain increment of order one (trial stress about
  900 f_y) ends on the yield surface for perfect, linear and nonlinear hardening (both models, implicit and
  explicit).
- **Results:** coupling tangents of `GeneralGradientEnhancedDisplacementFiniteElement` in plane stress. The blocks
  `kUK`, `kKU` and `kKK` used the derivatives of the material's 3D update without the out-of-plane strain, which the
  plane stress iteration adjusts with the strain and the non-local field. They now include it by static condensation
  of dσ₃₃ = 0, which changes the Newton iterations (not the converged results) of gradient-enhanced plane stress
  analyses whose out-of-plane stress depends on the non-local field or whose local driving variable depends on the
  out-of-plane strain, e.g., `RankineGradientSofteningPlasticity`. A finite difference test of the stiffness in plane
  stress with that material covers it.
- **Results:** energies of `VonMises`. The elastic energy density is ½ σ : C⁻¹ : σ; it was accumulated as
  ½ Δε : C : Δε − 2G Δλ n : Δε per increment, without the work of the stress at the start of the increment. The
  dissipation increment is the plastic work f_y(κₙ₊₁) Δκ of the return mapping instead of 2G Δλ n : Δε. The yield
  stress is floored at zero, so that a softening law flows without stress beyond zero strength instead of reversing
  the stress, and a trial stress without deviatoric part is elastic. This applies to `computeStressExplicit()` too.
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
