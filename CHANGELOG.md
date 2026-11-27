# Changelog

All notable changes to Marmot are documented in this file. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Changes that break existing input or code are marked
**Breaking**; changes that alter numerical results are marked **Results**.

## [Unreleased]

### Added
- Marmot's own automatic differentiation types, without dependencies (`Math::Differentiation::Automatic`):
  - `ScalarMode::Dual< T >` / `ScalarMode::HyperDual< order >` (`MarmotAutomaticDifferentiationScalarMode.h`):
    nested hyper-dual numbers of any order, and `ScalarMode::Taylor< order >`, univariate truncated Taylor numbers.
    They back `dual`, `dual2nd`, `dual3rd`, `HyperDual` and `Taylor`. As fast as or faster than autodiff (GCC 14,
    relative to autodiff: energy Hessian 1.0, third derivative 1.07, phase-field degradation 0.58, B4 and power-law
    Kelvin moduli (Taylor order 5 and 7) 0.92 and 0.53).
  - `VectorMode::Dual< nDirections >` (`MarmotAutomaticDifferentiationVectorMode.h`): vector-mode numbers with
    expression templates that compute values and local partial derivatives once, when the expression is formed, and
    all partial derivatives in one vectorizable loop on assignment. It is the `JacobianDual` of the first-derivative
    drivers `df_dT`, `dF_dT` and `dF_dX< n >`, which therefore compute all directions in one evaluation: tangents of
    `ADVonMises` 3.1x, `FiniteStrainOrthotropicBiotViscoelasticity` 2.6x and `ADCompressibleNeoHooke` 1.9x faster
    than before (one direction per evaluation with autodiff).
- `VectorMode::Dual2nd< nDirections >` (`MarmotAutomaticDifferentiationVectorModeSecondOrder.h`): second-order
  vector-mode numbers carrying value, gradient and Hessian, tracking the directions each number depends on so that
  operations on the seeded inputs only update their Hessian rows. It is the new `HessianDual` of
  `SecondOrder::d2f_dT2`, which therefore evaluates the function once instead of once per pair of directions: the
  Pence-Gou energy Hessian takes 0.65 us instead of 2.4 us (symmetric) or 3.6 us (full), 4x the analytical one.
- `ThirdOrder::d3f_dT3` evaluates the function once per unordered triple of directions (165 instead of 729 for a
  3x3 tensor; with `isSymmetric`, transposed triples are skipped as well), since mixed partial derivatives do not
  depend on the order of differentiation: the Yeoh energy's third derivative, which
  `CompressibleFiniteStrainLinearViscoelasticity` computes in every stress update, takes 11.7 us instead of 55 us.
  Results are unchanged.
- Documentation: the page "Automatic Differentiation" (Features) on the number types, the functions for derivatives
  with respect to scalars, Eigen vectors and Fastor tensors, the AD material base classes and the pitfalls.
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
- Shared hyperelastic strain energy density functions in `MarmotEnergyDensityFunctions.h`, split into
  `EnergyDensityFunctions::Isochoric` (`NeoHookePotential`, `YeohPotential`, `MooneyRivlinPotential`,
  `ArrudaBoyce8ChainPotential`) and `EnergyDensityFunctions::Volumetric` (`PenaltyPotential` and the functions of
  Doll & Schweizerhof (2000), https://doi.org/10.1115/1.321146: `SussmanBathePotential`, `SimoTaylorPotential`,
  `OgdenPotential`, `LiuHofstetterMangPotential`, `DollSchweizerhofPotential`,
  `GeneralizedDollSchweizerhofPotential`, `DollSchweizerhofExponentialPotential`,
  `DollSchweizerhofLogarithmicPotential`), each with a `FirstOrder` overload returning the energy and ∂Ψ/∂C.
  `DecoupledPotential` and `makeDecoupledPotential` combine an isochoric and a volumetric part. The compressible
  `ArrudaBoyce8ChainPotential(C, mu, lambdaL, D1)` is new, as are the `FirstOrder` overloads of the compressible
  `MooneyRivlinPotential` and `YeohPotential`.
- `Math::Differentiation::Automatic::dF_dSymmetricT`: gradient of a tensor-valued function with respect to a
  symmetric square tensor by first-order dual numbers along the dim (dim + 1) / 2 symmetric directions. The output
  may have any shape and the function does not have to be transpose-equivariant; e.g. the closed-form second
  derivative of a potential w.r.t. C yields its third derivative with 6 evaluations.
- Single-pass `DecoupledPotential`: every isochoric and volumetric potential also exists as a scalar function type
  (`Isochoric::NeoHooke`, `Yeoh`, `MooneyRivlin`, `ArrudaBoyce8Chain`; `Volumetric::Penalty`, `SussmanBathe`,
  `SimoTaylor`, `Ogden`, `LiuHofstetterMang`, `DollSchweizerhof`, `GeneralizedDollSchweizerhof`,
  `DollSchweizerhofExponential`, `DollSchweizerhofLogarithmic`), e.g.
  `SecondOrder::DecoupledPotential( C, Isochoric::Yeoh{ C10, C20, C30 }, Volumetric::SimoTaylor{ kappa } )`. The
  plain, `FirstOrder::` and `SecondOrder::` overloads add the scalar derivatives of both parts before the chain rule,
  so the invariants of C are computed once and the tensors assembled once, for any pair. `makeDecoupledPotential`
  accepts these types too. User-defined parts qualify through the concepts `IsochoricScalarFunction` and
  `VolumetricScalarFunction`. The compressible `MooneyRivlinPotential`, `YeohPotential` and
  `ArrudaBoyce8ChainPotential` are such pairs with `Volumetric::Ogden{ 2 / D1, -2 }`.
- `CompressibleFiniteStrainLinearViscoelasticity`: hyperelastic base `ArrudaBoyce` (`baseModel = 4`, elastic
  properties mu, lambdaL, D1); a locking stretch lambdaL <= 1 is rejected.
- `MarmotRightCauchyGreenInvariants.h` (`ContinuumMechanics::Invariants::RightCauchyGreen`): the principal invariants
  (I1, I2, J) and the isochoric invariants (Ī1, Ī2, J) of the right Cauchy-Green tensor with their first and second
  derivatives w.r.t. C (`computeInvariants`, `computeIsochoricInvariants`), and `composeWithInvariants`, the chain
  rule from the scalar derivatives of a function of these invariants to ∂Ψ/∂C and ∂²Ψ/∂C∂C. A hyperelastic
  potential then only states its scalar function and scalar derivatives. The derivatives hold for a general,
  non-symmetric C, as automatic and complex-step differentiation produce them.
- `SecondOrder` overloads (energy, ∂Ψ/∂C and ∂²Ψ/∂C∂C in closed form) of every potential in
  `EnergyDensityFunctions::Isochoric` and `EnergyDensityFunctions::Volumetric`, and of the compressible
  `MooneyRivlinPotential`, `YeohPotential` and `ArrudaBoyce8ChainPotential`. The plain, `FirstOrder` and
  `SecondOrder` overloads of a potential evaluate the same scalar function through the invariant helper.
  The second derivatives of the invariants are kept as a few scalar coefficients (`SecondDerivativeCoefficients`),
  from which the fourth-order tensor is assembled in one pass, and the compressible potentials process C once
  for both of their parts.
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
- `materials/MicroplaneM7`, `materials/GM7`: optional property `stressUnit`, the value of 1 MPa in the user's stress
  unit (default 1.0, e.g. 1e6 for Pa). It converts the dimensional fixed parameters of Caner and Bažant (2013)
  (`c18`, `c21`, the reference strength and modulus), the default values of `fc` and `transitionWidth` and the
  stiffness tolerance of the unloading transition, so that the models can be used in any consistent unit system.
  Results with the default are unchanged.
- `materials/GM7`: optional property 16, `nonlocalMicroInertia` (after the nonlocal viscosity): makes the six nonlocal
  fields second order in time, admissible up to eta^2/4; absent, they stay first order.
- `GM7/COROTATIONAL` (`Materials::GradientEnhancedCorotationalWrapper`): GM7 at finite strains, for the particles
  with six nonlocal variables. The wrapped material is driven by the increment of the logarithmic strain in the frame
  of the polar rotation and returns the rotated Kirchhoff stress, so a rigid rotation rotates the stress exactly also
  for a material that recomputes its stress from total strains (which the Hughes-Winget wrapper cannot); consistent
  analytic tangent through the spectral derivatives of ln C and C^(1/2).
- `Registry::ParticleFactory::isRegistered`: whether a particle name is registered, i.e. whether its module is built
  (replaces the declared but undefined private `checkIfParticleIsRegistered`).
- Explicit dynamics of the gradient-enhanced finite-strain particles (point, SQCNI, SQCNIxNSNI; any number of
  nonlocal variables): `updatePhysicsExplicit`, `computePhysicsKernelsExplicit` (nonlocal balance in total form),
  lumped inertia (density; micro-inertia on the nonlocal rows, zero for a first-order field), lumped damping (nonlocal
  viscosity) and lumped momentum (including the nonlocal rate). `Particle::computeLumpedDamping` (default: no
  damping).
- `Particle::describeStiffness`: the particle describes its stiffness (the tangent of a zero increment of the last
  accepted configuration and the working state, without the inertia) node by node to a `StiffnessVisitor`, in the
  factored form of the chain rule: the shape functions at its evaluation point, then per node the derivatives of its
  force with respect to the interpolated increment, its gradient and its second gradient (`StiffnessOfNode`). It
  leaves the particle as it was found (the state variables the material point was assigned,
  `MaterialPoint::getAssignedStateVars`, and the response and tangents of the material point, `SavedMaterialPoint`),
  also if the visitor throws; the default throws. Implemented in `DisplacementParticle` and
  `GradientEnhancedFiniteStrainParticle` for all derived particles (`addStiffnessOfNode`). The host's visitors, in
  `MarmotParticleStiffnessVisitors.h` and compiled in Marmot: the bound of |K| W, the absolute value of the stiffness
  times nonnegative weights, for the stable time increment of an explicit host (`addBoundOfAbsoluteStiffnessTimes`
  with `AbsoluteStiffnessBound`, in O(nodes) from weighted moments of the shape functions), and the stiffness itself
  (`assembleStiffness` with `StiffnessAssembler`). The particles' tangent is assembled from the same factors
  (`addRowsOfNode`): equal to the former one to round-off, residuals bitwise. The NSNI stabilization rows are shared by
  both particle families (`Meshfree::NSNI`). The material points take their density once, at `initializeYourself`.
- `Particle::getInterpolationVectorLinearizedAboutCenter`: the shape functions linearized about the center of the
  particle, `G_B = N_B(Y_c) + dN_B/dY . (Y - Y_c)`, the particle's own interpolation at its own points (its center, a
  vertex, a face center) for every field: for the displacement the motion of the point in the next increment, for a
  nonlocal field the value the particle carries there. `getInterpolationVector` stays the exact shape functions at an
  arbitrary location. A constraint that measures a field at a point of the particle and distributes its force there
  with `G_B` is symmetric and conservative.
- `Particle::updateShapeFunctions`: evaluate the shape functions anew where the particle is now, with the kernel
  functions it holds, for a host whose search found the same kernel functions as before (the particle has moved, so
  its shape functions have to be evaluated anew in every increment). `assignMeshfreeKernelFunctions` stores the kernel
  functions and calls it, so the two give the same shape functions, bitwise; it is `final` in `GenericParticle`,
  `GenericSDIParticle` and `GradientEnhancedFiniteStrainParticle`, and the SQCNI and NSNI particles override
  `updateShapeFunctions` instead. **Breaking** for implementers of `Particle` outside Marmot: it is pure virtual, and a
  particle derived from one of the three bases that overrode `assignMeshfreeKernelFunctions` moves its evaluation into
  `updateShapeFunctions`.
- `ReproducingKernelApproximation` and `ReproducingKernelApproximationImplicit`: an optional limit on the condition of
  the moment matrix (constructor argument `momentMatrixConditionLimit`, default none). The completeness order of a
  point whose moment matrix, scaled to a unit diagonal, has a larger condition number (1-norm estimate,
  `conditionOfMomentMatrix()`) is lowered one order at a time, down to the order 0, the Shepard functions; so is that
  of a singular moment matrix of order ≥ 1. Such a point lies on one side of nearly coplanar covering kernels, e.g. a
  face of a torn particle reaching into a void, where the linear shape functions grow without bound. Without a limit
  nothing changes; with one, every point below it is evaluated bit for bit as without one. Both approximations choose
  the order of a point in one place (`chooseCompletenessOrderAndAssembleMomentMatrix()`).
- `ReproducingKernelApproximationImplicit` returns the Shepard gradients (the derivatives of its values) at a point of
  order 0, as the explicit approximation does, instead of implicit gradients of order 0, which vanish identically. This
  also applies, without a limit, to a point covered by fewer kernels than the linear basis has entries.
- Axisymmetric RKPM particles, `DisplacementSQCNI/Axisymmetric/Quad` and
  `GradientEnhancedFiniteStrainSQCNIxNSNI/Axisymmetric/Quad`: the twins of the plane-strain SQCNI particles for a body
  of revolution under a load of revolution, in the meridian half plane (r, z) (feature page "Axisymmetric Particles").
  Every integral is taken over the ring the particle sweeps (weight 2πr, ring volume V_0), the hoop stretch
  increment 1 + Δu_r / r of the center enters the material as ΔF_33, and the smoothed gradient gains the area term of
  the ring (two Gauss points per face), so that the linear patch tests pass exactly; a `PRESSURE` is integrated with
  the same two points per face, so that the patch test with a traction on the boundary passes exactly as well; NSNI
  uses the ring's second moments. A particle may touch the axis with a vertex or a face (the ring weight vanishes
  there); one with a vertex at r < 0 is refused. Only the update `DeformationGradient` is registered; VCI and
  `CWFCORRECTION` are refused. The section is a template parameter of the particle classes, `Meshfree::Section`
  (`Cartesian`, the default, or `Axisymmetric`), and of `ParticleSectionIntegrals< section, nDim, nVertices >`, which
  also integrates the face loads (`faceQuadrature`, `addPressure`); the hoop terms of the weak form are
  `Meshfree::Hoop`. `Particle::getSection()` (virtual, default `Section::Cartesian`) tells a host what the integrals of
  a particle measure, e.g. to integrate contact points on the faces of the particles over the same rings. The 2D
  material points (`DisplacementMaterialPoint2D`, `GradientEnhancedFiniteStrainMaterialPoint< 2, nNonlocal >`) take an
  out-of-plane stretch increment ΔF_33 (`addToOutOfPlaneStretchIncrement`, `outOfPlaneStretchIncrement`) and return
  the out-of-plane stress τ_33 and tangents (`OutOfPlaneTangents`, `MarmotOutOfPlaneTangents.h`); plane strain is
  unchanged.

### Changed
- `StressMeasures::FirstOrder::KirchhoffStressAndTangentFromPK2( PK2, dPK2_dC, F )`: the Kirchhoff stress and its
  total derivative dτ/dF for a 2nd Piola-Kirchhoff stress that depends on F through C, evaluated without the
  fourth-order tensors dτ/dS and dC/dF and their two double contractions. `CompressibleHyperelasticity`,
  `CompressibleFiniteStrainLinearViscoelasticity` and `GradientEnhancedCompressibleNeoHookeDamage` use it for their
  tangent: a `computeStress` call of `CompressibleHyperelasticity` is about 1.7 times as fast (about 630 ns to
  375 ns for NeoHooke + DollSchweizerhof) and one of `GradientEnhancedCompressibleNeoHookeDamage` about 1.9 times
  (about 740 ns to 390 ns); for `CompressibleFiniteStrainLinearViscoelasticity` the third derivative of the
  potential dominates, so the gain is below 10 %. Results agree to round-off.
- **Breaking**: the functions `standardNeoHooke`, `SecondOrder::standardNeoHooke` and
  `ThirdOrder::standardNeoHooke` (and with it the namespace `EnergyDensityFunctions::ThirdOrder`) are removed. The
  potential is the scalar function type `Coupled::StandardNeoHooke{ K, G }` of the principal invariants (I1, J),
  evaluated with the new `CoupledPotential( C, potential )`, `FirstOrder::CoupledPotential` and
  `SecondOrder::CoupledPotential` (concept `PrincipalInvariantScalarFunction`), the counterpart of
  `DecoupledPotential` for potentials that are not the sum of an isochoric and a volumetric part. The second
  derivative is the one w.r.t. the nine independent components of C and no longer its minor-symmetric part
  (identical when contracted with symmetric increments).
- **Results**: the `NeoHooke` base of `CompressibleFiniteStrainLinearViscoelasticity` takes its third derivative
  from first-order dual numbers on the closed-form second derivative, like the other bases; results agree to
  round-off and a `computeStress` call takes about 4.4 µs instead of about 70 µs.
- **Breaking:** the material `CompressibleNeoHooke` (`COMPRESSIBLENEOHOOKE`) is replaced by
  `CompressibleHyperelasticity` (`COMPRESSIBLEHYPERELASTICITY`), a hyperelastic material with a decoupled potential
  whose first two material properties select the isochoric part (0 `NeoHooke`, 1 `Yeoh`, 2 `MooneyRivlin`,
  3 `ArrudaBoyce8Chain`) and the volumetric part (1 `SussmanBathe`, 2 `Penalty`, 3 `SimoTaylor`, 4 `Ogden`,
  5 `LiuHofstetterMang`, 6 `DollSchweizerhof`, 7 `DollSchweizerhofExponential`, 8 `DollSchweizerhofLogarithmic`,
  numbered as U1 to U8 in Doll & Schweizerhof (2000)), followed by the parameters of the isochoric part, those of the
  volumetric part and, optionally, the density. Inputs with `COMPRESSIBLENEOHOOKE` and the properties `K, G[, rho]`
  keep their results (to round-off) as `COMPRESSIBLEHYPERELASTICITY` with `0, 6, G, K, 1, 3[, rho]`. The header is
  `Marmot/CompressibleHyperelasticity.h`, the module `CompressibleHyperelasticity`; the examples moved to
  `examples/*/CompressibleHyperelasticity`.
- Faster reproducing kernel shape functions, with bitwise identical results: the kernels give value and gradient from
  one evaluation, the basis of every covering kernel is evaluated once per point (it was evaluated for the moment
  matrix and again for the shape functions), and the moment matrix is summed in local matrices.
- **Breaking:** Automatic differentiation no longer uses autodiff. Marmot code uses the vocabulary of
  `Marmot/MarmotAutomaticDifferentiationTypes.h` in `Math::Differentiation::Automatic`:
  `dual`, `dual2nd`, `dual3rd`, `HyperDual< order >` (was `autodiff::HigherOrderDual< order, double >`),
  `Taylor< order >` (was `autodiff::Real< order, double >`), `VectorXdual`, `VectorXdual2nd`, `MatrixXdual`, and
  `seed< level >`, `derivative< level >`, `value`, `taylorDerivatives< order >` (was `autodiff::derivatives( f,
  along( 1. ), at( x ) )`). `Math::Differentiation::Automatic` no longer re-exports the `autodiff` and `Eigen`
  namespaces, `valnode` is removed (`coefficients< order >` gives the coefficient storage), and
  `Math::Tensors::makeReal` accepts hyper-dual numbers of any order through `isHyperDual`. Results are unchanged.
- **Breaking:** `SecondOrder::d2f_dT2` takes functions of `HessianDual< dim * dim >` numbers instead of `dual2nd`.
  `isSymmetric` of `d2f_dT2`, `df_dT` and `dF_dT` has no effect any more, since all directions are computed in one
  evaluation.
- **Breaking:** autodiff is no longer a dependency: all automatic differentiation uses Marmot's own types (see
  Added), and the conda package no longer depends on autodiff. Results agree with autodiff to about 1e-15
  (relative).
- **Breaking:** The first-derivative drivers `df_dT` and `dF_dT` take functions of `JacobianDual< nDirections >`
  numbers, and `MaterialFiniteStrainAD` and `MaterialHypoElasticAD` pass `ADScalar` numbers
  (`JacobianDual< nDim * nDim >`, `JacobianDual< 6 >`; `Vector6AD`) to `computeStressAD`. These are vector-mode
  numbers, so code differentiated by these drivers must use the new names instead of `dual`.
  `ADLinearElastic::C` is a `Matrix6d`. New: `dF_dX< n >` for fixed-size Eigen vectors (used by
  `MaterialHypoElasticAD`). `df_dT` and `SecondOrder::d2f_dT2` deduce the dimensions from the tensor argument only, so
  they also accept lambdas; `dF_dT` takes the number type as its first template parameter. (Deducing the dimensions
  through the number type of the function crashes MSVC 19.44.)
- All of Marmot's automatic differentiation numbers convert explicitly to `double` (their value), so Fastor's
  `Tensor::cast` converts tensors of them in both directions; `makeReal`, `makeDual`, `makeHigherOrderDual`,
  `makeOtherScalarType` and `fastorTensorFromDoubleTensor` use it.
- **Breaking:** `Element::getStateView` is no longer virtual: it checks the quadrature point index and calls the
  protected `Element::getStateViewAtQuadraturePoint`, which an element implements instead (rename the override).
  Callers are unchanged.
- **Breaking:** Marmot requires Eigen 5; Eigen 3.4 is no longer supported (CMake rejects Eigen < 5.0, the conda package
  requires `eigen >=5`). CI builds and tests with Eigen 5.0.1.
- **Breaking:** `T2D2` is now the 2D 2-node element of the `TrussElement` module instead of the spatial wrapper of
  the 1D `DisplacementFiniteElement` (uniaxial stress). Same property (the cross section area), same small strain
  response for hypoelastic materials; but 1 quadrature point instead of 2 (a deck selecting `quadraturePoint=1`
  must use 0), and the quadrature point states are those of the truss element (`stress`, `strain`,
  `normal force`, ...). `ElementSpatialWrapper` and the uniaxial-stress section of `DisplacementFiniteElement` stay as
  generic facilities (tested), now without a registered element.
- `CompressibleFiniteStrainLinearViscoelasticity` obtains the third derivative of the energy density for the bases
  `Yeoh`, `MooneyRivlin`, `PenceGouNeoHooke` and `ArrudaBoyce` by differentiating the closed-form `SecondOrder`
  potentials with first-order dual numbers (6 evaluations) instead of applying third-order dual numbers to the plain
  potentials (6 x 81 evaluations), using `dF_dSymmetricT`. Results agree to round-off; this evaluation is about 25
  to 35 times faster (about 2.5 µs instead of 60 to 95 µs per material point).
- `CompressibleFiniteStrainLinearViscoelasticity` takes the initial shear modulus for `onlyShearCreep = 1` from one
  function, `initialShearModulus()`, with the analytical expression of each hyperelastic base; the Arruda-Boyce one
  is written as `mu (3 lambdaL^2 - 1) / (3 (lambdaL^2 - 1))`. For full creep the initial compliance is the closed
  form `1/(2 G0) Dev + 1/(9 K0) I(x)I` with the analytical `initialBulkModulus()`, instead of the numerically
  inverted initial stiffness from automatic differentiation. Results agree to round-off.
- `FirstOrder::PenceGouPotentialB` and `SecondOrder::PenceGouPotentialB` evaluate their derivatives through the
  invariant helper. Results agree with the previous closed forms to round-off; the `FirstOrder` overload takes about
  half the time, the `SecondOrder` overload the same. The plain `PenceGouPotentialB` is unchanged.
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
- **Breaking:** the gradient-enhanced finite-strain stack carries `nNonlocal` nonlocal variables:
  `Materials::MaterialGradientEnhancedFiniteStrain<nNonlocal>` and `Registry::MaterialGradientEnhancedFiniteStrainFactory
  <nNonlocal>` (instantiated for 1 and 6), the material point `GradientEnhancedFiniteStrainMaterialPoint<nDim,
  nNonlocal>` (one class; `GradientEnhancedFiniteStrainMaterialPoint2D`/`3D` are removed) and the particles
  `GradientEnhancedFiniteStrainParticle<nDim, nNonlocal>`, `...SQCNI<nDim, nVertices, nNonlocal>`,
  `...SQCNIxNSNI<nDim, nVertices, nNonlocal>`. Nonlocal quantities are Eigen vectors (`L`, `nonLocalRadius`, `N`,
  `getNonlocalViscosity`, `getNonlocalMicroInertia`), the mixed tangents one tensor per nonlocal variable
  (`dTau_dN[K]`, `dL_dF[K]`). The material point keeps the total
  gradient and the rate of the nonlocal fields (new states `nonlocal damage gradient`, `nonlocal damage rate`).
  Registered with six nonlocal variables (`strain symmetric`, 3D): `GradientEnhanced6FiniteStrain/3D`,
  `GradientEnhanced6FiniteStrain/3D/Point`, `GradientEnhanced6FiniteStrainSQCNI/3D/Hexa`,
  `GradientEnhanced6FiniteStrainSQCNIxNSNI/3D/Hexa`. The cell and the finite element carry one nonlocal variable. The
  Hughes-Winget wrapper shares its base `GradientEnhancedHypoElasticWrapper` with the corotational one. Results with
  one nonlocal variable are unchanged (bitwise).
- **Breaking:** `Particle::updatePhysicsExplicit`, `Particle::computePhysicsKernelsExplicit`,
  `Particle::computeBodyLoadExplicit` and `Particle::computeDistributedLoadExplicit` throw by default instead of doing
  nothing, so that a particle without the explicit path is refused by an explicit host instead of being integrated with
  its mass and without its force, and never leaves a load out silently. The explicit path is in the NSNI particles of
  the displacement family and in all gradient-enhanced particles; `DisplacementParticleSQCNI` implements the explicit
  distributed load only. The defaults of the lumped inertia and momentum name the method in their message.
- **Breaking** for implementers of `Particle` outside Marmot: `Particle::getInterpolationVectorLinearizedAboutCenter` is pure
  virtual, so a particle must implement it to compile (there is no safe default: the shape functions at the point are
  not the particle's own interpolation).
- `Math::MonomialBasis` multiplies only the factors of nonzero exponent of an entry, listed at construction, and sets
  the derivative with respect to a coordinate that is not a factor to zero: the values and nonzero derivatives are
  unchanged (bitwise; a zero derivative is +0).
- **Results:** `ReproducingKernelApproximation` and `ReproducingKernelApproximationImplicit` evaluate each kernel once
  at a point (its gradient only if needed), with one assembly of the moment matrix (and its gradients) for both
  classes, and types of a fixed largest size (up to the cubic basis in 3D): no allocation per evaluation. A larger
  basis, a dimension beyond 3 or a negative order is refused at construction (`std::invalid_argument`). The monomial
  basis is `Math::MonomialBasis`, a table of the exponents of its entries made once and evaluated by products
  (`Math::computeMonomialBasis` and `computeMonomialBasisGradient` use it, take any contiguous Eigen vector or matrix,
  `Eigen::Ref`, and give the same values); the moment matrix is assembled entry by entry. Shape functions change by
  round-off.
- `GM7` overrides `computeStressExplicit`: the stress update without homogenizing the tangents of the microplanes, one
  update with the tangents an optional output (`stress_update`); stress and state variables are those of
  `computeStress`, bitwise. The explicit path of `GeneralGradientEnhancedDisplacementFiniteElement` takes it.
- One name for one quantity: the limit 2 / omega_max of the central difference method is the critical time step
  everywhere (documentation, messages, and `CriticalTimeStep::IsotropicElasticStiffnessBound::criticalTimeStep`, was
  `stableTimeIncrement`; the test helper `Testing::exactCriticalTimeStepOfTheElement`, was
  `stableTimeIncrementOfTheElement`). **Breaking** for code that calls the two functions. The momentum buffer of
  `Particle::computeLumpedMomentum` is named `pLumped` (was `mLumped`, the name of the mass buffer).

### Deprecated
- `BoundaryElement::computeSurfaceNormalVectorialLoadVectorForAxisymmetricElements()`: construct the boundary element
  with `axisymmetric = true` and call `computeSurfaceNormalVectorialLoadVector()`.

### Removed
- **Breaking:** `Approximation::computeShapeFunctionGradients()`, which both reproducing kernel approximations only
  threw from; and their protected `computeHVector`, `computeHVectorGradient`, `computeMMatrix`,
  `computeMMatrixAndGradient`, `H0Vector` and `findCoveringKernelFunctionIndices`, replaced by
  `evaluateKernelFunctions`, `evaluateBasisOfKernels` and `assembleMomentMatrix`.
- **Breaking:** `Meshfree::KernelFunction::computeKernelFunctionGradient()`, replaced by
  `computeKernelFunctionAndGradient()`, which returns the value with the gradient (the B-spline kernels evaluate their
  one-dimensional factors once for both; `computeBSpline2ndOrder/3rdOrder` give the derivative on request, their
  `...Gradient` functions are gone).
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
- **Results** The axisymmetric finite strain elements `CX8UL` and `CX8RUL` integrate their lumped and consistent mass
  and their body force over the ring, with the weight 2 pi r of their internal force, through one section factor
  (`DisplacementFiniteStrainULElement::sectionFactorAt`, the ring instead of the thickness). The mass and the body
  force were the plane ones, so force and mass referred to different measures and every dynamic analysis with these
  elements had the wrong inertia, by a factor 2 pi r that varies over the mesh; the critical time step, built from
  the mass the element assembles, changes with it. The stored elastic energy and dissipation, read back divided by
  the plane weight, were off by the same factor. The thickness property is no longer used by them.
- **Results** The distributed loads of `CX8UL` and `CX8RUL` (`DisplacementFiniteStrainULElement::computeDistributedLoad`)
  are integrated over the ring 2 pi r of the face, as those of `CAX*` and `GCAX*`: a pressure on the current face,
  with the derivative of 2 pi r in its stiffness, a surface traction on the reference face. They were integrated over
  the face times the thickness property, so a surface load on these elements was off by the factor 2 pi r / thickness
  against their internal force. `BoundaryElement::computeDSurfaceNormalVectorialLoadVector_dCoordinates` supports
  axisymmetric boundary elements.
- The critical time step of the continuum elements (`DisplacementFiniteElement`, `DisplacementFiniteStrainULElement`,
  `GeneralGradientEnhancedDisplacementFiniteElement`) refuses a wave speed or a density that is not a positive number
  at any quadrature point (`CriticalTimeStep::pWaveModulusOfQuadraturePoint`), with the element and the point in the
  message. A NaN passed the guard `c <= 0` and dropped out of the largest P-wave modulus, which made the time
  increment too large. `criticalTimeStep` checks every entry of the lumped mass (Eigen's `minCoeff` may skip a
  NaN) and the convergence of the eigenvalue solver; the default `MaterialFiniteStrain::getMaximumWaveSpeed` returns
  NaN for a NaN stiffness of one direction instead of the wave speed of the others.
- `CX8UL` and `CX8RUL` evaluate the wave speed of the critical time step at their hoop stretch
  F_tt = 1 + u_r / r, the deformation their kernels integrate, instead of at F_33 = 1; they refuse a node at r < 0 or
  a quadrature point at r <= 0 when their coordinates are assigned, and a hoop stretch that is not positive in an
  increment (`std::runtime_error`).
- `Element::getStateView` throws `std::invalid_argument` (naming the state, the index and the number of quadrature
  points of the element) for a quadrature point outside `0` to `getNumberOfQuadraturePoints() - 1`, instead of
  reading past the end of the element's quadrature points; e.g. an output asking for quadrature points 0 to 8 of a
  4-point element segfaulted.
- **Results** `HughesWingetWrapper` hands the stored elastic energy density and dissipation (per unit reference
  volume) back to the wrapped material (per unit current volume) converted, as an energy of the state divided by
  J_n = det(F_n) and as a dissipation accumulated along the path divided by J = det(F_{n+1}), so that
  D_{n+1} = D_n + J * (the wrapped material's increment). Without it, a wrapped material that adds to the value it
  is given had its stored energy and dissipation multiplied by J in every increment, compounding:
  `VONMISES/HUGHES-WINGET` (and its exact tangent variant), `LINEARELASTIC/HUGHES-WINGET`,
  `RANKINEDAMAGEPLASTICITY/HUGHES-WINGET` and `RANKINESOFTENINGPLASTICITY/HUGHES-WINGET`. The conversions
  (`Materials::toDensityPerVolumeAtStartOfIncrement`, `toDensityPerVolumeAtEndOfIncrement`,
  `toDensityPerReferenceVolume`) are shared with `GradientEnhancedHughesWingetWrapper`, which divided both by J_n:
  its dissipation is now also accumulated along the path, `D_{n+1} = D_n + J * increment`, instead of rescaled by
  J/J_n in every increment. Stresses, tangents and state variables are unchanged; only the reported energies
  change, at J != 1. The interface documents both densities as per unit reference volume.
- **Breaking:** the gradient-enhanced SQCNI/SNNI particles (`GradientEnhanced[6]FiniteStrainSQCNI...`,
  `...SNNI...`, `..._R...`, `..._RU...`) report the vertices and face centers of their geometry, mapped with the
  deformation gradient, as the displacement particles do, not those of the smoothing domain (which differ for SNNI,
  `_R` and `_RU`): the vertices then move with `getInterpolationVectorLinearizedAboutCenter`. The geometry is stored in the
  state (one more block of `nDim x nVertices` values); the state views are `vertex displacements` (geometry) and
  `smoothing vertex displacements` (smoothing domain, formerly `vertex displacements`). The physics is unchanged.
- **Results** The critical time step of the continuum elements (`DisplacementFiniteElement`,
  `DisplacementFiniteStrainULElement` incl. its axisymmetric variant, and the displacement field of
  `GeneralGradientEnhancedDisplacementFiniteElement`) is the element's own eigenvalue bound
  (`Discretizations::Elements::CriticalTimeStep`), from the strain operator the element integrates with and the lumped
  mass it assembles. The former characteristic length over wave speed overestimated the stable increment of quadratic
  elements by a factor of 2.3 to 3.4 (explicit runs with them diverged at courant numbers of 0.5 and the default 0.8)
  and of compact linear elements by up to 1.7. Explicit runs take smaller increments; quadratic elements now run at
  the default courant number.
- **Breaking** `GeometryElement::characteristicLengthForCriticalTimeStepAt` is now `smallestElementExtentAt` (the
  artificial bulk viscosity is its only user); `MassLumping::manifoldMassFractions` and
  `MassLumping::timeStepFactorFromMassDistribution` are removed.
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
- **Results:** `CompressibleFiniteStrainLinearViscoelasticity` with `baseModel = MooneyRivlin` and
  `onlyShearCreep = 1` takes the reference shear modulus 2 (C1 + C2) of the Mooney-Rivlin potential; it took 2 C1
  before. Results change for inputs with C2 ≠ 0.
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
