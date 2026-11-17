# Automatic differentiation backends for Marmot

Marmot's automatic differentiation (AD) goes through one header, `Marmot/MarmotAutomaticDifferentiationTypes.h`. Its
backend-neutral vocabulary in `Marmot::Math::Differentiation::Automatic` is:

* types: `dual`, `dual2nd`, `dual3rd`, `HyperDual< order >`, `Taylor< order >`
* functions: `seed< level >`, `derivative< level >`, `value`, `taylorDerivatives< order >`

The header is backed by [autodiff](https://github.com/autodiff/autodiff). autodiff's last commit was in January 2025,
and 1.1.2 needs the unmerged patch autodiff#397 to compile with Eigen 5.

This benchmark compares every AD library that can technically and legally replace that backend, on Marmot's AD code
paths:

* Each backend is wrapped in an adapter with the vocabulary above.
* The autodiff adapter *is* the Marmot types header.
* The kernels and drivers in `ADKernels.h` are written once for all backends. They call the actual Marmot templates:
  the Pence-Gou and Yeoh energies, `rightCauchyGreen`, and `KirchhoffStressFromPK2`.

Three differentiation modes appear in the results:

* **scalar mode**: one direction per pass, nested for higher orders. This is what Marmot does today (`dual`, `dual2nd`,
  `dual3rd`).
* **vector mode** (`-vec`): all 9 (or 6) directions in one pass, nested for the Hessian and the third derivative, e.g.
  `Sacado::Fad::SFad< double, 9 >`. Using it in Marmot would mean vector-mode versions of the drivers `dF_dT`,
  `d2f_dT2`, `d3f_dT3` and `dF_dX`; `ADKernels.h` contains them.
* **reverse mode** (`-rev`): one taped evaluation, then the Jacobian from sweeps over the tape (Adept).

## Candidates

Marmot is LGPL-2.1-or-later and is linked into proprietary hosts (e.g. Abaqus). A backend therefore has to be under a
permissive license, MPL-2.0 or LGPL. It also has to compile with GCC, clang and MSVC, and work with Fastor tensors and
Eigen.

| library | license | maintained | verdict |
|---|---|---|---|
| autodiff 1.1.2 + autodiff#397 | MIT | last commit Jan 2025 | **reference** (current backend) |
| [Boost.Math `fvar`](https://www.boost.org/doc/libs/release/libs/math/doc/html/math_toolkit/autodiff.html) (Boost 1.87) | BSL-1.0 | yes | benchmarked |
| [Sacado](https://github.com/trilinos/Trilinos/tree/master/packages/sacado) (Trilinos 17.3.0) | BSD-3 | yes | benchmarked |
| [Ceres `Jet`](http://ceres-solver.org/automatic_derivatives.html) 2.2.0 | BSD-3 | yes | benchmarked, **needs a patched `jet.h`** |
| [Stan Math](https://github.com/stan-dev/math) v5.4.0 (`fvar`) | BSD-3 | yes | benchmarked |
| [TinyAD](https://github.com/patr-schm/TinyAD) (4b48d1a, no releases) | MIT | slow (last commit May 2025) | benchmarked |
| [Adept-2](https://github.com/rjhogan/Adept-2) v2.1.3 | Apache-2.0 | yes | benchmarked, **wrong Jacobians through Fastor `einsum`** |
| Eigen `AutoDiffScalar` (unsupported module) | MPL-2.0 | with Eigen | excluded: does not compile with Marmot's Fastor code (see below) |
| [Enzyme](https://github.com/EnzymeAD/Enzyme) | Apache-2.0 + LLVM exception | yes | excluded: clang plugin for one exact LLVM version, no GCC/MSVC; differentiates functions, not number types |
| [Clad](https://github.com/vgvassilev/clad) | LGPL-3.0 | yes | excluded: same plugin and architecture limitations as Enzyme |
| [CasADi](https://github.com/casadi/casadi) | LGPL-3.0 | yes | excluded: symbolic graph per evaluation, not a number type for per-integration-point code |
| CppAD | EPL-2.0 or GPL-2.0+ | yes | excluded by license: EPL adds restrictions LGPL-2.1 forbids; GPL would make libMarmot GPL |
| ADOL-C | EPL-1.0 or GPL-2.0+ | yes | excluded by license (as CppAD) |
| CoDiPack | GPL-3.0 | yes | excluded by license |
| audi | GPL-3.0 | yes | excluded by license |
| XAD | AGPL-3.0 | yes | excluded by license |
| cppduals, FADBAD++ | – | no (2019 / 2012) | excluded: unmaintained, so no improvement over autodiff |

## Build and run

```bash
# fetch all backends at the benchmarked versions; prints the CMake arguments
args=$(benchmarks/automatic-differentiation/fetchBackends.sh /tmp/ad-backends $CONDA_PREFIX/include)
cmake -S benchmarks/automatic-differentiation -B build-ad-bench -DCMAKE_BUILD_TYPE=Release \
      -DEigen3_DIR=$CONDA_PREFIX/share/eigen3/cmake $args
cmake --build build-ad-bench
./build-ad-bench/ad-backend-benchmark [batches=15] [min-seconds-per-batch=0.05]   # kernels + correctness
./build-ad-bench/ad-op-benchmark                                                  # per-operation cost
./build-ad-bench/adept-fastor-einsum-repro                                        # Adept + Fastor einsum issue
```

Notes:

* Everything is header-only from source trees. The exception is Adept's small runtime, which is compiled from its
  sources.
* Sacado needs `sacado-config/Sacado_config.h` instead of the configured one.
* CMake applies `ceres-2.2.0-jet-fastor.patch` to a copy of Ceres' headers.
* Boost.Math comes from Stan's bundled Boost 1.87, so the executable contains one Boost version only.

## Results

Setup:

* AMD Threadripper 7970X, one pinned core, `-O3 -march=native`.
* Eigen 5.0.1, Fastor 0.6.4.post6 (conda), Marmot `upstream/master` plus the backend-neutral types.
* Raw output: `results/`.

Each value is the time of one kernel call relative to autodiff scalar mode, which is the current Marmot
implementation; **< 1 is faster**. The autodiff column gives absolute times. *n/a* means the library cannot express
the kernel: no Taylor mode, no nesting, or no vector mode.

All backends agree with autodiff to a relative deviation ≤ 7.4e-14 (≤ 2e-15 except for Taylor order 5). The one
exception is Adept's finite-strain tangent, which is wrong and not run (see below). Run-to-run noise is about ±10 %.

Kernels and the Marmot code paths they mirror:

| kernel | Marmot code path |
|---|---|
| `tangent` | dτ/dF of `ADCompressibleNeoHooke` (`MaterialFiniteStrainAD`, `dF_dT`, 9 directions) |
| `hypoJ2` | tangent of an ADVonMises-like return mapping (`MaterialHypoElasticAD`, `dF_dX`, Eigen, 6 directions) |
| `hessian` | d²ψ/dC² of the Pence-Gou energy (`SecondOrder::d2f_dT2`) |
| `third` | d³ψ/dC³ of the Yeoh energy (`ThirdOrder::d3f_dT3`; `CompressibleFiniteStrainLinearViscoelasticity` calls it every stress update) |
| `phasefield` | generic phase-field degradation and its derivatives (`dual2nd` scalar) |
| `b4` | B4 drying-creep Kelvin moduli, Post-Widder of order 5 (Taylor mode; every B4 stress update) |
| `powerlaw` | power-law Kelvin moduli, Post-Widder of order 7 (Taylor mode) |

### GCC 14.3 (Marmot's main compiler)

| backend | tangent | hypoJ2 | hessian | third | phasefield | b4 | powerlaw |
|---|---:|---:|---:|---:|---:|---:|---:|
| autodiff [ns] | 1401 | 380 | 3537 | 71198 | 16 | 898 | 993 |
| boost | 1.26 | 1.15 | 1.99 | 7.64 | 4.70 | 3.35 | 2.11 |
| sacado | 1.13 | 1.04 | 2.08 | 3.52 | 3.54 | 3.76 | 2.40 |
| **sacado-vec** | **0.62** | **0.44** | **0.75** | **0.74** | n/a | n/a | n/a |
| ceres | 1.15 | 0.96 | 1.93 | 9.02 | 5.66 | n/a | n/a |
| ceres-vec | 3.04 | 1.68 | 2.08 | 3.08 | n/a | n/a | n/a |
| stan | 1.08 | 1.09 | 1.33 | 1.76 | 0.69 | n/a | n/a |
| tinyad | 1.44 | 1.16 | n/a | n/a | n/a | n/a | n/a |
| tinyad-vec | 2.85 | 0.39 | 0.41 | n/a | n/a | n/a | n/a |
| adept-rev | wrong | 5.18 | n/a | n/a | n/a | n/a | n/a |

### clang 22.1

| backend | tangent | hypoJ2 | hessian | third | phasefield | b4 | powerlaw |
|---|---:|---:|---:|---:|---:|---:|---:|
| autodiff [ns] | 1312 | 293 | 3921 | 105837 | 13 | 883 | 929 |
| boost | 1.06 | 1.20 | 1.07 | 3.49 | 4.23 | 1.35 | 1.17 |
| sacado | 1.38 | 1.03 | 2.13 | 1.75 | 10.60 | 3.30 | 2.51 |
| sacado-vec | 1.06 | 0.64 | 0.99 | 0.64 | n/a | n/a | n/a |
| ceres | 1.07 | 1.03 | 1.30 | 1.15 | 1.00 | n/a | n/a |
| ceres-vec | 1.68 | 0.59 | 0.78 | 0.92 | n/a | n/a | n/a |
| stan | 1.23 | 1.06 | 1.35 | 1.05 | 1.38 | n/a | n/a |
| tinyad | 0.88 | 1.06 | n/a | n/a | n/a | n/a | n/a |
| **tinyad-vec** | **0.78** | **0.51** | **0.30** | n/a | n/a | n/a | n/a |
| adept-rev | wrong | 6.71 | n/a | n/a | n/a | n/a | n/a |

`ad-op-benchmark` (in `results/`) shows the cost per elementary operation. All alternatives match autodiff on `*` and
`sqrt`. Their nested (2nd/3rd order) `pow( x, double )` is 2–5× slower, and Marmot's energies use it heavily, e.g.
`pow( J, -2. / 3 )`.

## Findings per library

**Sacado**
* Vector mode is the only option that is faster than autodiff with GCC on every first-, second- and third-order kernel
  (0.44–0.75×). With clang it ranges from 0.64× to 1.06×.
* Scalar mode is slower than autodiff.
* Its Taylor type has a run-time degree and lives on the heap, which makes it 2.4–3.8× slower than autodiff's `Real`.
* Works with Fastor and Eigen as is. The adapter only adds Eigen `NumTraits`.
* Shipped only inside Trilinos. Header-only use needs `packages/sacado/src` (2.1 MB) and a config stub.

**TinyAD**
* Very fast in vector mode on Eigen code (`hypoJ2` 0.39× with GCC) and for full Hessians, which it computes directly
  in one pass (0.41× GCC, 0.30× clang).
* The Fastor-based `tangent` is slow with GCC (2.85×) but fast with clang (0.78×).
* Does not nest: no mixed hyper-duals, no third derivatives, no Taylor mode. It cannot replace `dual2nd`/`dual3rd`
  scalar code.
* Needs a unary `operator+` for Fastor's `inverse()`, which the adapter adds. Maintenance is slow, and there are no
  releases.

**Stan Math**
* The closest scalar-mode drop-in: 1.08–1.76× with GCC, 1.05–1.38× with clang, and the phase-field kernel 0.69× with
  GCC.
* Works with Fastor and Eigen without glue.
* No vector mode and no Taylor mode.
* Brings Stan's Boost and TBB headers along; the repository is 242 MB.

**Ceres `Jet`**
* Upstream `jet.h` does not compile with Fastor (explicit value constructor). Nested Jets also lack `double` arithmetic.
  The benchmark carries a patch plus adapter operators.
* Vector mode is much slower with Eigen 5 than it was with Eigen 3.4 (GCC tangent 3.04×, previously 1.54×).

**Boost.Math**
* Slow nested hyper-duals with GCC (2–8×), and Taylor mode 2–3× slower with GCC.
* No vector mode.

**Adept-2**
* Recording through Fastor's `einsum` gives correct values but **wrong Jacobians**. Even a two-operand
  `einsum< iI, IJ, to_ij >` fails, see `repro/adept-fastor-einsum.cpp`. Marmot's finite-strain code contracts with
  `einsum` everywhere, so this rules Adept out.
* On Eigen code the Jacobian is correct but 5–7× slower than autodiff: taping costs more than forward mode for 6×6
  Jacobians.

**Eigen `AutoDiffScalar`** (excluded)
* It declares `void operator+() const`, so Fastor's `inverse()` (`dst[0] = + a * b - …`) does not compile.
* Nested scalars have no `double` arithmetic.
* Eigen's own `pow( AutoDiffScalar, double )` calls `pow( x.value(), y )` after `using std::pow`, which GCC 14's
  libstdc++ rejects with a hard error for nested types.

## Recommendation

**Sacado in vector mode** is the only candidate that beats the current implementation with GCC on all first-, second-
and third-order kernels. It needs no patches, has a long-term maintainer (Sandia; deal.II uses it for constitutive
models), and is permissively licensed.

With the backend-neutral types in place, the next steps are:

1. Add vector-mode drivers (`dF_dT`, `dF_dX`, `d2f_dT2`, `d3f_dT3`) on `Sacado::Fad::SFad< …, N >` to the types header.
   The material code only sees the scalar type.
2. Keep Taylor mode on autodiff's `Real`, vendored. It is faster than every alternative, and none of the
   vector-mode candidates has a fixed-order Taylor type.
3. Decide how to ship Sacado: vendor `packages/sacado/src` with the config stub, or depend on an installed Trilinos.

Alternatives:

* **TinyAD** is attractive for first- and second-order vector mode, especially full Hessians. It cannot cover the
  nested third-order code, though, and its GCC performance on Fastor code is poor.
* **Stan** is the only near-drop-in replacement for today's scalar-mode code if vector mode is not adopted. It costs
  about 10–75 % speed and brings a large dependency.
