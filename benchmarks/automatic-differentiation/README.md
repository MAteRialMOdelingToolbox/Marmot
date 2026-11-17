# Forward-mode AD backends for Marmot

Compares candidate replacements for [autodiff](https://github.com/autodiff/autodiff) on Marmot's AD code paths.
autodiff's last commit was in January 2025, and its "Update to Eigen 5" issue (#392) has been open since February 2026.
The candidates:

* [Boost.Math `autodiff`](https://www.boost.org/doc/libs/release/libs/math/doc/html/math_toolkit/autodiff.html) (`fvar`), Boost 1.88
* [Sacado](https://github.com/trilinos/Trilinos/tree/master/packages/sacado) (`Fad::SFad`, `Tay::Taylor`), Trilinos 17.3.0, header-only
* [Ceres `Jet`](http://ceres-solver.org/automatic_derivatives.html) (`ceres::Jet<T, N>`), ceres-solver 2.2.0, header-only, **patched** (see below)

Every backend is wrapped in an adapter with one vocabulary.
`ADKernels.h` therefore contains the physics and the drivers only once.
The kernels call the actual Marmot templates: Pence-Gou/Yeoh energies, `rightCauchyGreen` and `KirchhoffStressFromPK2`.

Two differentiation modes are measured:

* **scalar mode**: one direction per pass, nested for higher orders. This is what Marmot does today (`autodiff::dual`, `dual2nd`, `dual3rd`).
* **vector mode** (`-vec`): all directions in one pass.
  Examples: `SFad<double, 9>` for the 9 components of F; `SFad<SFad<double, 9>, 9>` for a full Hessian in one evaluation.
  Sacado and Ceres support this; autodiff and Boost do not.
  Using it means changing Marmot's driver functions (`dF_dT`, `d2f_dT2`, …); `ADKernels.h` shows the vector-mode versions.

## Build and run

```bash
cmake -S benchmarks/automatic-differentiation -B build-ad-bench -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_PREFIX_PATH="<eigen-3.4 prefix>;<conda env>" \
      -DSACADO_SOURCE_DIR=<Trilinos>/packages/sacado \
      -DCERES_SOURCE_DIR=<ceres-solver 2.2.0>
cmake --build build-ad-bench
./build-ad-bench/ad-backend-benchmark [batches=15] [min-seconds-per-batch=0.05]   # kernels + correctness
./build-ad-bench/ad-op-benchmark                                                  # per-operation cost
```

* Sacado: only its source tree is needed. `sacado-config/Sacado_config.h` stands in for the CMake-generated config (no Kokkos, Teuchos or MPI).
  A sparse checkout is enough: `git clone --depth 1 --branch trilinos-release-17-3-0 --filter=blob:none --sparse … && git sparse-checkout set packages/sacado`.
* Ceres: CMake copies `include/ceres` into the build tree and applies `ceres-2.2.0-jet-fastor.patch`.
* Eigen 3.4 is required, because autodiff 1.1.2 does not compile against Eigen 5.

## Results

Setup: AMD Threadripper 7970X, one pinned core, `-O3 -march=native`, Eigen 3.4.0.
Each value is the kernel time relative to autodiff scalar mode, the current Marmot implementation. **< 1 is faster.**
n/a means the library has no Taylor mode, or vector mode does not apply.
All results agree with autodiff to ≤ 7.4e-14 relative (≤ 2e-15 except for Taylor 5).
Run-to-run noise is about ±15 %.

### GCC 14.3 (Marmot's main compiler)

| kernel (Marmot path)                                       | autodiff [ns] | boost | sacado | **sacado-vec** | ceres | ceres-vec |
|------------------------------------------------------------|------:|-----:|------:|------:|------:|------:|
| finite-strain tangent dTau/dF (`MaterialFiniteStrainAD`)   |  1716 | 1.36 | 1.24 | **0.75** | 1.31 | 1.54 |
| hypoelastic J2 tangent, Eigen (`MaterialHypoElasticAD`)   |   480 | 1.15 | 0.94 | **0.44** | 0.95 | 1.27 |
| energy Hessian d²ψ/dC² (`d2f_dT2`)                         |  4101 | 1.97 | 3.06 | **0.82** | 2.45 | 1.22 |
| energy 3rd derivative d³ψ/dC³ (`d3f_dT3`, every stress update) | 84466 | 9.35 | 4.47 | **0.94** | 7.40 | 1.74 |
| phase-field degradation (dual2nd scalar)                   |    25 | 3.97 | 3.59 | n/a | 5.79 | n/a |
| B4 drying Kelvin moduli (Taylor 5, every stress update)    |  1294 | 2.98 | 3.59 | n/a | n/a | n/a |
| power-law Kelvin moduli (Taylor 7)                         |  1114 | 2.81 | 3.04 | n/a | n/a | n/a |

### clang 22.1

| kernel                                   | autodiff [ns] | boost | sacado | sacado-vec | ceres | **ceres-vec** |
|------------------------------------------|------:|-----:|------:|------:|------:|------:|
| finite-strain tangent dTau/dF            |  1738 | 0.97 | 1.54 | 1.15 | 1.13 | **0.48** |
| hypoelastic J2 tangent, Eigen            |   406 | 1.21 | 1.01 | 0.65 | 1.07 | **0.57** |
| energy Hessian d²ψ/dC²                   |  5265 | 0.91 | 2.26 | 0.88 | 1.12 | **0.35** |
| energy 3rd derivative d³ψ/dC³            | 141880 | 3.27 | 1.93 | 0.65 | 1.25 | **0.36** |
| phase-field degradation (dual2nd scalar) |    16 | 1.04 | 10.90 | n/a | 0.96 | n/a |
| B4 drying Kelvin moduli (Taylor 5)       |  1022 | 1.67 | 3.28 | n/a | n/a | n/a |
| power-law Kelvin moduli (Taylor 7)       |  1266 | 1.16 | 2.41 | n/a | n/a | n/a |

### Per operation (GCC, ns per chained op, `ad-op-benchmark`)

| type    | backend  | mul | div | sqrt | exp∘log | pow∘pow |
|---------|----------|----:|----:|-----:|--------:|--------:|
| dual    | autodiff | 3.6 | 3.6 | 10.7 | 21.0 | 41.9 |
| dual    | sacado   | 3.7 | 6.7 | 10.9 | 24.5 | 97.7 |
| dual    | ceres    | 3.7 | 4.3 | 11.4 | 19.6 | 55.7 |
| dual2nd | autodiff | 5.2 | 4.9 | 14.4 | 20.8 | 47.8 |
| dual2nd | sacado   | 4.8 | 9.1 | 16.6 | 36.4 | 117.0 |
| dual2nd | ceres    | 5.2 | 4.7 | 7.6 | 18.6 | 103.7 |
| dual3rd | autodiff | 8.8 | 13.3 | 20.6 | 25.9 | 54.4 |
| dual3rd | sacado   | 7.6 | 16.2 | 25.6 | 86.0 | 366.7 |
| dual3rd | ceres    | 8.7 | 13.0 | 16.3 | 28.4 | 281.5 |
| Taylor5 | autodiff | 7.0 | 10.1 | 22.1 | 26.8 | 70.7 |
| Taylor5 | sacado   | 40.6 | 74.2 | 94.4 | 96.8 | 282.4 |

The Boost rows and the clang table are printed by `ad-op-benchmark`.

## Findings per library

**Sacado**
* Vector mode is the only option that is faster than autodiff with GCC, and it is faster on every kernel it covers:
  0.44–0.94× (0.65–1.15× with clang).
* Scalar mode is slower than autodiff, mainly because `pow(x, double)` is 2–7× slower and Marmot's energies are full of `pow(J, -2./3)`.
* `Tay::Taylor` has a run-time degree and allocates on the heap, so it is 2.4–3.6× slower than autodiff's `Real<N>`.
  It is not a replacement for the Kelvin-chain Taylor mode.
* Works with Fastor and the Marmot templates unchanged.
  Eigen needs `NumTraits` and `ScalarBinaryOpTraits` (in `ADBackendSacado.h`).
* Dependency: Sacado ships only inside Trilinos. Header-only use needs the 2.1 MB `packages/sacado/src` tree and a config stub.
  It could be vendored (BSD-3) or taken from a Trilinos install.

**Ceres `Jet`**
* **Upstream Ceres does not work with Fastor.**
  Jet's value constructor is `explicit`, and Fastor does `T x = 0.` and `x = 0` internally.
  Nested Jets (2nd/3rd order) also lack `double ⊕ Jet<Jet>` arithmetic and `Jet<Jet>(double)` construction.
  The numbers above use a patched `jet.h`: the constructor made implicit, plus one constructor from arithmetic types
  (`ceres-2.2.0-jet-fastor.patch`).
  The mixed `double` operators are added in `ADBackendCeres.h`.
  Using Ceres would therefore mean carrying a fork of `jet.h`.
* Strongly compiler-dependent.
  Vector mode is the fastest option with clang (0.35–0.57×) but slower than autodiff with GCC (1.2–1.7×),
  because Jet's derivative part is an Eigen vector and GCC vectorizes it poorly here.
* No Taylor mode.

**Boost.Math** (detailed in the previous revision of this README)
* No vector mode, slow nested hyper-duals with GCC (2–9×), and Taylor mode 2.8–3× slower with GCC.

**Common to all three**
* None of them beats autodiff's `Real<N>` for the Kelvin-chain Taylor mode.
  autodiff's `forward/real` is a small, self-contained header and could be vendored or replaced by an in-house `Taylor<N>`.
* autodiff 1.1.2 does not compile with Eigen 5 (`IndexedView` mismatch in `autodiff/common/eigen.hpp`).

## Recommendation

**Sacado in vector mode** is the only tested library that beats the current implementation with GCC, on every first-, second- and third-order kernel.
It also needs no patching, and it has a long-term maintainer (Sandia; used by deal.II for exactly this purpose).

Suggested plan:

1. Introduce the backend-neutral AD types header and route all `autodiff::` uses through it, as proposed earlier.
2. Port the drivers `dF_dT`, `dF_dX`, `d2f_dT2`, `d3f_dT3` and `df_dT` to vector mode on `Sacado::Fad::SFad<…, N>`.
   `ADKernels.h` contains the reference implementations.
   The material code only sees the scalar type, so most call sites stay unchanged.
3. Keep Taylor mode (`Real<N>`) separate:
   vendor autodiff's `forward/real` header, or write a small in-house `Taylor<N>`, and benchmark it against today's.
4. Decide how to ship Sacado: vendor `packages/sacado/src` with the config stub, or depend on an installed Trilinos.
5. Re-measure in full simulations, particularly `CompressibleFiniteStrainLinearViscoelasticity` (3rd order) and the AD finite-strain materials.

Ceres is only worth reconsidering if Marmot moves to clang and accepts maintaining a patched `jet.h`.
