# autodiff vs. Boost.Math forward-mode AD

Evaluates replacing [autodiff](https://github.com/autodiff/autodiff) (last commit Jan 2025) with
[Boost.Math `autodiff`](https://www.boost.org/doc/libs/release/libs/math/doc/html/math_toolkit/autodiff.html)
(`boost::math::differentiation::fvar`) as Marmot's AD backend.

Both backends are wrapped in adapters with one identical vocabulary
(`HyperDual<N>`, `Taylor<N>`, `seed<level>`, `derivative<level>`, `value`, `taylorDerivatives<N>`),
so the drivers and kernels in `ADKernels.h` are shared verbatim between backends.
The kernels mirror Marmot's AD code paths and call the actual Marmot templates where possible.

## Build and run

```bash
cmake -S benchmarks/automatic-differentiation -B build-ad-bench -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_PREFIX_PATH="<eigen-3.4 prefix>;<conda env>"
cmake --build build-ad-bench
./build-ad-bench/ad-backend-benchmark [batches=15] [min-seconds-per-batch=0.05]   # kernels + correctness
./build-ad-bench/ad-op-benchmark                                                  # per-operation cost
```

Eigen 3.4 is required for the comparison, because autodiff 1.1.2 does not compile against Eigen 5 (see below).

## Results

Setup: AMD Threadripper 7970X, pinned to one core, `-O3 -march=native`, Boost 1.88, autodiff 1.1.2, Eigen 3.4.0, Fastor (conda).
Times are the median ns per kernel call. Ratio = boost / autodiff.

| kernel (Marmot path)                                                  | GCC 14.3 autodiff | GCC boost | **GCC ratio** | clang 22 autodiff | clang boost | **clang ratio** |
|-----------------------------------------------------------------------|------:|------:|------:|------:|------:|------:|
| finite-strain tangent, dual, 9 passes (`MaterialFiniteStrainAD::dF_dT`) |  1848 |  2441 | 1.32 |  1798 |  1892 | 1.05 |
| hypoelastic J2 tangent, Eigen dual, 6 passes (`MaterialHypoElasticAD`) |   492 |   541 | 1.10 |   336 |   487 | 1.45 |
| energy Hessian, dual2nd, 81 passes (`SecondOrder::d2f_dT2`)           |  4397 |  9809 | 2.23 |  5557 |  4912 | 0.88 |
| energy 3rd derivative, dual3rd, 729 passes (`ThirdOrder::d3f_dT3`)    | 90323 | 661931 | **7.33** | 147597 | 480431 | 3.26 |
| phase-field degradation, dual2nd scalar                               |    24 |   106 | 4.42 |    20 |    20 | 1.00 |
| B4 drying Kelvin moduli, Taylor 5, per stress update                  |  1315 |  4085 | 3.11 |  1176 |  1881 | 1.60 |
| power-law Kelvin moduli, Taylor 7                                     |  1136 |  3014 | 2.65 |  1274 |  1536 | 1.21 |

All derivatives agree between the backends (max relative difference ≤ 1.5e-15, ≤ 7.4e-14 for Taylor 5).
Run-to-run noise is about ±15 % on the first-order kernels.

Per-operation cost (GCC, ns per chained op, see `ad-op-benchmark`):

| type    | backend  | mul | div | sqrt | exp∘log | pow∘pow |
|---------|----------|----:|----:|-----:|--------:|--------:|
| dual    | autodiff | 3.8 | 3.8 | 10.9 | 20.2 | 43.6 |
| dual    | boost    | 4.1 | 9.8 | 12.6 | 22.3 | 52.6 |
| dual2nd | autodiff | 5.0 | 5.3 | 16.0 | 22.8 | 48.9 |
| dual2nd | boost    | 11.6 | 18.0 | 35.1 | 102.2 | 96.6 |
| dual3rd | autodiff | 7.9 | 10.5 | 20.7 | 26.8 | 53.0 |
| dual3rd | boost    | 27.6 | 32.9 | 112.2 | 446.1 | 227.6 |
| Taylor5 | autodiff | 7.4 | 10.2 | 22.5 | 25.6 | 73.8 |
| Taylor5 | boost    | 21.5 | 35.8 | 110.1 | 220.4 | 256.3 |

### Why Boost is slower for higher orders

`fvar` treats a nested `fvar<fvar<double,1>,1>` as a multivariate truncated Taylor polynomial.
Every elementary function is evaluated generically:
it computes all derivatives up to the *summed* order,
then runs a Horner scheme of full nested `fvar` multiplications with `factorial()` scaling (`apply_derivatives`).
autodiff's nested `Dual<T, G>` applies the chain rule level by level, which is much cheaper for hyper-duals.
Clang optimizes Boost's generic coefficient loops considerably better than GCC.
Raising GCC's inlining limits did not help.

## Compatibility findings

* **autodiff 1.1.2 does not compile against Eigen 5.0** (`autodiff/common/eigen.hpp`: `IndexedView` template mismatch).
  CI currently works only because it pins Eigen 3.4.0.
  This is the most concrete cost of autodiff being unmaintained.
* **Boost `fvar` and Fastor**: the Marmot energy, deformation and stress templates compile unchanged with `fvar`.
* **Boost `fvar` and Eigen**: `fvar` has an unconstrained converting constructor.
  Without help, Eigen therefore treats any matrix as a promotable scalar, which makes `A * B` ambiguous.
  This happens in both Eigen 3.4 and 5.
  It is fixed in `ADBackendBoost.h` by specializing `Eigen::internal::promote_scalar_arg` to promote only arithmetic types,
  together with `NumTraits`/`ScalarBinaryOpTraits`.
  With that, the Boost backend works with Eigen 3.4 and Eigen 5.0.
* **No mutable coefficient access in `fvar`**: seeding needs the coefficient storage.
  Its layout is identical to autodiff's nested `Dual` (2^N doubles, level-k epsilon at index 2^(N-k)).
  The adapter accesses it through one guarded `coefficients()` helper (`static_assert` on size and standard layout).
  Marmot's `increaseDualOrderWithShift` and related helpers already rely on the same kind of layout access with autodiff.
* Header parse time: autodiff (dual + real + eigen) 0.83 s, Boost.Math autodiff 1.07 s.

## Suggested changes for a Boost.Math migration

1. **Backend-neutral AD types header** (do this regardless of the backend).
   Add `Marmot/MarmotAutomaticDifferentiationTypes.h` with the adapter vocabulary from `ADBackendBoost.h`:
   `dual`, `dual2nd`, `dual3rd`, `HyperDual<N>`, `Taylor<N>`, `seed<level>`, `derivative<level>`, `value`, `taylorDerivatives<N>`.
   Then replace the `autodiff::` spellings in the 39 files that use them.
   After that, swapping the backend touches one header.
2. **Port the core helpers** (`MarmotAutomaticDifferentiation{,ForFastor}.h/.cpp`).
   Reimplement `valnode`, `increaseDualOrderWithShift`, `decreaseDualOrder*` and `shiftTo2ndOrderDual` on `coefficients()`.
   The drivers (`dF_dT`, `df_dT`, `d2f_dT2`, `d3f_dT3`, `dF_dX`, `dF_dX_2nd`) port 1:1, as shown in `ADKernels.h`.
3. **Kelvin chain**: replace `autodiff::derivatives( phi, along( 1. ), at( tau_ ) )[k]` with `taylorDerivatives<k>( phi, tau * k )[k]`,
   and change the `std::function< autodiff::Real<k> (...) >` signatures to `Taylor<k>`.
4. **Eigen glue**: move the `NumTraits` / `ScalarBinaryOpTraits` / `promote_scalar_arg` specializations into the types header.
   `Vector6dual`, `MatrixXdual` and similar become `Eigen::Matrix<dual, …>`.
5. **Utilities**: `Math::makeReal`, `Testing::getString` and `Testing::checkIfEqual` get `fvar` overloads.
   Also check the `autodiff::dual` workaround in `MarmotTensorExponential.h`.
6. **CMake**: replace the `autodiff::autodiff` dependency with Boost.Math headers.
   Boost.Math can be used standalone; the `autodiff.hpp` include set was not tested in standalone mode here.
7. **Mitigate the higher-order slowdown** (backend-independent, but needed before switching):
   * `d3f_dT3` evaluates all 729 index triples, although the third derivative is fully symmetric.
     Computing only the 165 unique triples and permuting cuts the cost about 4.4×.
     This matters because `CompressibleFiniteStrainLinearViscoelasticity` calls it in every stress update.
   * `d2f_dT2` likewise needs only 45 of 81 passes.

## Recommendation

A full switch to Boost.Math is not a clear win at the moment.
It is correct everywhere and actively maintained, and it works with Eigen 5.
But with GCC, Marmot's main compiler, it costs:

* about 1.1–1.3× on the first-order tangent paths, which are the most common AD use;
* 2–4× on `dual2nd` and Taylor mode, including B4's per-update Kelvin-chain moduli;
* about 7× on `dual3rd`.

With clang it is close to parity, except for `dual3rd` (about 3×).

The suggested order:

* Do steps 1–5 (the backend-neutral layer) now; they are useful whatever backend Marmot ends up with.
* Keep autodiff behind that layer for now, vendored, with a patched `eigen.hpp` for Eigen 5.
* Use step 7 to remove most of the higher-order cost.
* Then re-run this benchmark through the layer. Boost becomes the drop-in fallback if autodiff breaks again,
  and the right choice if the first-order overhead turns out to be acceptable in full simulations.
