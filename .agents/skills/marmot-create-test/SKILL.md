---
name: marmot-create-test
description: >-
  Step-by-step instructions for constructing, validating, and registering unit and regression tests in Marmot.
  Use when adding new test executables for bug fixes, new features, or material/element verifications.
---

# Creating Tests in Marmot

All unit and regression tests are standalone C++ executables registered in CMake via `test.cmake` and driven by CTest.

---

## 1. Testing Rules
1. **Analytical Validation**: Validate constitutive responses and element modes against closed-form analytical solutions.
2. **Tangent Verification**: Validate consistent tangents against `MarmotMathCore` numerical differentiation (`MarmotNumericalDifferentiation.h`) or automatic differentiation (`autodiff` / `MarmotAutomaticDifferentiation.h`).
3. **Complex-Step Method**: Use `Marmot::NumericalAlgorithms::Differentiation::Complex::forwardDifference` for machine-precision derivatives ($h \sim 10^{-20}$) without subtractive cancellation.
4. **Performance & Hygiene**: Tests must run fast (< 0.1s), remain deterministic, and update `stateVars` in place.

---

## 2. Test Skeleton with Analytical & Tangent Checks (`test/test.cpp`)

```cpp
#include "Marmot/MarmotNumericalDifferentiation.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MyMaterial.h"
#include <iostream>
#include <vector>

using namespace Marmot;
using namespace Marmot::Testing;
namespace NumDiff = Marmot::NumericalAlgorithms::Differentiation;

int testMaterial() {
  std::vector< double > props{ 210000.0, 0.3, 7.85e-9 }; // E, nu, rho
  Materials::MyMaterial mat( props.data(), static_cast< int >( props.size() ), 1 );

  MarmotMaterialHypoElastic::state3D state;
  double stateVars[10] = { 0.0 };
  state.stateVars = stateVars;

  Vector6d dStrain = Vector6d::Zero();
  dStrain( 0 ) = 1e-3;

  Matrix6d tangent = Matrix6d::Zero();
  MarmotMaterialHypoElastic::timeInfo tInfo{ 0.0, 1.0 };
  mat.computeStress( state, tangent, dStrain, tInfo );

  // 1. Check against analytical solution: sigma_11 = 210.0 MPa
  if ( !checkIfEqual( state.stress( 0 ), 210.0, 1e-10 ) ) {
    std::cerr << "Analytical stress mismatch: " << state.stress( 0 ) << std::endl;
    return 1;
  }

  // 2. Check tangent against MarmotMathCore numerical differentiation
  NumDiff::vector_to_vector_function_type stressFunc = [&]( const VectorXd& dEps ) -> VectorXd {
    MarmotMaterialHypoElastic::state3D tmpState;
    double tmpVars[10] = { 0.0 };
    std::copy( stateVars, stateVars + 10, tmpVars );
    tmpState.stateVars = tmpVars;
    Matrix6d tmpTangent = Matrix6d::Zero();
    mat.computeStress( tmpState, tmpTangent, dEps, tInfo );
    return tmpState.stress;
  };

  MatrixXd numTangent = NumDiff::centralDifference( stressFunc, dStrain );
  if ( !checkIfEqual( tangent, numTangent, 1e-6 ) ) {
    std::cerr << "Tangent verification failed against numerical differentiation!" << std::endl;
    return 1;
  }

  return 0;
}

int main() {
  int res = testMaterial();
  if ( res == 0 ) std::cout << "All tests passed." << std::endl;
  return res;
}
```

---

## 3. CMake Registration & Execution

`test.cmake`:
```cmake
SET(CURR_TEST_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/test")
add_marmot_test("TestMarmotMyModule" "${CURR_TEST_SOURCE_DIR}/test.cpp")
```

Execution:
```bash
cd build && make -j$(nproc)
ctest -R TestMarmotMyModule --output-on-failure
ctest --output-on-failure
```
