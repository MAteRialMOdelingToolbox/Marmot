# Standalone benchmark: forward-mode AD backends on Marmot's AD code paths.
#
#   cmake -S benchmarks/automatic-differentiation -B build-ad-bench -DCMAKE_BUILD_TYPE=Release \
#         -DSACADO_SOURCE_DIR=<Trilinos>/packages/sacado -DCERES_SOURCE_DIR=<ceres-solver 2.2.0>
#   cmake --build build-ad-bench && ./build-ad-bench/ad-backend-benchmark
#
# Only Marmot headers are used (plus MarmotTensor.cpp for the CommonTensors initializers).
# Sacado and Ceres are used header-only straight from their source trees; no build/install of either is needed.
cmake_minimum_required(VERSION 3.20)
project(MarmotADBackendBenchmark CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

find_package(Eigen3 3.4 REQUIRED NO_MODULE)
find_package(autodiff REQUIRED)
find_path(FASTOR_INCLUDE_DIR Fastor/Fastor.h REQUIRED)
find_path(BOOST_MATH_INCLUDE_DIR boost/math/differentiation/autodiff.hpp REQUIRED)

set(SACADO_SOURCE_DIR "" CACHE PATH "Trilinos packages/sacado directory")
set(CERES_SOURCE_DIR "" CACHE PATH "ceres-solver 2.2.0 source directory")
if(NOT EXISTS "${SACADO_SOURCE_DIR}/src/Sacado.hpp")
  message(FATAL_ERROR "Set SACADO_SOURCE_DIR to Trilinos' packages/sacado")
endif()
if(NOT EXISTS "${CERES_SOURCE_DIR}/include/ceres/jet.h")
  message(FATAL_ERROR "Set CERES_SOURCE_DIR to a ceres-solver 2.2.0 checkout")
endif()

# ceres::Jet's value constructor is explicit, which Fastor's `T x = 0.` cannot work with. Apply the minimal change
# documented in ceres-2.2.0-jet-fastor.patch to a copy of the headers in the build tree.
set(CERES_PATCHED_INCLUDE_DIR ${CMAKE_BINARY_DIR}/ceres-patched/include)
file(COPY ${CERES_SOURCE_DIR}/include/ceres DESTINATION ${CERES_PATCHED_INCLUDE_DIR})
file(READ ${CERES_SOURCE_DIR}/include/ceres/jet.h JET_H)
set(JET_CTOR_UPSTREAM "  explicit Jet(const T& value) {\n    a = value;\n    v.setConstant(Scalar());\n  }\n")
string(FIND "${JET_H}" "${JET_CTOR_UPSTREAM}" JET_CTOR_POS)
if(JET_CTOR_POS EQUAL -1)
  message(FATAL_ERROR "ceres/jet.h does not match ceres-2.2.0-jet-fastor.patch (expected ceres-solver 2.2.0)")
endif()
string(REPLACE "${JET_CTOR_UPSTREAM}"
  "  Jet(const T& value) {\n    a = value;\n    v.setConstant(Scalar());\n  }\n\n  template <typename S,\n            std::enable_if_t<std::is_arithmetic_v<S> && !std::is_same_v<S, T>, int> = 0>\n  Jet(S value) : Jet(T(value)) {}\n"
  JET_H "${JET_H}")
file(WRITE ${CERES_PATCHED_INCLUDE_DIR}/ceres/jet.h "${JET_H}")

set(MARMOT_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../..)
file(GLOB MARMOT_CORE_INCLUDE_DIRS LIST_DIRECTORIES true ${MARMOT_ROOT}/modules/core/*/include)

set(AD_BACKEND_INCLUDE_DIRS
    ${FASTOR_INCLUDE_DIR}
    ${BOOST_MATH_INCLUDE_DIR}
    ${SACADO_SOURCE_DIR}/src
    ${SACADO_SOURCE_DIR}/src/new_design
    ${SACADO_SOURCE_DIR}/src/mpl
    ${CERES_PATCHED_INCLUDE_DIR})

add_executable(ad-backend-benchmark main.cpp SuiteAutodiff.cpp SuiteBoost.cpp SuiteSacado.cpp SuiteCeres.cpp
  ${MARMOT_ROOT}/modules/core/MarmotMathCore/src/MarmotTensor.cpp) # CommonTensors initializers
target_include_directories(ad-backend-benchmark PRIVATE ${MARMOT_CORE_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/sacado-config)
target_include_directories(ad-backend-benchmark SYSTEM PRIVATE ${AD_BACKEND_INCLUDE_DIRS})
target_link_libraries(ad-backend-benchmark PRIVATE Eigen3::Eigen autodiff::autodiff)
target_compile_options(ad-backend-benchmark PRIVATE -O3 -march=native)

add_executable(ad-op-benchmark opbench.cpp)
target_include_directories(ad-op-benchmark PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/sacado-config)
target_include_directories(ad-op-benchmark SYSTEM PRIVATE ${AD_BACKEND_INCLUDE_DIRS})
target_link_libraries(ad-op-benchmark PRIVATE Eigen3::Eigen autodiff::autodiff)
target_compile_options(ad-op-benchmark PRIVATE -O3 -march=native)
