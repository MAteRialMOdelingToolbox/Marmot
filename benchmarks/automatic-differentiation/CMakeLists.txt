# Standalone benchmark: automatic differentiation backends on Marmot's AD code paths.
#
#   cmake -S benchmarks/automatic-differentiation -B build-ad-bench -DCMAKE_BUILD_TYPE=Release \
#         -DAUTODIFF_INCLUDE_DIR=<autodiff 1.1.2 + autodiff#397> -DFASTOR_INCLUDE_DIR=<dir containing Fastor/> \
#         -DSACADO_SOURCE_DIR=<Trilinos>/packages/sacado -DCERES_SOURCE_DIR=<ceres-solver 2.2.0> \
#         -DSTAN_SOURCE_DIR=<stan math v5.4.0> -DTINYAD_SOURCE_DIR=<TinyAD> -DADEPT_SOURCE_DIR=<Adept-2 v2.1.3>
#   cmake --build build-ad-bench && ./build-ad-bench/ad-backend-benchmark
#
# Only Marmot headers are used (plus MarmotTensor.cpp for the CommonTensors initializers). All backends are used
# straight from their source trees; none needs to be installed. Boost.Math comes from Stan's bundled Boost, so that a
# single Boost version ends up in the executable.
cmake_minimum_required(VERSION 3.20)
project(MarmotADBackendBenchmark CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

find_package(Eigen3 REQUIRED NO_MODULE)
if(Eigen3_VERSION VERSION_LESS 5.0)
  message(FATAL_ERROR "Marmot requires Eigen 5, found ${Eigen3_VERSION}")
endif()

set(AUTODIFF_INCLUDE_DIR "" CACHE PATH "autodiff 1.1.2 with the Eigen 5 patch autodiff#397 (directory containing autodiff/)")
set(FASTOR_INCLUDE_DIR "" CACHE PATH "directory containing Fastor/ (and no other autodiff or Boost)")
set(SACADO_SOURCE_DIR "" CACHE PATH "Trilinos packages/sacado directory")
set(CERES_SOURCE_DIR "" CACHE PATH "ceres-solver 2.2.0 source directory")
set(STAN_SOURCE_DIR "" CACHE PATH "Stan Math v5.4.0 source directory")
set(TINYAD_SOURCE_DIR "" CACHE PATH "TinyAD source directory")
set(ADEPT_SOURCE_DIR "" CACHE PATH "Adept-2 v2.1.3 source directory")
foreach(_check
    "AUTODIFF_INCLUDE_DIR;autodiff/forward/dual.hpp"
    "FASTOR_INCLUDE_DIR;Fastor/Fastor.h"
    "SACADO_SOURCE_DIR;src/Sacado.hpp"
    "CERES_SOURCE_DIR;include/ceres/jet.h"
    "STAN_SOURCE_DIR;stan/math/fwd.hpp"
    "TINYAD_SOURCE_DIR;include/TinyAD/Scalar.hh"
    "ADEPT_SOURCE_DIR;include/adept.h")
  list(GET _check 0 _var)
  list(GET _check 1 _file)
  if(NOT EXISTS "${${_var}}/${_file}")
    message(FATAL_ERROR "Set ${_var} (${${_var}}/${_file} not found)")
  endif()
endforeach()

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

# Adept's run-time part (stack, Jacobian sweeps); BLAS/LAPACK-dependent parts compile to nothing without HAVE_BLAS
file(GLOB ADEPT_SOURCES ${ADEPT_SOURCE_DIR}/adept/*.cpp)
add_library(adept STATIC ${ADEPT_SOURCES})
target_include_directories(adept SYSTEM PUBLIC ${ADEPT_SOURCE_DIR}/include)
target_compile_options(adept PRIVATE -O3 -march=native -w)

set(MARMOT_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../..)
file(GLOB MARMOT_CORE_INCLUDE_DIRS LIST_DIRECTORIES true ${MARMOT_ROOT}/modules/core/*/include)

set(AD_BACKEND_INCLUDE_DIRS
    ${AUTODIFF_INCLUDE_DIR}
    ${FASTOR_INCLUDE_DIR}
    ${STAN_SOURCE_DIR}
    ${STAN_SOURCE_DIR}/lib/boost_1.87.0
    ${STAN_SOURCE_DIR}/lib/tbb_2020.3/include
    ${STAN_SOURCE_DIR}/lib/sundials_6.1.1/include
    ${SACADO_SOURCE_DIR}/src
    ${SACADO_SOURCE_DIR}/src/new_design
    ${SACADO_SOURCE_DIR}/src/mpl
    ${CERES_PATCHED_INCLUDE_DIR}
    ${TINYAD_SOURCE_DIR}/include)

foreach(_target ad-backend-benchmark ad-op-benchmark adept-fastor-einsum-repro)
  if(_target STREQUAL "ad-backend-benchmark")
    add_executable(${_target} main.cpp SuiteAutodiff.cpp SuiteBoost.cpp SuiteSacado.cpp SuiteCeres.cpp SuiteStan.cpp
      SuiteTinyAD.cpp SuiteAdept.cpp
      ${MARMOT_ROOT}/modules/core/MarmotMathCore/src/MarmotTensor.cpp) # CommonTensors initializers
  elseif(_target STREQUAL "ad-op-benchmark")
    add_executable(${_target} opbench.cpp)
  else()
    add_executable(${_target} repro/adept-fastor-einsum.cpp ${MARMOT_ROOT}/modules/core/MarmotMathCore/src/MarmotTensor.cpp)
  endif()
  target_include_directories(${_target} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${MARMOT_CORE_INCLUDE_DIRS}
                                                ${CMAKE_CURRENT_SOURCE_DIR}/sacado-config)
  target_include_directories(${_target} SYSTEM PRIVATE ${AD_BACKEND_INCLUDE_DIRS})
  target_compile_definitions(${_target} PRIVATE BOOST_DISABLE_ASSERTS _REENTRANT)
  target_link_libraries(${_target} PRIVATE Eigen3::Eigen adept)
  target_compile_options(${_target} PRIVATE -O3 -march=native)
endforeach()
