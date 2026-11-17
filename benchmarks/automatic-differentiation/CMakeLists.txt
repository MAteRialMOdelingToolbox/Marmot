# Standalone benchmark: autodiff vs. Boost.Math forward-mode AD on Marmot's AD code paths.
#
#   cmake -S benchmarks/automatic-differentiation -B build-ad-bench -DCMAKE_BUILD_TYPE=Release
#   cmake --build build-ad-bench && ./build-ad-bench/ad-backend-benchmark
#
# Only Marmot headers are used (no Marmot libraries are linked).
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

set(MARMOT_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../..)
file(GLOB MARMOT_CORE_INCLUDE_DIRS LIST_DIRECTORIES true ${MARMOT_ROOT}/modules/core/*/include)

add_executable(ad-backend-benchmark main.cpp SuiteAutodiff.cpp SuiteBoost.cpp
  ${MARMOT_ROOT}/modules/core/MarmotMathCore/src/MarmotTensor.cpp) # CommonTensors initializers
target_include_directories(ad-backend-benchmark PRIVATE ${MARMOT_CORE_INCLUDE_DIRS})
target_include_directories(ad-backend-benchmark SYSTEM PRIVATE ${FASTOR_INCLUDE_DIR} ${BOOST_MATH_INCLUDE_DIR})
target_link_libraries(ad-backend-benchmark PRIVATE Eigen3::Eigen autodiff::autodiff)
target_compile_options(ad-backend-benchmark PRIVATE -O3 -march=native)

add_executable(ad-op-benchmark opbench.cpp)
target_include_directories(ad-op-benchmark SYSTEM PRIVATE ${BOOST_MATH_INCLUDE_DIR})
target_link_libraries(ad-op-benchmark PRIVATE Eigen3::Eigen autodiff::autodiff)
target_compile_options(ad-op-benchmark PRIVATE -O3 -march=native)
