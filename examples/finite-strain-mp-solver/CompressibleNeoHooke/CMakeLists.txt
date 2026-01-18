cmake_minimum_required(VERSION 3.16)
project(MarmotSolverExample LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Python3 COMPONENTS Interpreter REQUIRED)

execute_process(
    COMMAND "${Python3_EXECUTABLE}" -c "import sys; print(sys.prefix, end='')"
    OUTPUT_VARIABLE PYTHON_PREFIX
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

set(MARMOT_ROOT "${PYTHON_PREFIX}")
set(EIGEN_DIR   "${PYTHON_PREFIX}/include/eigen3")
set(FASTOR_DIR  "${PYTHON_PREFIX}/include/Fastor")

add_executable(solver_example example.cpp)

target_include_directories(solver_example PRIVATE
    "${MARMOT_ROOT}/include"
    "${EIGEN_DIR}"
    "${FASTOR_DIR}"
)

target_link_directories(solver_example PRIVATE "${MARMOT_ROOT}/lib")

target_link_libraries(solver_example PRIVATE Marmot)

target_compile_options(solver_example PRIVATE -Wall -Wextra -fPIC)

set_target_properties(solver_example PROPERTIES
    BUILD_RPATH   "${MARMOT_ROOT}/lib"
    INSTALL_RPATH "${MARMOT_ROOT}/lib"
)
