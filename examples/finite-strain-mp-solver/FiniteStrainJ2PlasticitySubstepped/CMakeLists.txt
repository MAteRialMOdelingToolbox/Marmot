cmake_minimum_required(VERSION 3.16)
project(MarmotSolverExample LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Marmot and its dependencies (Eigen, autodiff, Fastor) are resolved through
# find_package. If Marmot is installed into a non-default prefix, point
# CMAKE_PREFIX_PATH at it; an activated conda environment is searched
# automatically (CMake derives search prefixes from PATH).
find_package(Marmot REQUIRED)

add_executable(solver_example example.cpp)
target_link_libraries(solver_example PRIVATE Marmot::Marmot)
target_compile_options(solver_example PRIVATE -Wall -Wextra)
