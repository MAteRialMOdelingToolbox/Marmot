# Marmot's tests/consumer/TestExportedAPI.cpp, built as an external project against the installed package.
cmake_minimum_required(VERSION 3.16)
project(MarmotConsumer LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(MSVC)
  add_compile_options(/bigobj /permissive- /Zc:__cplusplus)
endif()

find_package(Marmot REQUIRED CONFIG)

add_executable(TestExportedAPI ../tests/consumer/TestExportedAPI.cpp)
target_link_libraries(TestExportedAPI PRIVATE Marmot::Marmot)
