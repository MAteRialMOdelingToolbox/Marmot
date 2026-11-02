# The consumer test links the shared library only, as EdelweissFE and the Abaqus interface do, and so checks that
# the exported API (MARMOT_API, see MarmotPortability.h) is complete. All other tests link Marmot's object files and
# cannot detect a missing export; see add_marmot_test in the top-level CMakeLists.txt.
add_executable(TestExportedAPI TestExportedAPI.cpp)
target_link_libraries(TestExportedAPI PRIVATE ${PROJECT_NAME})
add_test(NAME TestExportedAPI COMMAND TestExportedAPI WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
if(MARMOT_ENABLE_COVERAGE)
  target_compile_options(TestExportedAPI PRIVATE --coverage -O0)
  target_link_options(TestExportedAPI PRIVATE --coverage)
endif()
