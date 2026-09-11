# Prefer active virtualenv / Conda environment
set(Python_FIND_VIRTUALENV FIRST)
set(Python_FIND_STRATEGY LOCATION)

# Find Python first, so that the very same interpreter is used to patch nanobind
find_package(Python 3.8 COMPONENTS Interpreter Development.Module REQUIRED)

# Fetch nanobind
include(FetchContent)
FetchContent_Declare(
  nanobind
  GIT_REPOSITORY https://github.com/wjakob/nanobind.git
  GIT_TAG        v1.9.2
  GIT_SHALLOW    TRUE
  EXCLUDE_FROM_ALL TRUE
  PATCH_COMMAND "${Python_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/patch_nanobind.py" "include/nanobind/nb_attr.h"
)
FetchContent_MakeAvailable(nanobind)

# Define the python extension module
nanobind_add_module(_marmot
    src/bind_main.cpp
    src/bind_solver_finite_strain.cpp
    src/bind_solver_hypo_elastic.cpp
)

# Link against the Marmot C++ library
target_link_libraries(_marmot PRIVATE Marmot)

# Set RPATH and output directory for in-tree execution
set_target_properties(_marmot PROPERTIES
    BUILD_RPATH "${CMAKE_BINARY_DIR}/lib"
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/python/marmot"
)

# Copy python package files to build directory so tests can run in-place
file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/marmot/"
     DESTINATION "${CMAKE_BINARY_DIR}/python/marmot"
     FILES_MATCHING PATTERN "*.py"
)

# Add include directory for fastor_nanobind_caster.h
target_include_directories(_marmot PRIVATE src)

execute_process(
    COMMAND "${Python_EXECUTABLE}" -c "import sysconfig, os; print(os.path.relpath(sysconfig.get_path('platlib'), sysconfig.get_path('data')))"
    OUTPUT_VARIABLE PYTHON_RELATIVE_INSTALL_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

install(TARGETS _marmot DESTINATION ${PYTHON_RELATIVE_INSTALL_DIR}/marmot)
install(DIRECTORY marmot/ DESTINATION ${PYTHON_RELATIVE_INSTALL_DIR}/marmot)

# Register Python example test suite with CTest
file(GLOB_RECURSE PYTHON_TEST_SCRIPTS CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/examples/python/*.py")
foreach(SCRIPT ${PYTHON_TEST_SCRIPTS})
    get_filename_component(SCRIPT_NAME ${SCRIPT} NAME_WE)
    get_filename_component(DIR_PATH ${SCRIPT} DIRECTORY)
    get_filename_component(MAT_NAME ${DIR_PATH} NAME)
    add_test(NAME "Python_${MAT_NAME}_${SCRIPT_NAME}"
             COMMAND "${Python_EXECUTABLE}" "${SCRIPT}")
    set_tests_properties("Python_${MAT_NAME}_${SCRIPT_NAME}" PROPERTIES
             ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python")
endforeach()
