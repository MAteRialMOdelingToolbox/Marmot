Installation
============

Requirements
************

Marmot itself requires the `Eigen <https://eigen.tuxfamily.org/>`_ library,
`autodiff <https://github.com/autodiff/autodiff>`_,
and `Fastor <https://github.com/romeric/Fastor>`_.

These are header-only libraries, so no compilation is required.

Both Eigen 3.4 and Eigen 5 are supported.
Eigen 5 requires an autodiff version with Eigen 5 support;
autodiff 1.1.2 and older do not compile against Eigen 5 (see `autodiff#397 <https://github.com/autodiff/autodiff/pull/397>`_).

Building with Anaconda
**********************

Building with anaconda is the easiest way to get a working version of Marmot.

Assuming that you are in an empty directory,
you can quickly get a working version of Marmot in a Linux based
environment:

Installation steps
__________________

If necessary, get Miniforge:

.. code-block:: console
   :caption: Step 1

    curl -L -O \
        "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
    bash "Miniforge3-$(uname)-$(uname -m).sh" -b -p ./miniforge3

Add conda to your environment:

.. code-block:: console
   :caption: Step 2

    export MARMOTROOT=$PWD
    export PATH=$MARMOTROOT/miniforge3/bin:$PATH
    conda init --all
    exit

Restart shell and activate conda

.. code-block:: console
   :caption: Step 3

    export MARMOTROOT=$PWD
    conda activate
    mamba install cmake make compilers

Get Eigen:

.. code-block:: console
   :caption: Step 4

    cd $MARMOTROOT
    git clone --branch 3.4.0  https://gitlab.com/libeigen/eigen.git
    cd eigen
    mkdir build
    cd build
    cmake \
        -DBUILD_TESTING=OFF  \
        -DINCLUDE_INSTALL_DIR=$CONDA_PREFIX/include \
        -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX \
        ..
    make install

Get autodiff:

.. code-block:: console
   :caption: Step 5

    cd $MARMOTROOT
    git clone --branch v1.1.0 https://github.com/autodiff/autodiff.git
    cd autodiff
    mkdir build
    cd build
    cmake \
        -DAUTODIFF_BUILD_TESTS=OFF \
        -DAUTODIFF_BUILD_PYTHON=OFF \
        -DAUTODIFF_BUILD_EXAMPLES=OFF \
        -DAUTODIFF_BUILD_DOCS=OFF \
        -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX \
        ..
    make install

Get Fastor:

.. code-block:: console
   :caption: Step 6

    cd $MARMOTROOT
    git clone https://github.com/romeric/Fastor.git
    cd Fastor
    cmake -DBUILD_TESTING=OFF -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX .
    make install
    cd ../

Get Marmot:

.. code-block:: console
   :caption: Step 7

    cd $MARMOTROOT
    git clone https://github.com/MAteRialMOdelingToolbox/Marmot.git
    cd Marmot
    mkdir build
    cd build
    cmake \
        -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX \
        -DMARMOT_BUILD_PYTHON_BINDINGS=ON \
        ..
    make install
    ctest --output-on-failure

Build options
*************

The following CMake options adjust the build:

* ``-DBUILD_TESTING=OFF`` skips the test executables, e.g., for install-only builds.
* ``-DMARMOT_MARCH_NATIVE=ON`` compiles for the host CPU (``-march=native``, GCC/Clang only), which lets Fastor
  vectorize with AVX2/FMA. The resulting library does not run on older CPUs.
* ``-DMARMOT_ENABLE_COVERAGE=ON`` instruments the library and the tests for ``gcov`` (with ``-DCMAKE_BUILD_TYPE=Debug``).
* ``-DMARMOT_BUILD_PYTHON_BINDINGS=ON`` builds the Python bindings.

Modules
*******

Marmot consists of modules, one directory ``modules/<category>/<Name>/`` each, with the categories
``core``, ``materials``, ``elements``, ``particles``, ``materialpoints``, ``cells`` and ``cellelements``.
Every such directory containing a ``module.cmake`` is found automatically; a module of your own is added by placing
its directory there. The variables ``CORE_MODULES``, ``MATERIAL_MODULES``, ``ELEMENT_MODULES``, ``PARTICLE_MODULES``,
``MATERIALPOINT_MODULES``, ``CELL_MODULES`` and ``CELLELEMENT_MODULES`` (default ``all``) select a subset, e.g.,
``-DMATERIAL_MODULES="LinearElastic;VonMises"``.

A ``module.cmake`` declares the module and the modules whose headers it includes:

.. code-block:: cmake

    marmot_add_module(MyMaterial
        REQUIRES MarmotFiniteStrainMechanicsCore)

A module whose required module is not built is skipped with a warning; if it was selected explicitly, configuring
fails. A module sees only the headers of the modules it requires, so a missing ``REQUIRES`` shows as a compile error.
All modules are compiled into the one library ``libMarmot``.

Installing copies the headers of every built module to ``<prefix>/include/Marmot``.
Other CMake projects use the installed library through ``find_package``:

.. code-block:: cmake

    find_package(Marmot REQUIRED)
    target_link_libraries(<target> PRIVATE Marmot::Marmot)

``Marmot::Marmot`` brings along Eigen, autodiff and Fastor; autodiff and Fastor are found either through their
CMake packages or, if they were installed as plain headers, through a header search.

Building on Windows
*******************

Marmot builds as a DLL with MSVC (Visual Studio 2022), in the ``Release`` configuration.
Install Eigen, autodiff and Fastor into a common prefix as above (``cmake --install`` instead of ``make install``),
then build Marmot from a *Developer PowerShell for VS 2022*:

.. code-block:: console

    cmake -S Marmot -B Marmot/build -DCMAKE_PREFIX_PATH=<prefix> -DCMAKE_INSTALL_PREFIX=<prefix>
    cmake --build Marmot/build --config Release --parallel
    ctest --test-dir Marmot/build -C Release --output-on-failure
    cmake --install Marmot/build --config Release

This installs ``Marmot.dll`` into ``<prefix>/bin`` and its import library ``Marmot.lib`` into ``<prefix>/lib``.
Programs linking Marmot must find ``Marmot.dll`` at run time, e.g., through ``PATH``.

A Windows DLL exports only what is marked for export, which in Marmot is ``MARMOT_API``
(defined in ``Marmot/MarmotPortability.h``).
The rule for what is marked: the exported API is the interface layer through which a consumer drives Marmot,
that is, the element and material factories, the interface classes they hand out (``MarmotElement``,
``MarmotMaterialSection``, ``ElementProperties``), and the few non-virtual functions a consumer calls on them,
marked at the smallest granularity that links (a single member rather than its class, where the class is otherwise
reached through virtual functions only).
Everything else is used through the virtual functions of the objects the factories create.
Exporting everything instead is not an option: it exceeds the limit of 65535 exported symbols of a Windows DLL,
since that would include every Eigen, Fastor and autodiff template instantiated in Marmot.
Code that needs more of Marmot must mark it ``MARMOT_API``, following the rule above.
To check this without Windows, configure with ``-DMARMOT_EXPORT_API_ONLY=ON``,
which exports only the ``MARMOT_API`` symbols on Linux and macOS, too; the Ubuntu CI builds this way as well.
The module tests link Marmot's object files directly and are not affected.
The one test that is, ``TestExportedAPI`` in ``tests/consumer``, links the shared library only and uses Marmot
as a consumer does; it is the test that fails when the exported API is incomplete.

Marmot's global constants are defined ``inline const`` in the headers, not ``extern const`` in a source file:
exported data would have to be marked ``MARMOT_API`` as well, and could not be used in constant expressions
across the DLL boundary.
An ``inline`` variable is instantiated in every translation unit that includes its header, so everything its
initializer calls must be defined in a header, too (or be marked ``MARMOT_API``). New modules must follow this.
Likewise, include ``Marmot/MarmotPortability.h`` (e.g., through ``Marmot/MarmotJournal.h``) before using
``__PRETTY_FUNCTION__``, which MSVC does not provide.
The Python bindings are not yet supported on Windows, nor with ``MARMOT_EXPORT_API_ONLY``;
configuring with ``-DMARMOT_BUILD_PYTHON_BINDINGS=ON`` fails there.

Building with Python Bindings
*****************************

Marmot optionally provides a Python interface using `nanobind <https://github.com/wjakob/nanobind>`_.
To enable the Python module, configure CMake with `-DMARMOT_BUILD_PYTHON_BINDINGS=ON`:

.. code-block:: console

    cmake -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX -DMARMOT_BUILD_PYTHON_BINDINGS=ON ..
    make install
    ctest --output-on-failure

After installation, the ``marmot`` package is available directly in Python:

.. code-block:: python

    import marmot
    import numpy as np

    props = np.array([20000.0, 0.25])
    solver = marmot.solvers.HypoElasticSolver("LINEARELASTIC", props)


