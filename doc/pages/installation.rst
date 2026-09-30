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

Marmot's global constants are defined ``inline const`` in the headers, not ``extern const`` in a source file,
since a Windows DLL does not export data (only functions and class members, via ``WINDOWS_EXPORT_ALL_SYMBOLS``).
New modules must follow this.
Likewise, include ``Marmot/MarmotPortability.h`` (e.g., through ``Marmot/MarmotJournal.h``) before using
``__PRETTY_FUNCTION__``, which MSVC does not provide.
The Python bindings are not yet supported on Windows.

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


