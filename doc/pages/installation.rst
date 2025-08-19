Installation
============

Requirements
************

Marmot itself requires the `Eigen <https://eigen.tuxfamily.org/>`_ library,
autodiff `autodiff <github.com/autodiff/autodiff>`_,
and potentially `Fastor <https://github.com/romeric/Fastor>`_, depending on the requested modules.

These are header-only libraries, so no compilation is required.

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
        https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-aarch64.sh
    bash Miniforge3-Linux-aarch64.sh -b -p ./miniforge3

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
        ..
    make install
    ctest --output-on-failure 

