# Installation

## Third-party dependencies

Marmot requires the Eigen (>3.3.8) and the autodiff (>0.6.0) libraries:

```bash
git clone https://gitlab.com/libeigen/eigen.git
cd eigen
mkdir build
cd build
cmake ..
sudo make install
```

```bash
git clone https://github.com/autodiff/autodiff.git
cd autodiff
mkdir build
cd build
cmake -DAUTODIFF_BUILD_TESTS=OFF -DAUTODIFF_BUILD_PYTHON=OFF -DAUTODIFF_BUILD_EXAMPLES=OFF -DAUTODIFF_BUILD_DOCS=OFF ..
sudo make install
```

For anaconda users, those libraries can also be installed from the channels.

## How to install Marmot

Marmot including all submodules can be installed with the following steps:

```bash
git clone --recurse-submodules https://github.com/MAteRialMOdelingToolbox/Marmot/
cd Marmot
mkdir build
cd build
cmake ..
make
sudo make install
```

CMake options CORE_MODULES, ELEMENT_MODULES and MATERIAL_MODULES
allow to specify the modules which should be compiled, either by passing a semicolon seperated list,
option **none** or option **all** (default).
For instance:

```bash
git clone --recurse-submodules https://github.com/MAteRialMOdelingToolbox/Marmot/
cd Marmot
mkdir build
cd build
cmake \
-DCORE_MODULES='MarmotMechanicsCore;MarmotFiniteElementCore' \
-DELEMENT_MODULES='none' \
-DMATERIAL_MODULES='all' \
..
make
sudo make install
```

CMake option CMAKE_INSTALL_PREFIX allows to specify the installation directory.
For instance:

```bash
git clone --recurse-submodules https://github.com/MAteRialMOdelingToolbox/Marmot/
cd Marmot
mkdir build
cd build
cmake \
-DCMAKE_INSTALL_PREFIX=/your/special/installationdirectory \
..
make
make install
```
