[![build (Ubuntu)](https://github.com/MAteRialMOdelingToolbox/Marmot/actions/workflows/build_ubuntu.yml/badge.svg?branch=master)](https://github.com/MAteRialMOdelingToolbox/Marmot/actions/workflows/build_ubuntu.yml)
[![build (macOS)](https://github.com/MAteRialMOdelingToolbox/Marmot/actions/workflows/build_macos.yml/badge.svg?branch=master)](https://github.com/MAteRialMOdelingToolbox/Marmot/actions/workflows/build_macos.yml)
[![build (Windows)](https://github.com/MAteRialMOdelingToolbox/Marmot/actions/workflows/build_windows.yml/badge.svg?branch=master)](https://github.com/MAteRialMOdelingToolbox/Marmot/actions/workflows/build_windows.yml)
[![conda](https://github.com/MAteRialMOdelingToolbox/Marmot/actions/workflows/conda.yml/badge.svg?branch=master)](https://github.com/MAteRialMOdelingToolbox/Marmot/actions/workflows/conda.yml)
[![codecov](https://codecov.io/gh/MAteRialMOdelingToolbox/Marmot/branch/master/graph/badge.svg)](https://codecov.io/gh/MAteRialMOdelingToolbox/Marmot)
[![clang-format](https://github.com/MAteRialMOdelingToolbox/Marmot/actions/workflows/indent.yml/badge.svg?branch=master)](https://github.com/MAteRialMOdelingToolbox/Marmot/actions/workflows/indent.yml)
[![documentation](https://github.com/MAteRialMOdelingToolbox/Marmot/actions/workflows/sphinx.yml/badge.svg?branch=master)](https://materialmodelingtoolbox.github.io/Marmot/)
[![license](https://img.shields.io/badge/license-LGPL--2.1--or--later-blue.svg)](LICENSE.md)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17509327.svg)](https://doi.org/10.5281/zenodo.17509327)

# Marmot

`Marmot` (MAteRialMOdellingToolbox) is a C++-library aiming to provide robust and efficient implementations of finite elements
and state-of-the-art constitutive models for different materials.
In particular for quasi-brittle materials such as (sprayed) concrete, rock and soils.
It uses modern, object-oriented programming techniques and provides a generic interface that can be easily wrapped into different numerical solvers.
It can be  seamlessly used with the [EdelweissFE](https://github.com/EdelweissFE/EdelweissFE) finite element code, but also with commercial and open source finite element codes.
Standard interfaces for commercial finite element codes, such as Abaqus, Plaxis and open source codes like [MOOSE](https://github.com/idaholab/moose) or OpenSees, are ready-to-use available.

![Truss in compression using a micropolar von Mises plasticity model](share/truss.gif)

Truss in compression using a micropolar von Mises plasticity model.

## Quick start

Prebuilt packages of released versions (Linux, macOS, Windows):

```bash
conda create -n marmot -c https://repo.prefix.dev/matthiasneuner/edelweiss -c conda-forge marmot
```

To develop Marmot or add modules, build it yourself as described below, in a separate environment.

**Upgrading from v25.11?** Release v26.11 reorganized all namespaces (e.g., `MarmotLibrary::MarmotElementFactory` is
now `Marmot::Registry::ElementFactory`) without compatibility aliases. The
[changelog](CHANGELOG.md) lists every renamed class and namespace.

`Marmot` requires the
[Eigen](https://gitlab.com/libeigen/eigen) (>3.3.8, including Eigen 5),
[autodiff](https://github.com/autodiff/autodiff) (>0.6.0)
and [Fastor](https://github.com/romeric/Fastor) (≥0.6.4) libraries to be installed.
Detailed instructions on how to install these dependencies can be found [here](https://materialmodelingtoolbox.github.io/Marmot/pages/installation.html).

For a quick start we recommend to use conda.
A conda environment including all needed dependencies can be created by running

```bash
conda create -n marmot-env -c https://repo.prefix.dev/matthiasneuner/edelweiss -c conda-forge cmake compilers eigen autodiff fastor
```
To activate the environment use

```bash
conda activate marmot-env
```

Alternatively, Fastor can be installed manually by executing
```bash
git clone --branch V0.6.4 https://github.com/romeric/Fastor.git
cd Fastor
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX ..
make install
```

To finally build and install `Marmot` execute
```bash
git clone https://github.com/MAteRialMOdelingToolbox/Marmot.git
cd Marmot
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX ..
make install
```

Within this build directory, the tests can be executed by running

```bash
ctest --output-on-failure
```

On Windows, Marmot builds as a DLL with MSVC (Release configuration), see the
[installation instructions](https://materialmodelingtoolbox.github.io/Marmot/pages/installation.html#building-on-windows).

### Python bindings

`Marmot` optionally provides a Python interface to the material point solvers. It requires the
Python development headers and fetches [nanobind](https://github.com/wjakob/nanobind) at configure
time, and is therefore disabled by default. To build and install it as well, configure with

```bash
cmake -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX -DMARMOT_BUILD_PYTHON_BINDINGS=ON ..
make install
```

The examples in `examples/python` can then be run directly, e.g.

```bash
python examples/python/LinearElastic/example_strain_controlled.py
```

Each example compares its results against a stored reference solution. Pass
`--writeReferenceSolution` to regenerate that reference instead of checking against it.

## How to use Marmot with EdelweissFE
The [EdelweissFE](https://github.com/Edelweiss-Numerics/EdelweissFE) finite element code is designed to work seamlessly with `Marmot`.
A released EdelweissFE including `Marmot` can be installed with

```bash
conda create -n edelweissfe -c https://repo.prefix.dev/matthiasneuner/edelweiss -c conda-forge edelweissfe
```

To build EdelweissFE against your own `Marmot` build, e.g., to develop `Marmot` modules, follow the development
installation in the [EdelweissFE README](https://github.com/Edelweiss-Numerics/EdelweissFE#installation).
EdelweissFE looks for `Marmot` in the active conda environment, or in the directory given by the environment variable
`MARMOT_INSTALL_DIR`. The `Marmot` test cases can then be run with

```bash
run_tests_edelweissfe testfiles/marmot
```

## How to use Marmot with Abaqus

The [Abaqus-MarmotInterface](https://github.com/MAteRialMOdelingToolbox/Abaqus-MarmotInterface) allows to use `Marmot` in Abaqus simulations.

## How to use Marmot with MOOSE

The [chamois App](https://github.com/matthiasneuner/chamois) allows to use `Marmot` directly in [MOOSE](https://github.com/idaholab/moose).
A singularity container recipe is [available](https://github.com/matthiasneuner/chamois-singularity).

## Documentation

The documentation can be found under [https://materialmodelingtoolbox.github.io/Marmot/](https://materialmodelingtoolbox.github.io/Marmot/).

You can also built the documentation locally using Doxygen and Sphinx/Breathe following the instructions below.

```bash
   # install the requirements
   conda install --file doc/requirements.txt
   # build the documentation
   python scripts/buildDocumentation.py
```

## How to cite

If you use Marmot in your research, please cite the version you used. For v26.11:

```bibtex
@software{marmot_v26_11,
  author       = {Dummer, Alexander and Mader, Thomas and Basche, Konstantin and Stathas, Alexandros and Hofer, Paul and Neuner, Matthias},
  title        = {Marmot: MAteRialMOdellingToolbox},
  version      = {v26.11},
  month        = nov,
  year         = 2026,
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.23895565},
  url          = {https://doi.org/10.5281/zenodo.23895565}
}
```

To cite Marmot independently of a version, use the concept DOI
[10.5281/zenodo.17509327](https://doi.org/10.5281/zenodo.17509327), which always resolves to the latest release.
GitHub's *Cite this repository* button provides the same information in other formats.

## License

This library is freely available under the GNU Lesser General Public License, version 2.1 or later
(LGPL-2.1-or-later). Please find the details in the ```LICENSE.md``` file.

## Contributing

Contributions are very welcome! Please see the [contributing guidelines](CONTRIBUTING.md) for more information.
