![build](https://github.com/MAteRialMOdelingToolbox/Marmot/actions/workflows/build_ubuntu.yml/badge.svg)
![clang-format](https://github.com/MAteRialMOdelingToolbox/Marmot/actions/workflows/indent.yml/badge.svg)
[![documentation](https://github.com/MAteRialMOdelingToolbox/Marmot/actions/workflows/sphinx.yml/badge.svg)](https://materialmodelingtoolbox.github.io/Marmot/)
[![license](https://img.shields.io/badge/license-LGPLv2-blue.svg)](LICENSE.md)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

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

`Marmot` requires the
[Eigen](https://gitlab.com/libeigen/eigen) (>3.3.8),
[autodiff](https://github.com/autodiff/autodiff) (>0.6.0)
and [Fastor](https://github.com/romeric/Fastor) (>6.4.0) libraries to be installed.
Detailed instructions on how to install these dependencies can be found [here](https://materialmodelingtoolbox.github.io/Marmot/pages/installation.html).
For anaconda users both `autodiff` and `Eigen` can be installed via conda-forge.

For a quick start we recommend to use anaconda.
A conda environment including all needed dependencies (excluding Fastor) can be created by running

```bash
conda create -n marmot-env -c conda-forge cmake compilers eigen autodiff
```
To activate the environment use

```bash
conda activate marmot-env
```

Fastor needs to be installed manually by executing
```bash
git clone https:\\github.com/romeric/Fastor.git
cd Fastor
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX ..
make install
```

To finally build and install `Marmot` execute
```bash
git clone https:\\github.com/MAteRialMOdelingToolbox/Marmot.git
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

## How to use Marmot with EdelweissFE
The [EdelweissFE](https://github.com/EdelweissFE/EdelweissFE) finite element code is designed to work seamlessly with `Marmot`.
After the installation of `Marmot`, EdelweissFE can be built with `Marmot` support by executing

```bash
git clone https:\\github.com/Edelweiss-FE/EdelweissFE.git
cd EdelweissFE
conda install --file requirements.txt
pip install .
```
The build can be tested by running
```bash
run_tests_edelweissfe testfiles
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

## License

This library is freely available under the LGPLV2 license. Please find the details in the ```LICENSE.md``` file.

## Contributing

Contributions are very welcome! Please see the [contributing guidelines](CONTRIBUTING.md) for more information.
