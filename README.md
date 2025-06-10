# Marmot

:warning: [WIP] Note: This project is still in the transfer stage to an open-source project. Thus, many changes will be made in the near future.

## What is Marmot?

```Marmot``` (MAteRialMOdellingToolbox) is a C++-library aiming to provide robust and efficient implementations of state-of-the-art constitutive models for different materials, in particular for quasi-brittle materials such as (sprayed) concrete, rock and soils. It uses modern, object-oriented programming techniques and provides a generic interface that can be easily wrapped into your application. Standard interfaces for commercial finite element codes, such as Abaqus, Plaxis and open source codes like MOOSE or OpenSees, are ready-to-use available.

## Gallery

![Truss in compression using a micropolar von Mises plasticity model](share/truss.gif)

Truss in compression using a micropolar von Mises plasticity model.

![Plane strain compression using a micropolar Drucker-Prager plasticity model](share/plane_strain_gmdruckerprager.gif)

Plane strain compression using a micropolar Drucker-Prager plasticity model

![Triaxial compression using an orthotropic jointed rock plasticity model](share/MultiJoint_Rock.gif)

Triaxial compression using an orthotropic jointed rock plasticity model

## Third-party dependencies

```Marmot``` requires the Eigen (>3.3.8), autodiff (>0.6.0) and Fastor (>6.4.0) libraries:

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

```bash
git clone https://github.com/romeric/Fastor.git
cd Fastor
mkdir build
cd build
cmake ..
sudo make install
```

For ```anaconda``` users, eigen and autodiff can also be installed from the channels.

## How to install Marmot

```Marmot``` can be installed with the following steps:

```bash
git clone https://github.com/MAteRialMOdelingToolbox/Marmot/ 
cd Marmot
mkdir build
cd build
cmake ..
make
sudo make install
```

CMake option ```CMAKE_INSTALL_PREFIX``` allows to specify the installation directory.
For instance:

```bash
git clone https://github.com/MAteRialMOdelingToolbox/Marmot/ 
cd Marmot
mkdir build
cd build
cmake \
    -DCMAKE_INSTALL_PREFIX=/your/special/installationdirectory \
    ..
make
make install
```

## How to use Marmot with Abaqus

The [Abaqus-MarmotInterface](https://github.com/MAteRialMOdelingToolbox/Abaqus-MarmotInterface) allows to use ```Marmot``` in Abaqus simulations.

## How to use Marmot with MOOSE

The [chamois App](https://github.com/matthiasneuner/chamois) allows to use ```Marmot``` directly in [MOOSE](https://github.com/idaholab/moose).
A singularity container recipe is [available](https://github.com/matthiasneuner/chamois-singularity).

## Documentation

The documentation can be found under [https://materialmodelingtoolbox.github.io/Marmot/](https://materialmodelingtoolbox.github.io/Marmot/).

## License

This library is freely available under the LGPLV2 license. Please find the details in the ```LICENSE.md``` file.

## Authors

The principal developers are (in alphabetical order):
* [Alexander Dummer](https://www.uibk.ac.at/bft/mitarbeiter/dummer.html.de) [@alexdummer](https://github.com/alexdummer) (since 2019), University of Innsbruck
* [Paul Hofer](https://www.uibk.ac.at/bft/mitarbeiter/paulhofer.html.de) [@advktEntnschdl](https://github.com/advktEntnschdl) (since 2020), University of Innsbruck
* [Thomas Mader](https://www.uibk.ac.at/bft/mitarbeiter/mader.html) [@maderthomas](https://github.com/maderthomas) (since 2019), University of Innsbruck
* [Matthias Neuner](https://www.uibk.ac.at/bft/mitarbeiter/neuner.html) [@matthiasneuner](https://github.com/matthiasneuner) (since 2015), University of Innsbruck, Stanford University
* [Magdalena Schreter](https://www.uibk.ac.at/bft/mitarbeiter/schreter.html) [@mschreter](https://github.com/mschreter) (since 2015), University of Innsbruck

Contributors are (in alphabetical order):
* [Andreas Brugger](https://www.uibk.ac.at/bft/mitarbeiter/brugger.html.de) [@theBruegge](https://github.com/theBruegge), University of Innsbruck
* [Peter Gamnitzer](https://www.uibk.ac.at/bft/mitarbeiter/gamnitzer.html.de), University of Innsbruck

## Publications (selected)
The results of the following publications were obtained using ```Marmot```:

* On discrepancies between time-dependent nonlinear 3D and 2D finite element simulations of deep tunnel advance: A numerical study on the Brenner Base Tunnel.
M Neuner, M Schreter, P Gamnitzer, G Hofstetter - Computers and Geotechnics, 2020
[https://doi.org/10.1016/j.compgeo.2019.103355](https://www.sciencedirect.com/science/article/abs/pii/S0266352X19304197)

* On the importance of advanced constitutive models in finite element simulations of deep tunnel advance.
M Schreter, M Neuner, D Unteregger, G Hofstetter - Tunnelling and Underground Space Technology, 2018
[https://doi.org/10.1016/j.tust.2018.06.008](https://www.sciencedirect.com/science/article/abs/pii/S0886779818301950)

* A 3D gradient-enhanced micropolar damage-plasticity approach for modeling quasi-brittle failure of cohesive-frictional materials.
M Neuner, P Gamnitzer, G Hofstetter - Computers & Structures, 2020
[https://doi.org/10.1016/j.compstruc.2020.106332](https://www.sciencedirect.com/science/article/pii/S0045794920301358)

* On the prediction of complex shear dominated concrete failure by means of classical and higher order damage-plasticity continuum models.
M Neuner, P Hofer, G Hofstetter - Engineering Structures, 2022
[https://10.1016/j.engstruct.2021.113506](https://www.sciencedirect.com/science/article/pii/S0141029621016072)

* An extended gradient-enhanced damage-plasticity model for concrete considering nonlinear creep and failure due to creep.
A Dummer, M Neuner, G Hofstetter - International Journal of Solids and Structures, 2022
[https://doi.org/10.1016/j.ijsolstr.2022.111541](https://www.sciencedirect.com/science/article/pii/S0020768322000907)

* A gradient enhanced transversely isotropic damage plasticity model for rock - formulation and comparison of different approaches.
T Mader, M Schreter, G Hofstetter - International Journal for Numerical and Analytical Methods in Geomechanics, 2022; 46: 933-960
[https://doi.org/10.1002/nag.3327](https://doi.org/10.1002/nag.3327)
