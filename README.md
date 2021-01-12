# Marmot

## What is Marmot?

:warning: [WIP] Note: This project is still in the transfer stage to an open-source project. Thus, many changes will be made in the near future.

```Marmot``` (MAteRialMOdellingToolbox) is a C++-library aiming to provide robust and efficient implementations of state-of-the-art constitutive models for different materials, in particular for quasi-brittle materials such as (sprayed) concrete, rock and soils. It uses modern, object-oriented programming techniques and provides a generic interface that can be easily wrapped into your application. Standard interfaces for commercial finite element codes, such as Abaqus, Plaxis and OpenSees, are ready-to-use available.

## Gallery

![Truss in compression using a micropolar von Mises plasticity model](share/truss.gif)

Truss in compression using a micropolar von Mises plasticity model.

![Plane strain compression using a micropolar Drucker-Prager plasticity model](share/plane_strain_gmdruckerprager.gif)

Plane strain compression using a micropolar Drucker-Prager plasticity model

## Third-party dependencies

If you would like to have the full experience with ```Marmot``` the Eigen (>3.3.8) library has to be installed. 

## How to install Marmot

```Marmot``` including all submodules can be installed with the following steps:

```bash
git clone --recurse-submodules https://github.com/MAteRialMOdelingToolbox/Marmot/ 
cd Marmot
mkdir build
cd build
cmake ..
make
sudo make install
```

CMake options ```CORE_MODULES```, ```ÈLEMENT_MODULES``` and ```MATERIAL_MODULES``` 
allow to specify the modules which should be compiled, either by passing a 
```semicolon seperated list```, option ```none``` or option ```all``` (default).
For instance:

```bash
git clone --recurse-submodules https://github.com/MAteRialMOdelingToolbox/Marmot/ 
cd Marmot
mkdir build
cd build
cmake \
    -DCORE_MODULES='MarmotMechanicsCore;MarmotMicromorphicCore' \
    -DELEMENT_MODULES='none' \
    -DMATERIAL_MODULES='all' \
    ..
make
sudo make install
```

## How to use Marmot with Abaqus

The [Abaqus-MarmotInterface](https://github.com/MAteRialMOdelingToolbox/Abaqus-MarmotInterface) allows to use ```Marmot``` in Abaqus simulations.

## Documentation

The documentation can be found under https://materialmodelingtoolbox.github.io/marmot/.

## License

This library is freely available under the LGPLV2 license. Please find the details in the ```LICENSE.md``` file.

## Authors

The principal developers are (in alphabetical order):
* [Alexander Dummer](https://www.uibk.ac.at/bft/mitarbeiter/dummer.html.de) [@alexdummer](https://github.com/alexdummer) (since 2019), University of Innsbruck
* [Thomas Mader](https://www.uibk.ac.at/bft/mitarbeiter/mader.html) [@maderthomas](https://github.com/maderthomas) (since 2019), University of Innsbruck
* [Matthias Neuner](https://www.uibk.ac.at/bft/mitarbeiter/neuner.html) [@matthiasneuner](https://github.com/matthiasneuner) (since 2015), University of Innsbruck
* [Magdalena Schreter](https://www.uibk.ac.at/bft/mitarbeiter/schreter.html) [@mschreter](https://github.com/mschreter) (since 2015), University of Innsbruck

Contributors are (in alphabetical order):
* [Andreas Brugger](https://www.uibk.ac.at/bft/mitarbeiter/brugger.html.de) [@theBruegge](https://github.com/theBruegge), University of Innsbruck
* [Peter Gamnitzer](https://www.uibk.ac.at/bft/mitarbeiter/gamnitzer.html.de), University of Innsbruck

## Publications (selected)
The results of the following publications were obtained using ```Marmot```:

On discrepancies between time-dependent nonlinear 3D and 2D finite element simulations of deep tunnel advance: A numerical study on the Brenner Base Tunnel.
M Neuner, M Schreter, P Gamnitzer, G Hofstetter - Computers and Geotechnics, 2020
[https://doi.org/10.1016/j.compgeo.2019.103355](https://www.sciencedirect.com/science/article/abs/pii/S0266352X19304197)

On the importance of advanced constitutive models in finite element simulations of deep tunnel advance.
M Schreter, M Neuner, D Unteregger, G Hofstetter - Tunnelling and Underground Space Technology, 2018
[https://doi.org/10.1016/j.tust.2018.06.008](https://www.sciencedirect.com/science/article/abs/pii/S0886779818301950)

A 3D gradient-enhanced micropolar damage-plasticity approach for modeling quasi-brittle failure of cohesive-frictional materials.
M Neuner, P Gamnitzer, G Hofstetter - Computers & Structures, 2020
[https://doi.org/10.1016/j.compstruc.2020.106332](https://www.sciencedirect.com/science/article/pii/S0045794920301358)
