# Marmot

:warning: [WIP] Note: This project is still in the transfer stage to an open-source project. Thus, many changes will be made in the near future.

## What is Marmot?

```Marmot``` (MAteRialMOdellingToolbox) is a C++-library aiming to provide robust and efficient implementations of state-of-the-art constitutive models for different materials, in particular for quasi-brittle materials such as (sprayed) concrete, rock and soils. It uses modern, object-oriented programming techniques and provides a generic interface that can be easily wrapped into your application. Standard interfaces for commercial finite element codes, such as Abaqus, Plaxis and OpenSees, are ready-to-use available.

## Gallery

![Truss in compression using a micropolar von Mises plasticity model](share/truss.gif)

Truss in compression using a micropolar von Mises plasticity model.

![Plane strain compression using a micropolar Drucker-Prager plasticity model](share/plane_strain_gmdruckerprager.gif)

Plane strain compression using a micropolar Drucker-Prager plasticity model

## Third-party dependencies

If you would like to have the full experience with ```Marmot``` the Eigen (>3.3.8) library has to be installed. 

## How to install marmot

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
## How to use marmot with Abaqus

todo

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
The results of the following three publications were obtained using ```marmot```. 
*
* 
*
<!--
## Internal

This repository depends on several sub-repositories which can be cloned by

`git clone https://github.com/MAteRialMOdelingToolbox/marmot --recurse-submodules`

A new submodule, e.g. a material, can be added by 

`cd modules/materials &&
git submodule add https://github.com/MAteRialMOdelingToolbox/LinearElastic`

Submodules can be updated in your local repo by 

`git submodule foreach git pull origin master`

### Add a UIBK gitlab project to github

```cd existingRepoUIBK
git remote add github https://github.com/MAteRialMOdelingToolbox/NewRepoGithub 
push -u github --all
git push -u github --tags
```

### Update submodules

Submodules are automatically updated every day at noon (12:00). Otherwise, you may re-run the `update_submodules.yml` action.


### Documentation

The documentation must be created currently by hand

`cd doc/doxygen && doxygen dconfig`
-->
