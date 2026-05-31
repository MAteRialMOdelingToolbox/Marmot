# Interfaces

There are several convenient ways to interface with Marmot.
If you want to use Marmot in a finite element framework,
you may want to try one of the existing available interfaces.

## How to use Marmot with EdelweissFE

The [EdelweissFE](https://github.com/Edelweiss-Numerics/EdelweissFE) finite element code is designed to work seamlessly with ```Marmot```.
After the installation of `Marmot`, EdelweissFE can be installed via

```bash
git clone https://github.com/Edelweiss-Numerics/EdelweissFE.git
cd EdelweissFE
conda install --file requirements.txt
pip install .
```

## How to use Marmot with Abaqus

The [Abaqus-MarmotInterface](https://github.com/MAteRialMOdelingToolbox/Abaqus-MarmotInterface) allows to use ```Marmot``` in Abaqus simulations.

## How to use Marmot with MOOSE

The [chamois App](https://github.com/matthiasneuner/chamois) allows to use ```Marmot``` directly in [MOOSE](https://github.com/idaholab/moose).
A singularity container recipe is [available](https://github.com/matthiasneuner/chamois-singularity).

## Create your custom interface

To create your custom interface to Marmot, you can leverage the `MarmotLibrary::MarmotMaterialHypoElasticFactory` and `MarmotLibrary::MarmotElementFactory` factories to create
instances of registered material models and finite elements.
Marmot provides two main base classes for materials:

- `MarmotMaterialHypoElastic` — for small-strain materials formulated in rate form (hypoelastic constitutive relations).
- `MarmotMaterialFiniteStrain` — for large-strain materials formulated in terms of the deformation gradient.

Finite elements are derived from the general parent class `MarmotElement`.

For instance, the following procedure may be employed to create an instance of a `Marmot::Materials::LinearElastic` material model,
which belongs to the class of `MarmotMaterialHypoElastic` materials, and compute the current stress state:
```cpp
#include "Marmot/MarmotMaterialHypoElasticFactory.h"

// Create the instance by material name
auto material = std::unique_ptr<MarmotMaterialHypoElastic>(
        MarmotLibrary::MarmotMaterialHypoElasticFactory::createMaterial(
            "LINEARELASTIC", materialProperties, nMaterialProperties, anArbitraryIdentificationCode ));

// assign a characteristic length parameter (specific to hypoelastic materials)
material->setCharacteristicElementLength( theCharacteristicElementLength );

// set up a material state (stress, strain energy, state variables)
const int             nStateVars = material->getNumberOfRequiredStateVars();
std::vector< double > stateVars( nStateVars, 0.0 );

MarmotMaterialHypoElastic::state3D state = MarmotMaterialHypoElastic::state3D()
state.stateVars           = stateVars.data();

// compute stresses
Marmot::Matrix6d                          dStress_dStrain;
const Marmot::Vector6d                    dStrain  = Marmot::Vector6d::Zero(); // strain increment
const double                              time     = 0.0;                      // current time
const double                              dTime    = 0.01;                     // time increment
const MarmotMaterialHypoElastic::timeInfo timeInfo = { time, dTime };
material->computeStress( state, dStress_dStrain, dStrain, timeInfo );
```

## Registering a material or finite element

Materials are registered using the `MarmotMaterialHypoElasticFactory::registerMaterial` interface.
The function takes a unique name (a string chosen by the user) and the concrete material type as a template parameter.

For instance, the `LinearElastic` material is registered as
```cpp
#include "Marmot/MarmotMaterialHypoElasticFactory.h"

const static bool LinearElasticIsRegistered =
    MarmotLibrary::MarmotMaterialHypoElasticFactory::registerMaterial< Marmot::Materials::LinearElastic >(
        "LINEARELASTIC" );
```

Finite elements are registered using the `MarmotElementFactory::registerElement` interface.
The function takes a unique name (a string chosen by the user) and a factory function.

For instance, the 4-node plane strain displacement finite element "CPE4" is registered as
```cpp
#include "Marmot/MarmotElementFactory.h"

template < class T,
           Marmot::FiniteElement::Quadrature::IntegrationTypes integrationType,
           typename T::SectionType                             sectionType >
MarmotLibrary::MarmotElementFactory::elementFactoryFunction makeFactoryFunction()
{
    return []( int elementID ) -> MarmotElement* { return new T( elementID, integrationType, sectionType ); };
}

const static bool CPE4_isRegistered = MarmotLibrary::MarmotElementFactory::
    registerElement( "CPE4",
                     makeFactoryFunction< DisplacementFiniteElement< 2, 4 >,
                                          FullIntegration,
                                          DisplacementFiniteElement< 2, 4 >::PlaneStrain >() );
```
making use of an auxiliary factory function generator.

Creating a specific element instance and computing the kernel is straightforward:
```cpp
#include "Marmot/MarmotElementFactory.h"
#include "Marmot/MarmotElementProperty.h"

int elementNumber = 1;

// create the instance by element name
auto theElement = std::unique_ptr<MarmotElement>(
        MarmotLibrary::MarmotElementFactory::createElement( "CPE4", elementNumber ) );

// assign nodal coordinates
theElement->assignNodeCoordinates( coordinates );

// assign properties (e.g., thickness)
theElement->assignProperty( ElementProperties( propertiesElement, nPropertiesElement ) );

// assign a material section by material name
// MarmotElement takes care to create the required instances of MarmotMaterial
theElement->assignProperty( MarmotMaterialSection( "LINEARELASTIC", propertiesMaterial, nPropertiesMaterial ) );

// assign state vars to the element
// MarmotElement automatically assigns the state vars also to MarmotMaterial instances
theElement->assignStateVars( stateVars, nStateVars );

// initialization
theElement->initializeYourself();

// compute the kernel
theElement->computeYourself( U, dU, rightHandSide, KMatrix, time, dTime, pNewDT );
```
