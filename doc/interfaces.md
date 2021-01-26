# Interfaces

\page interfaces Interfacing with Marmot

There are several conventient ways to interface with Marmot.
If you want to use Marmot in a finite element framework,
you may want to try one of the existing available interfaces.

## How to use Marmot with Abaqus

The [Abaqus-MarmotInterface](https://github.com/MAteRialMOdelingToolbox/Abaqus-MarmotInterface) allows to use ```Marmot``` in Abaqus simulations.

## How to use Marmot with MOOSE

The [chamois App](https://github.com/matthiasneuner/chamois) allows to use ```Marmot``` directly in [MOOSE](https://github.com/idaholab/moose).
A singularity container recipe is [available](https://github.com/matthiasneuner/chamois-singularity).

## Create your custom interface

To create your custom interface to Marmot, you can leverage the MarmotLibrary::MarmotMaterialFactory and MarmotLibrary::MarmotElementFactory factories to create
instances of registered material models and finite elements.
In general, specialized Marmot materials are derived from the general parent class MarmotMaterial, and
Marmot finite elements are derived from the general parent class MarmotElement. 
In particular for Marmot materials, Marmot strongly makes use of inheritance to specialize 
on specific continuum theories, e.g., MarmotMaterialHypoElastic for hypoelastic materials.

Base classes MarmotMaterial and MarmotElement are considered the most general description of materials and finite elements.
In particular, MarmotMaterial essentially is only defined by a a identification number, a set of material properties and an internal state (represented by a set of state variables) and may describe any material, e.g., a mechanical, a thermal or a couple thermo-mechanical material.
Similarly, MarmotElement makes no assumption on the actual physics described by the finite element, but may represent any finite element.

For materials, physical interpretation is introduced by derivation and specialization of MarmotMaterial.
For instance, the hypoelastic material MarmotMaterialHypoElastic (computes stresses in rate form) is derived from the more general MarmotMaterialMechanical (computes stresses from a deformation state), which itself is a specialization of MarmotMaterial (reads material properties and may have an internal state).

Instances of MarmotMaterial are casted down to their specialization.
For instance, following procedure may be employed to create an instance of a Marmot::Materials::LinearElastic material model,
which belongs to the class of MarmotMaterialHypoElastic materials, and compute the current stress state:
```
// Get the unique material ID of the registered material "LINEARELASTIC"
const auto materialCode = MarmotMaterialFactory::getMaterialCodeFromName( "LINEARELASTIC" );

// Create the instance and cast to a hypoelastic material
auto material = std::unique_ptr<MarmotMaterialHypoElastic> (
        dynamic_cast<MarmotMaterialHypoElastic*> (
            MarmotLibrary::MarmotMaterialFactory::createMaterial( materialCode, materialProperties, nMaterialProperties, anArbitraryIdentificationCode )));

// assign a set of state vars
material->assignStateVars(stateVars, nStateVars);

// assign a characteristic length parameter (specific to hypoelastic materials)
material->setCharacteristicElementLength(theCharacteristicElementLength);

// compute stresses
material->computeStress(stress_in_and_out, dstress_dstrain_out,  dstrain_in, time_in, dtime_in, pNewDT_in_and_out);
```

## Registering a material or finite element

Materials and finite elements are registered using the MarmotMaterialFactory::registerMaterial and MarmotElementFactory::registerElement interfaces.
Both functions take a unique name (a string chosen by the user), a unique integer identifier (chosen by the user), and factory functions.
Those Factory functions are used by Marmot to create specific instances of the material or element, 
providing a custom element number (MarmotElement) or a material number and a set of material properties (MarmotMaterial).

For instance, the 4 node, plane strain displacement finite element "CPE4" is registered as
```
template < class T,
           Marmot::FiniteElement::Quadrature::IntegrationTypes integrationType,
           typename T::SectionType                             sectionType >
MarmotLibrary::MarmotElementFactory::elementFactoryFunction makeFactoryFunction()
{
    return []( int elementID ) -> MarmotElement* { return new T( elementID, integrationType, sectionType ); };
}

MarmotElementFactory:: registerElement( "CPE4",
                         MarmotLibrary::ElementCode::CPE4,
                         makeFactoryFunction< DisplacementFiniteElement< 2, 4 >,
                                              FullIntegration,
                                              DisplacementFiniteElement< 2, 4 >::PlaneStrain >() );
```
making use of an auxiliary factory function generator.

Creating the specific instance and computing the kernel is easy as
```
const auto elementCode = Marmot::getElementCodeFromName( "CPE4" );

int elementNumber = 1;

// create the instance
auto theElement = std::unique_ptr<MarmotElement> ( MarmotLibrary::MarmotElementFactory::createElement(elementCode,  elementNumber) );

// assign properties (e.g, thickness)
theElement->assignProperty( ElementProperties( propertiesElement, nPropertiesElement ) );

// assign a material section
// in the present case: a MarmotMaterial setion
// MarmotElement takes care to create the required instances of MarmotMaterial
theElement->assignProperty( MarmotMaterialSection( materialID, propertiesMaterial, nPropertiesMaterial) );

// assign state vars to the element
// MarmotElement automatically assigns the state vars also to MarmotMaterial instances
theElement->assignStateVars(stateVars, nStateVars);

// initialization
theElement->initializeYourself(coordinates);

// compute the kernel
theElement->computeYourself(U , dU, rightHandSide, KMatrix, time, dTime, pNewDT); 
```
