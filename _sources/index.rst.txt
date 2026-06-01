Marmot Documentation
====================

What is Marmot?
---------------

Marmot (MAteRialMOdellingToolbox) is a C++-library aiming to provide robust and efficient
implementations of state-of-the-art constitutive models for different materials, in particular
for quasi-brittle materials such as (sprayed) concrete, rock and soils. It uses modern,
object-oriented programming techniques and provides a generic interface that can be easily
wrapped into your application. It can be seamlessly used with the
`EdelweissFE <https://github.com/Edelweiss-Numerics/EdelweissFE>`_ finite element code, but
also with commercial and open source finite element codes. Standard interfaces for commercial
finite element codes, such as Abaqus, Plaxis and open source codes like
`MOOSE <https://github.com/idaholab/moose>`_ or OpenSees, are ready-to-use available.

.. image:: ../share/truss.gif
   :alt: Truss in compression using a micropolar von Mises plasticity model

Truss in compression using a micropolar von Mises plasticity model.

.. image:: ../share/plane_strain_gmdruckerprager.gif
   :alt: Plane strain compression using a micropolar Drucker-Prager plasticity model

Plane strain compression using a micropolar Drucker-Prager plasticity model.

.. image:: ../share/MultiJoint_Rock.gif
   :alt: Triaxial compression using an orthotropic jointed rock plasticity model

Triaxial compression using an orthotropic jointed rock plasticity model.

This library is freely available under the LGPLV2 license.
Please find the details in the LICENSE.md file.

.. toctree::
   :maxdepth: 2
   :hidden:

   pages/installation
   pages/interfaces
   pages/features
   pages/documentation
   pages/publications
