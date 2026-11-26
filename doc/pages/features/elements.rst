Finite elements
===============

This section contains the ready to use available Finite elements.

Every element gives access to the states of its quadrature points by name, through ``Element::getStateView( stateName,
quadraturePoint )``. The index of the quadrature point is checked there, once for all elements: an index outside
``0`` to ``getNumberOfQuadraturePoints() - 1`` throws ``std::invalid_argument``. An element provides the view itself
in the protected ``getStateViewAtQuadraturePoint``, which is only called with a valid index.

.. toctree::
  :maxdepth: 1

  displacementfiniteelement
  displacementfinitestrainelement
  generalgradientenhanceddisplacementfiniteelement
  gradientenhancedfinitestraindisplacementelement
  interfacefiniteelement
  beamelement
  trusselement
  embeddedbondelement
  explicitdynamicsdevices
  energygradients
