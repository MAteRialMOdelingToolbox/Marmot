Gradient-enhanced compressible Neo-Hooke damage
===============================================

A minimal, fully consistent reference model for ``MarmotMaterialGradientEnhancedFiniteStrain``: the compressible
Neo-Hookean solid of :doc:`compressibleneohooke` with isotropic damage driven by the nonlocal field, in the spirit of
the implicit-gradient damage model of Peerlings et al. (1996). It is registered as
``GRADIENTENHANCEDCOMPRESSIBLENEOHOOKEDAMAGE``.

Theory
------

.. math::

   \boldsymbol{\tau} = (1 - D(\kappa))\,\boldsymbol{\tau}_0(\boldsymbol{F}),\qquad
   \kappa = \max_t\left(\kappa_0,\ \bar N\right),\qquad
   L = \sqrt{2\,\psi_0(\boldsymbol F)/E},\qquad E = \frac{9KG}{3K+G},

where :math:`\psi_0` and :math:`\boldsymbol\tau_0` are the energy density and Kirchhoff stress of the Pence-Gou
potential (variant B). :math:`L` is an energy-equivalent strain, which reduces to the axial strain in small-strain
uniaxial stress, and :math:`\bar N` solves :math:`\bar N - l^2\nabla^2\bar N = L`. The damage law is

.. math::

   D(\kappa) = 1 - \frac{\kappa_0}{\kappa}\exp\left(-\frac{\kappa-\kappa_0}{\kappa_f-\kappa_0}\right)
   \quad\text{for}\quad \kappa > \kappa_0,\qquad D = 0\ \text{otherwise}.

All four tangents are analytic. The dissipation is cumulative: the incoming value is incremented by the energy
released by the damage increment, :math:`\psi_0\,\Delta D`.

Material parameters
-------------------

.. list-table::
   :header-rows: 1

   * - Index
     - Symbol
     - Description
   * - 0
     - :math:`K`
     - bulk modulus
   * - 1
     - :math:`G`
     - shear modulus
   * - 2
     - :math:`\kappa_0`
     - damage threshold (equivalent strain), :math:`\kappa_0 > 0`
   * - 3
     - :math:`\kappa_f`
     - softening parameter, :math:`\kappa_f > \kappa_0`
   * - 4
     - :math:`l`
     - nonlocal radius
   * - 5
     - :math:`\rho`
     - density (optional)

State variable: ``kappa``, the history maximum of the nonlocal field.

Implementation
--------------

.. doxygenclass:: Marmot::Materials::GradientEnhancedCompressibleNeoHookeDamage
   :allow-dot-graphs:
