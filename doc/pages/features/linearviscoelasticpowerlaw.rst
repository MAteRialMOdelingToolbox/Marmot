Linear Viscoelastic Power Law model
===================================

Theory
------

The present model is an implementation of an linear viscoelastic power law based on an approximation using a Kelvin chain.
The constitutive law is given in total form as

.. math::

   \sig = \Cel : \epsE  = \Cel : \left( \eps - \epsVE \right)

relating the nominal stress tensor :math:`\sig`
to the elastic strain tensor :math:`\epsE`
in terms of the isotropic fourth order stiffness tensor :math:`\Cel`.
Where
:math:`\eps` is the total strain and
:math:`\epsVE` is the viscoelastic strain.

The rate of the compliance function of the material is given as

.. math::

   \JRate(t) = n m \frac{t^{n-1}}{\tau^n},

where
:math:`n` is the power law exponent,
:math:`m` is the power law compliance parameter,
:math:`t` is the time,
and
:math:`\tau` is a reference time.
The viscoelastic compliance is approximated by a Kelvin Chain.

The evolution of the viscoelastic strain :math:`\epsVE` is defined as

.. math::

   \epsVERate(t) = \int_0^t \JRate(t-t') \DelNu : d\sig(t'),

in terms of the unit compliance tensor :math:`\DelNu`. The unit compliance tensor is defined as

.. math::

   \DelNu = \CelInv(E=1, \nu),

with the Young's modulus :math:`E` and the Poisson's ratio :math:`\nu`.



Implementation
--------------

.. doxygenclass:: Marmot::Materials::LinearViscoelasticPowerLaw
   :allow-dot-graphs:
