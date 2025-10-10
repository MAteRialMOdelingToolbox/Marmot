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
in terms of the fourth order stiffness tensor :math:`\Cel`.
Where
:math:`\eps` is the total strain and
:math:`\epsVE` is the viscoelastic strain.

The compliance function of the material is given as

.. math::

   J(t) = \frac{1}{E} + m \left( \frac{t}{\tau}\right)^n,

where
:math:`t` is the time,
:math:`E` is the Young's modulus,
:math:`m` is the power law compliance parameter,
:math:`n` is the power law exponent and
:math:`\tau` is a reference time.

The evolution of the viscoelastic strain :math:`\epsVE` is defined as

.. math::

   \epsVERate(t) = \int_0^t \JRate(t-t') \DelNu : d\sig(t'),

in terms of the unit compliance tensor :math:`\DelNu`.



Implementation
--------------

.. doxygenclass:: Marmot::Materials::LinearViscoelasticPowerLaw
   :allow-dot-graphs:
