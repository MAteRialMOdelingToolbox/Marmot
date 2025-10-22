.. _linearviscoelasticpowerlaw:

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

Here we approximate a target **power-law creep compliance** :math:`J(t)` by a **Kelvin chain** (generalized Kelvin–Voigt / Prony series) so the viscoelastic strain rate becomes a finite sum of exponentials that’s efficient to evaluate. The continuous form

.. math::

   J(t)=J_0+\int_0^\infty L(\tau)\,\bigl(1-e^{-t/\tau}\bigr)\,d\tau,

is replaced by

.. math::

   J_N(t)=J_0+\sum_{i=1}^N J_i\bigl(1-e^{-t/\tau_i}\bigr), \qquad
   \dot J_N(t)=\sum_{i=1}^N \frac{J_i}{\tau_i}\,e^{-t/\tau_i},

where :math:`J_i` and :math:`\tau_i` are the partial compliances and retardation times. For a specific compliance target  the :math:`\{J_i,\tau_i\}` are chosen so that :math:`J_N(t)` tracks :math:`J(t)` over the time window of interest, while the rate kernel :math:`\dot J_N` approximates

.. math::

   \dot J(t)=\int_0^\infty \frac{L(\tau)}{\tau}\,e^{-t/\tau}\,d\tau.

The **Post–Widder inversion formula** provides a way to recover an approximate retardation spectrum from a known compliance function :math:`J(t)`.

.. math::

   L(\tau) \approx \frac{(-1)^{k}}{(k)!}
   \left(\frac{k}{\tau}\right)^{k}
   J^{(k)}\!\left(\frac{\tau}{k+1}\right),
   \qquad k \to \infty.

The order :math:`k` defines the order of approximation of the actual continous spectrum.
If it tends to infinity the continous spectrum is recovered.
The discrete spectrum of the Kelvin chain is here used to calculate the stiffnesses or compliances for
chosen retardation times by fitting the approximation of the continous spectrum.

A detailed description of the employed approach can be found, e.g., in Bazant & Jirasek (2018) *Creep and Hygrothermal Effects in Concrete Structures*.


Implementation
--------------

.. doxygenclass:: Marmot::Materials::LinearViscoelasticPowerLaw
   :allow-dot-graphs:
