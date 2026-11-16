.. _bondsliplaws:

Bond-Slip Laws
==============

Bond-slip laws relate the slip of a reinforcement bar relative to its host to the bond stress on the bar surface,
in the local frame of the bar (tangential, normal 1, normal 2). They derive from ``Marmot::Materials::BondSlipLaw``,
live as modules in the materials category and register in the ``Marmot::Registry::BondSlipLawFactory``, from
which the
:ref:`embedded bond element <embeddedbondelement>` creates them by name.

``LINEARELASTICBONDSLIP``
-------------------------

Properties ``[K_t, K_n]``: :math:`\tau_t = K_t s_t`, :math:`\sigma_n = K_n s_n`. Large stiffnesses give a penalty
formulation of perfect bond.

``MODELCODE2010BONDSLIP``
-------------------------

Properties ``[tauMax, s1, s2, s3, alpha, tauF, K0, Kn, (mu, (pMax))]`` (8 to 10).

**From fib Model Code 2010** (Eqs. 6.1-1 to 6.1-4) is only the monotonic envelope of the tangential bond stress
over the slip:

.. math::

   \tau_{\mathrm{env}}(s) = \begin{cases}
     K_0\, s & 0 \le s \le s_0 \\
     \tau_{\max} (s / s_1)^{\alpha} & s_0 < s \le s_1 \\
     \tau_{\max} & s_1 < s \le s_2 \\
     \tau_{\max} - (\tau_{\max} - \tau_f) \frac{s - s_2}{s_3 - s_2} & s_2 < s \le s_3 \\
     \tau_f & s_3 < s
   \end{cases}

**Modelling choices, not from MC2010:**

- The linear branch :math:`K_0 s` up to :math:`s_0` (where :math:`K_0 s_0 = \tau_{\max}(s_0/s_1)^\alpha`)
  regularizes the infinite initial slope of the power law. For :math:`\alpha = 1` the power law is linear and
  :math:`s_0 = 0`.
- Unloading and reloading are elastic with :math:`K_0`, inside the yield condition
  :math:`|\tau_t| \le \tau_{\mathrm{env}}(\kappa) + \mu p` with :math:`\kappa = \max |s|`, the largest *total* slip
  reached in either direction. Beyond it, the slip is plastic (state ``plastic slip``). Because :math:`\kappa` is
  driven by the total slip, the yield stress of a step is known before the check and the return mapping is closed
  form. The damage is symmetric: after slip in one direction, the reversed direction is bounded by the same, already
  reduced bond stress. Cyclic degradation beyond that is not modelled.
- The normal response is elastic with the penalty stiffness :math:`K_n`,
  :math:`\mathbf{t}_n = K_n (\mathbf{s}_n - \mathbf{s}_n^p)`, optionally perfectly plastic with an isotropic cap
  :math:`|\mathbf{t}_n| \le p_{\max}` (radial return). The cap models local crushing of the concrete around the bar;
  :math:`\mathbf{s}_n^p` (state ``plastic normal slip``) is the permanent widening of the hole the bar sits in.
- Coulomb friction: the contact pressure :math:`p = |\mathbf{t}_n|` raises the bound by :math:`\mu p`, with the
  consistent coupling tangent :math:`\partial \tau_t / \partial \mathbf{s}_n`. The magnitude is used because a bar
  pushed sideways presses on the concrete on one side of its hole, whatever the sign of the normal slip. MC2010
  itself accounts for transverse pressure by a factor :math:`\Omega_p` on the envelope, which is not implemented.

With 8 properties, :math:`\mu = 0` and there is no cap (linear normal response, :math:`K_n \ge 0`). A cap
:math:`p_{\max} > 0` requires :math:`K_n > 0`.

Parameters from fib MC2010
~~~~~~~~~~~~~~~~~~~~~~~~~~

MC2010 Table 6.1-1 gives the envelope for ribbed bars, with :math:`f_{ck}` in MPa, slips in mm and
:math:`c_{\mathrm{clear}}` the clear distance between the ribs:

.. list-table::
   :header-rows: 1

   * - failure, bond condition
     - ``tauMax`` [MPa]
     - ``s1``
     - ``s2``
     - ``s3``
     - ``alpha``
     - ``tauF``
   * - pull-out, good bond
     - :math:`2.5\sqrt{f_{ck}}`
     - 1.0
     - 2.0
     - :math:`c_{\mathrm{clear}}`
     - 0.4
     - :math:`0.4\,\tau_{\max}`
   * - pull-out, all other bond
     - :math:`1.25\sqrt{f_{ck}}`
     - 1.8
     - 3.6
     - :math:`c_{\mathrm{clear}}`
     - 0.4
     - :math:`0.4\,\tau_{\max}`
   * - splitting, good bond, unconfined
     - :math:`7.0 (f_{ck}/25)^{0.25}`
     - :math:`s(\tau_{\max})`
     - :math:`s_1`
     - :math:`1.2\, s_1`
     - 0.4
     - 0
   * - splitting, good bond, stirrups
     - :math:`8.0 (f_{ck}/25)^{0.25}`
     - :math:`s(\tau_{\max})`
     - :math:`s_1`
     - :math:`0.5\, c_{\mathrm{clear}}`
     - 0.4
     - :math:`0.4\,\tau_{\max}`

For splitting, :math:`s_1 = s_2 = s(\tau_{\max})` is the slip at the peak, taken from the pull-out branch. The
law requires :math:`s_2 < s_3`, so :math:`s_3` must exceed :math:`s_2`.

:math:`K_0` and :math:`K_n` are not MC2010 quantities but regularization and penalty choices, in N/mm³ for N, mm
units. :math:`K_0` must exceed the secant to the peak, :math:`K_0 > \tau_{\max}/s_1`; a few times that is usual
(:math:`K_0` sets the stiffness at small slip and of unloading). :math:`K_n` should be large compared to the host
stiffness per bar diameter, without ruining the conditioning. Example, :math:`f_{ck} = 30` MPa, good bond, pull-out
failure, units N and mm:

.. code-block:: none

   ** tauMax, s1, s2, s3, alpha, tauF, K0, Kn: fck = 30 MPa, good bond, pull-out (fib MC2010 Tab. 6.1-1)
   13.69, 1.0, 2.0, 10.0, 0.4, 5.48, 200.0, 1000.0

Sign convention
---------------

The slip is the displacement of the bar minus that of the host, in the local frame
:math:`[\mathbf{t}, \mathbf{n}_1, \mathbf{n}_2]`: the tangent :math:`\mathbf{t}` points from the first to the last
bar node, in 2D :math:`\mathbf{n}_1` is :math:`\mathbf{t}` rotated by +90°, in 3D the normals are an arbitrary
right-handed pair (built from the global axis least aligned with :math:`\mathbf{t}`); laws must therefore be
isotropic in the normal plane. A positive bond stress acts on the host in the slip direction and on the bar against
it.

Writing a new bond-slip law
---------------------------

No change to the bond elements or to EdelweissFE is needed: the elements create the law by name from the factory.
Copy ``LinearElasticBondSlip`` as the template; a law is six files plus a section on this page.

1. Create a module ``modules/materials/<Name>`` with a ``module.cmake`` (``REQUIRES MarmotMechanicsCore``), the
   header ``include/Marmot/<Name>.h``, the source ``src/<Name>.cpp``, the registration
   ``src/<Name>Registration.cpp``, and the test ``test/test.cpp`` with its ``test.cmake``.
2. Derive from ``Marmot::Materials::BondSlipLaw``. In the constructor, check the number of properties *before*
   reading them, validate them (throw ``std::invalid_argument``) and declare the state variables in
   ``stateLayout``, then ``stateLayout.finalize()``.
3. Implement ``computeBondStress``: ``state`` holds the values of the beginning of the increment on entry; update
   ``bondStress``, the state variables, ``elasticEnergyDensity`` and ``dissipation``, and return the algorithmic
   tangent ``dBondStress_dSlip``. The slip increment and the time increment are passed for rate-dependent laws; the
   laws here are rate independent and ignore them. ``State`` and ``TimeIncrement`` are small structs of this
   interface (instead of ``state3D``/``timeInfo`` of the continuum materials) because a bond law has its own
   stress measure (a 3-vector traction per bar surface).
4. Register it in a registration source (see ``ModelCode2010BondSlipRegistration.cpp``) with
   ``BondSlipLawFactory::registerBondSlipLaw<Name>("NAME")``.

Verification
------------

Each law has a test in its module, run with ``ctest``: the closed-form envelope on a monotonic path, unloading,
reversal, the Coulomb envelope and the normal cap, the rejection of invalid properties, and the tangent compared to
central differences on every branch, starting from committed states. A new law should at least check its tangent this
way.

In a structure, a linear law can be checked against the analytical pull-out of an elastic bar (axial stiffness
:math:`EA`, bonded length :math:`L`, perimeter :math:`p`) from a rigid host: the end force over the end displacement
is

.. math::

   \frac{F}{u} = EA\, \beta \tanh(\beta L), \qquad \beta = \sqrt{\frac{K_t\, p}{EA}},

to which the embedded bond element converges as the bar elements are refined (EdelweissFE test
``tests/test_embedded_bond.py``).

Implementation
--------------

.. doxygenclass:: Marmot::Materials::BondSlipLaw
   :project: Marmot
   :members:

.. doxygenclass:: Marmot::Materials::LinearElasticBondSlip
   :project: Marmot
   :members:

.. doxygenclass:: Marmot::Materials::ModelCode2010BondSlip
   :project: Marmot
   :members:
