.. _gradientEnhancedHughesWingetWrapper:

Gradient-enhanced Hughes-Winget wrapper
========================================

The gradient-enhanced sibling of the :ref:`Hughes-Winget wrapper <hughesWingetWrapper>`. It makes a
small-strain **gradient-enhanced** (implicit-gradient / nonlocal damage) material --
``MarmotMaterialGeneralGradientEnhancedHypoElastic<1>`` -- usable wherever a
``MarmotMaterialGradientEnhancedFiniteStrain`` is expected: the gradient-enhanced updated-Lagrangian
element and the gradient-enhanced meshfree particle and material point.

It works the same way as its local sibling -- the wrapped material is driven by the objective strain
increment of the Hughes-Winget algorithm, evaluated on the mid-step configuration, with the carried-over
stress rotated forward by the Cayley transform of the incremental spin and pushed to the Kirchhoff
stress -- so the kinematics are not repeated here; see the local wrapper's Theory section. What this
page documents is what is specific to the gradient-enhanced case: how the nonlocal field crosses the
two interfaces, and the three limitations that follow from wrapping a small-strain rate law rather than
a total-strain law.

.. note::

   ``MarmotMaterialGradientEnhancedFiniteStrain``, the interface this wrapper implements, is itself a
   **reconstruction**: the original header was written by the author of the gradient-enhanced
   finite-strain material point / particle / displacement element, but was never published to any
   repository. Everything those consumers observe (member names, struct shapes, call signatures) is
   fixed by them; everything else -- in particular the *order* of the four scalars carried by
   ``ConstitutiveResponse`` -- is inferred from the sibling interfaces, since the consumers only ever
   pass zeros there and read the members back by name. Should the original resurface, it takes
   precedence over this reconstruction.

Usage
-----

The wrapper is registered per material under the name of the wrapped model, suffixed with
``/HUGHES-WINGET``, with ``/HUGHES-WINGET/EXACT-TANGENT`` and ``/HUGHES-WINGET/NUMERICAL-TANGENT``
variants for the two more expensive tangent modes (see :ref:`hughesWingetWrapper`'s Tangent modes
section -- the same three modes apply here). Currently only **GCDP** is registered this way. In an
EdelweissFE input file::

  *material, name=GCDP/HUGHES-WINGET, id=myMaterial
  <GCDP material properties, unchanged>

and in EdelweissMeshfree:

.. code-block:: python

  material = {"material": "GCDP/HUGHES-WINGET", "properties": np.array([...])}

All material properties are forwarded to the wrapped model unchanged; the wrapper consumes none of
them.

Registering a further material
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Add one line to that material's own ``<Name>Registration.cpp``:

.. code-block:: cpp

  #include "Marmot/MarmotMaterialGradientEnhancedHughesWinget.h"
  ...
  const static bool GCDPHughesWingetIsRegistered = MarmotMaterialGradientEnhancedFiniteStrainFactory::
    registerMaterial< GradientEnhancedHughesWingetWrapper< GCDPModel > >( "GCDP/HUGHES-WINGET" );

Before doing so for a new material, read the three limitations below -- they are the reason not every
``MarmotMaterialGeneralGradientEnhancedHypoElastic<1>`` model should be registered this way.

The nonlocal field across the two interfaces
---------------------------------------------

On top of the Hughes-Winget stress update, the wrapper carries the nonlocal field across the two
interfaces:

.. list-table::
   :header-rows: 1
   :align: left

   * - this interface reports
     - from the wrapped material
   * - ``response.L``
     - ``KLocal(0)``, a total, which is what both sides mean
   * - ``response.nonLocalRadius``
     - :math:`\sqrt{c(0)}`
   * - ``tangents.dTau_dN``
     - :math:`J\,\partial\boldsymbol{\sigma}/\partial\bar{N}`
   * - ``tangents.dL_dF``
     - :math:`\partial K^{\mathrm{local}}/\partial\Delta\boldsymbol{\varepsilon}`, chained through the
       Hughes-Winget kinematics
   * - ``tangents.dL_dN``
     - ``dKLocalddK(0,0)``

The wrapped material wants the nonlocal field **and its increment**, while the finite-strain interface
hands over only the total. The increment is therefore formed inside the wrapper, against a value of the
last accepted increment carried in the wrapper's own state (``HughesWinget_N_n``) -- one slot more than
the local sibling needs.

.. note::

   **The nonlocal balance is formulated in the material (reference) configuration** by the consumers of
   this interface: they take the gradients with respect to :math:`\boldsymbol{X}` and integrate over the
   reference volume. So :math:`c` is handed over as the material constant the wrapped model reports,
   neither pushed forward nor weighted by :math:`J`. The stress, by contrast, is the Kirchhoff stress
   :math:`\boldsymbol{\tau}=J\boldsymbol{\sigma}` -- which is what integrating against the reference
   volume requires. Do not push :math:`c` forward, and do not leave the stress in its Cauchy form.

Limitations
-----------

.. warning::

   **The wrapped material must be incremental in stress.** It has to update the Cauchy stress it is
   handed, not recompute one from a stored strain: the wrapper's whole mechanism is to hand over the
   forward-rotated stress of the last increment, and a material that ignores that argument discards the
   rotation along with it. ``GCDPModel`` qualifies -- it maps ``res.stress`` and updates it in place,
   and its entire state is four scalars.

   ``AT2PhaseField`` does **not**: it keeps a six-component ``strain`` state and returns
   :math:`g(\varphi)\,\mathbb{C}:\boldsymbol{\varepsilon}` from it, ignoring the incoming stress
   entirely. Wrapped, its Kirchhoff stress does not rotate at all under a rigid rotation -- measured.
   It is therefore deliberately **not** registered with this wrapper, and a total-strain-based model
   must not be either.

.. warning::

   **Only the stress is rotated.** Tensor-valued internal variables of the wrapped material pass
   through untouched and are therefore *not* objective under large incremental rotations. GCDP is
   unaffected: its whole state is four scalars (``alphaP``, ``alphaD``, ``omega``, ``I1p``).

.. warning::

   :math:`\partial c/\partial\bar{N}` **cannot be forwarded** -- the finite-strain interface has no slot
   for it. Harmless for a material whose :math:`c` is constant (GCDP and AT2PhaseField both are); a
   material built on ``MarmotDecreasingInteractions`` would keep a correct residual but lose that term
   of its consistent tangent.

What is verified
-----------------

``TestMarmotMaterialGradientEnhancedHughesWinget.cpp`` drives the wrapper with a minimal, genuinely
incremental test-only material and checks:

- Superimposed rigid rotation, on top of a preloaded anisotropic stretch, co-rotates the Kirchhoff
  stress exactly, to :math:`10^{-10}` relative to the stress norm, and leaves the local driving force
  (to the same tolerance) and the nonlocal radius (to :math:`10^{-12}`) unchanged.
- All four tangent blocks (:math:`\partial\boldsymbol{\tau}/\partial\boldsymbol{F}`,
  :math:`\partial\boldsymbol{\tau}/\partial\bar{N}`, :math:`\partial L/\partial\boldsymbol{F}`,
  :math:`\partial L/\partial\bar{N}`) against a forward-differenced (``Numerical``) oracle, to a
  relative tolerance of :math:`10^{-5}` -- once in a purely elastic regime with the ``Analytic`` mode,
  and once along a history that reaches a partially degraded, sheared state with the ``Exact`` mode.
- Small-strain agreement: driven along the same 20-increment history with a ramping nonlocal field, the
  wrapper's Cauchy stress, local driving force and nonlocal radius agree with the bare small-strain
  material to :math:`10^{-5}` (stress; the same order of deviation between the mid-step objective rate
  and the naive additive strain documented for the local wrapper), :math:`10^{-8}` (local driving
  force) and :math:`10^{-12}` (nonlocal radius).
- The ``dK = N - N_n`` bookkeeping: a recording test material confirms the wrapper hands it the
  increment since the last accepted call, not the total field, across three successive calls.

.. doxygenclass:: Marmot::Materials::GradientEnhancedHughesWingetWrapper
   :allow-dot-graphs:
