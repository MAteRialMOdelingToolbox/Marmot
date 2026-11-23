.. _hughesWingetWrapper:

Hughes-Winget wrapper
=====================

A decorator that makes any small-strain (hypoelastic) material usable wherever a finite-strain
material is expected -- updated-Lagrangian elements such as ``C3D8UL``, and the meshfree
``Displacement/...`` particles and material points.

``Marmot::Materials::MaterialHypoElastic`` and ``Marmot::Materials::MaterialFiniteStrain`` are unrelated interfaces: the former
returns the Cauchy stress and :math:`\partial\boldsymbol{\sigma}/\partial\boldsymbol{\varepsilon}`
for a linearized strain increment, the latter the Kirchhoff stress and
:math:`\partial\boldsymbol{\tau}/\partial\boldsymbol{F}`. ``HughesWingetWrapper`` bridges the two by
integrating the small-strain model along an objective, co-rotational path.

Usage
-----

The wrapper is registered per material under the name of the wrapped model, suffixed with
``/HUGHES-WINGET``. In an EdelweissFE input file::

  *material, name=VONMISES/HUGHES-WINGET, id=myMaterial
  210000, 0.3, 550, 1000, 200, 1400

and in EdelweissMeshfree:

.. code-block:: python

  material = {"material": "VONMISES/HUGHES-WINGET", "properties": np.array([...])}

All material properties are forwarded to the wrapped model unchanged; the wrapper consumes none of
them.

Registering a further material
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Add one line to that material's own ``<Name>Registration.cpp``:

.. code-block:: cpp

  #include "Marmot/MarmotMaterialHughesWinget.h"
  ...
  const static bool VonMisesHughesWingetIsRegistered = Marmot::Registry::MaterialFiniteStrainFactory::registerMaterial<
    HughesWingetWrapper< VonMisesModel > >( "VONMISES/HUGHES-WINGET" );

Theory
------

Let :math:`\boldsymbol{F}^{(n)}` and :math:`\boldsymbol{F}^{(n+1)}` be the deformation gradients at
the beginning and the end of the increment. The velocity gradient increment is evaluated on the
mid-step configuration,

.. math:: \Delta\boldsymbol{l}
   = \left(\boldsymbol{F}^{(n+1)}-\boldsymbol{F}^{(n)}\right)\boldsymbol{F}_\mathrm{mid}^{-1},
   \qquad
   \boldsymbol{F}_\mathrm{mid} = \tfrac{1}{2}\left(\boldsymbol{F}^{(n)}+\boldsymbol{F}^{(n+1)}\right),

and split into its symmetric and skew parts

.. math:: \Delta\boldsymbol{\varepsilon} = \tfrac{1}{2}\left(\Delta\boldsymbol{l}+\Delta\boldsymbol{l}^{T}\right),
   \qquad
   \Delta\boldsymbol{\Omega} = \tfrac{1}{2}\left(\Delta\boldsymbol{l}-\Delta\boldsymbol{l}^{T}\right).

The incremental rotation follows from the Cayley transform

.. math:: \Delta\boldsymbol{R}
   = \left(\boldsymbol{I}-\tfrac{1}{2}\Delta\boldsymbol{\Omega}\right)^{-1}
     \left(\boldsymbol{I}+\tfrac{1}{2}\Delta\boldsymbol{\Omega}\right),

which is orthogonal to machine precision. The stress carried over from the previous increment is
rotated forward,

.. math:: \boldsymbol{\sigma}_\mathrm{rot}
   = \Delta\boldsymbol{R}\,\boldsymbol{\sigma}^{(n)}\,\Delta\boldsymbol{R}^{T},

the wrapped material is evaluated with :math:`\left(\boldsymbol{\sigma}_\mathrm{rot},
\Delta\boldsymbol{\varepsilon}\right)`, and its updated Cauchy stress is pushed forward to the
Kirchhoff stress expected by the consumers,

.. math:: \boldsymbol{\tau} = J\,\boldsymbol{\sigma}^{(n+1)}, \qquad J = \det\boldsymbol{F}^{(n+1)}.

The elastic energy density and the dissipation are converted as well. The finite-strain interface
reports them per unit reference volume, and its hosts store and reload them so; the wrapped material
works per unit current volume. The updated values of the wrapped material are multiplied by :math:`J`.
On the way in, the two differ, because they are different kinds of quantity:

* The **elastic energy** is a function of the state. The wrapped materials update it as
  :math:`\psi^{(n+1)} = \psi_\mathrm{in} + \psi(\boldsymbol{\sigma}^{(n+1)}) - \psi(\boldsymbol{\sigma}^{(n)})`,
  so the stored value :math:`\Psi^{(n)}` is handed in as the density of the state the increment starts
  from, :math:`\Psi^{(n)} / J_n` with :math:`J_n = \det\boldsymbol{F}^{(n)}`. The stored result is
  :math:`\Psi^{(n+1)} = J\,\psi(\boldsymbol{\sigma}^{(n+1)})`, the energy of the current state per unit
  reference volume.
* The **dissipation** is accumulated along the path. Per unit reference volume its rate is :math:`J`
  times the rate per unit current volume, so with the increment :math:`\Delta d` of the wrapped
  material

  .. math:: D^{(n+1)} = D^{(n)} + J\,\Delta d .

  The stored value is therefore handed in as :math:`D^{(n)} / J`, divided by the same :math:`J` it is
  multiplied with on the way out: the round trip leaves the history unchanged and pushes forward only the
  increment. Dividing by :math:`J_n` instead would rescale the whole history by :math:`J / J_n` in every
  increment, and a compaction would decrease the dissipation.

Without any conversion on the way in, a wrapped material that adds its increment to the value it is
given (``VonMises``, ``LinearElastic``, the Rankine models) would see its history multiplied by
:math:`J` in every increment, compounding. In an inverted configuration (:math:`J \le 0`) the energies
are as meaningless as the stress; the values are then passed through unconverted.

Superimposed rigid rotation is reproduced *exactly*, in a single increment and for any rotation
angle below :math:`180^\circ`: for :math:`\boldsymbol{F}^{(n+1)}=\boldsymbol{Q}\boldsymbol{F}^{(n)}`
the increment :math:`\Delta\boldsymbol{l}` is exactly skew, hence
:math:`\Delta\boldsymbol{\varepsilon}=\boldsymbol{0}` and
:math:`\Delta\boldsymbol{R}=\boldsymbol{Q}`.

Algorithmic tangent
^^^^^^^^^^^^^^^^^^^

With :math:`\boldsymbol{M}=\boldsymbol{F}_\mathrm{mid}^{-1}`,
:math:`\boldsymbol{P}=\boldsymbol{I}-\tfrac{1}{2}\Delta\boldsymbol{l}`,
:math:`\boldsymbol{A}=\boldsymbol{I}-\tfrac{1}{2}\Delta\boldsymbol{\Omega}` and
:math:`\boldsymbol{G}=\left(\boldsymbol{I}+\Delta\boldsymbol{R}\right)\boldsymbol{\sigma}^{(n)}\Delta\boldsymbol{R}^{T}`,

.. math:: \frac{\partial\Delta l_{ij}}{\partial F_{kl}} = P_{ik}M_{lj},
   \qquad
   \frac{\partial\sigma_{\mathrm{rot},ij}}{\partial\Delta\Omega_{kl}}
   = \tfrac{1}{2}\left(A^{-1}_{ik}G_{lj}+A^{-1}_{jk}G_{li}\right),

from which, with :math:`\mathbf{C}` the tangent reported by the wrapped material,

.. math:: \frac{\partial\tau_{ij}}{\partial F_{kl}}
   = \sigma^{(n+1)}_{ij}\,J\,F^{-1}_{lk}
   + J\left(\frac{\partial\sigma_{\mathrm{rot},ij}}{\partial\Delta\Omega_{mn}}
            \frac{\partial\Delta\Omega_{mn}}{\partial F_{kl}}
          + C_{ijmn}\frac{\partial\Delta\varepsilon_{mn}}{\partial F_{kl}}\right).

Tangent modes
^^^^^^^^^^^^^

The second template parameter selects how the tangent is obtained.

.. list-table::
   :header-rows: 1
   :align: left

   * - **Mode**
     - **Wrapped-material evaluations**
     - **Consistency**
   * - ``Analytic`` (default)
     - 1
     - Exact while the wrapped model responds elastically. The expression above implicitly assumes
       :math:`\partial\boldsymbol{\sigma}^{(n+1)}/\partial\boldsymbol{\sigma}_\mathrm{rot}=\boldsymbol{I}`,
       so during plastic flow the relative tangent error is bounded by roughly
       :math:`\|\boldsymbol{\sigma}\|/\|\mathbf{C}^\mathrm{e}\|`.
   * - ``Exact``
     - 7
     - Recovers :math:`\boldsymbol{S}=\partial\boldsymbol{\sigma}^{(n+1)}/\partial\boldsymbol{\sigma}_\mathrm{rot}`
       by forward-differencing the six incoming stress components. Fully consistent. Worth the cost
       for softening models under large incremental rotation.
   * - ``Numerical``
     - 11
     - Forward-differences the complete update with respect to :math:`\boldsymbol{F}`. Fully
       consistent and model-agnostic; intended for verification and debugging.

Characteristic element length
-----------------------------

The wrapper forwards ``setCharacteristicElementLength`` to the wrapped model, but derives **no**
length from geometry. Until one is assigned the wrapped model sees ``NaN``, so a model that actually
depends on it (for a mesh-adjusted softening modulus) fails loudly rather than silently regularising
against a wrong length. Models that never read it are unaffected.

In EdelweissMeshfree, assign one per particle -- particles carry no mesh, so nothing else can:

.. code-block:: python

  for particle in model.particles.values():
      particle.setProperty("characteristic element length",
                           particle.getVolumeUndeformed() ** (1.0 / dimension))

.. warning::

   **Only the stress is rotated.** Tensor-valued internal variables of the wrapped model --
   kinematic-hardening back stresses, plastic strain tensors, anisotropic damage tensors -- are
   passed through untouched and are therefore *not* objective under large incremental rotation.
   Models whose internal state is scalar or isotropic-invariant (isotropic hardening, scalar damage)
   are unaffected.

.. doxygenclass:: Marmot::Materials::HughesWingetWrapper
   :allow-dot-graphs:
