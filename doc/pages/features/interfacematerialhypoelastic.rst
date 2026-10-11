.. _interfacematerialhypoelastic:

Interface material (hypoelastic bulk material)
==============================================

Theory
------

``MarmotInterfaceMaterialHypoElastic`` turns **any registered hypoelastic bulk material** into the
constitutive law of a thin interface layer of thickness :math:`h`. It is the material counterpart of the
:ref:`Interface Finite Element <interfacefiniteelement>`: the element provides the displacement jump and the
surface gradient, the material returns the conjugate interface quantities and their tangents.

Kinematics of the layer
^^^^^^^^^^^^^^^^^^^^^^^

For an interface with unit normal :math:`\mathbf{n}`, jump :math:`[\![\mathbf{u}]\!] = \mathbf{u}^{+} - \mathbf{u}^{-}`
and average surface gradient :math:`\bar{\nabla}_s \mathbf{u} = \tfrac{1}{2}(\nabla_s \mathbf{u}^{+} + \nabla_s \mathbf{u}^{-})`,
the strain of the layer is approximated as

.. math::

   \Delta\eps = \operatorname{sym}\Big( \tfrac{1}{h}\, [\![\Delta\mathbf{u}]\!] \otimes \mathbf{n}
                + \Delta\bar{\nabla}_s \mathbf{u} \Big).

This strain increment is passed to the bulk material together with the time increment.
The bulk material returns the stress :math:`\sig` and its tangent :math:`\mathbb{C} = \partial\Delta\sig/\partial\Delta\eps`.

Conjugate quantities
^^^^^^^^^^^^^^^^^^^^

The interface returns the traction-like **force** and the **surface stress**,

.. math::

   \mathbf{t} = \sig\, \mathbf{n}, \qquad \mathbf{s} = h\, \sig,

where the surface stress is stored scaled by the thickness :math:`h` and unscaled again when it is
read back as the bulk stress at the start of the next increment.

Tangents
^^^^^^^^

The bulk tangent :math:`\mathbb{C}` is condensed with respect to the normal direction
(``InterfaceMaterialHelperFunctions::calculateInterfaceMaterialParameters``). This yields the four operators
:math:`\hat{\mathbf{Q}}`, :math:`\hat{\mathbb{Z}}`, :math:`\hat{\mathbf{H}}` and :math:`\hat{\mathbb{Y}}`,
which couple the jump and the surface strain to the force and the surface stress. They are scaled with the
thickness to the quantities handed to the element:

.. math::

   \mathbf{Q} = \tfrac{1}{h}\, \hat{\mathbf{Q}}, \qquad
   \mathbb{Z} = h\, \hat{\mathbb{Z}}, \qquad
   \mathbf{H} = \hat{\mathbf{H}}, \qquad
   \mathbb{Y} = h\, \hat{\mathbb{Y}}.

Their roles in the element tangent are described in the :ref:`Interface Finite Element <interfacefiniteelement>`
page.

Properties
----------

The interface material is created with the **name of the bulk material** and the property array

.. code-block:: text

   [ E, nu, h, <remaining properties of the bulk material> ]

that is, the thickness :math:`h` is inserted as third entry after :math:`E` and :math:`\nu`, and removed again
before the array is passed to the bulk material. For example, the
:ref:`Linear Viscoelastic Wiechert model <linearviscoelasticwiechert>` is used as
``[E, nu, h, m, n, nMaxwell, minTau, timeToDays]``. The state variables of the bulk material are
stored in the state-variable block ``baseMaterialStateVars``.

Implementation
--------------

.. doxygenclass:: MarmotInterfaceMaterialHypoElastic
   :allow-dot-graphs:

.. doxygennamespace:: Marmot::Materials::InterfaceMaterialHelperFunctions
