.. _heatconduction:

Heat conduction
===============

Transport materials (:cpp:class:`Marmot::Materials::TransportMaterial`, created by name from the
:cpp:class:`Marmot::Registry::TransportMaterialFactory`) give the content, the flux and the source of the fields of
the :ref:`transport elements <transportfiniteelement>`, together with their derivatives with respect to the fields and
their gradients.

``HEATCONDUCTION`` is Fourier heat conduction with an anisotropic conductivity tensor and a conductivity and a heat
capacity depending linearly on the temperature, :math:`\theta = T - T_\mathrm{ref}`:

.. math::

   \mathbf{q} = -\mathbf{k}_0 (1 + \beta_k\theta)\,\nabla T, \qquad
   H(T) = \rho c_0 \big(\theta + \tfrac{1}{2}\beta_c\theta^2\big).

Properties: ``[rho, c0, k11, k22, k33, k12, k13, k23, beta_k, beta_c, T_ref]``, or the short forms
``[rho, c0, k]`` (isotropic, constant) and ``[rho, c0, k, beta_k, beta_c, T_ref]`` (isotropic). A zero heat capacity
gives steady state heat conduction.

.. doxygenclass:: Marmot::Materials::TransportMaterial
   :project: Marmot
   :members:

.. doxygenclass:: Marmot::Materials::HeatConduction
   :project: Marmot
   :members:
