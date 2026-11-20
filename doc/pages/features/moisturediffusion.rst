.. _moisturediffusion:

Moisture diffusion
==================

``MOISTUREDIFFUSION`` is a :cpp:class:`Marmot::Materials::TransportMaterial` for the pore relative humidity
:math:`h` of concrete (the field ``relative humidity`` of the :ref:`transport elements <transportfiniteelement>`
``DC*H``):

.. math::

   \frac{\partial w}{\partial h}\dot h - \nabla\cdot\big(D(h)\nabla h\big) = -\dot w_s, \qquad
   D(h) = D_1\Big(\alpha_0 + \frac{1-\alpha_0}{1 + \big(\frac{1-h}{1-h_c}\big)^n}\Big)

(Bažant and Najjar 1972), with the sorption isotherm :math:`w(h)`: linear, :math:`w = K h` (:math:`K = 1`: the
classical form :math:`\dot h = \nabla\cdot(C(h)\nabla h)`; :math:`K = 0`: steady state), or BSB (Brunauer, Skalny and
Bodor 1969; Xi, Bažant and Jennings 1994), :math:`w = C k V_m h / ((1-kh)(1+(C-1)kh))`; and the optional
self-desiccation :math:`w_s(t) = \Delta w_\infty(1 - e^{-t/\tau_s})`, integrated exactly over each increment.

Properties: ``[D1, alpha0, hc, n]`` (linear isotherm, K = 1), ``[D1, alpha0, hc, n, 0, K]``,
``[D1, alpha0, hc, n, 1, C, k, Vm]`` or ``[..., dw_inf, tau_s]`` (10 properties). The tangent includes
:math:`\mathrm{d}D/\mathrm{d}h` and the moisture capacity :math:`\mathrm{d}w/\mathrm{d}h` at the end of the increment.

.. doxygenclass:: Marmot::Materials::MoistureDiffusion
   :project: Marmot
   :members:
