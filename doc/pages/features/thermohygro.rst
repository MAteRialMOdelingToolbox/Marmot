.. _thermohygro:

Coupled heat and moisture transport
===================================

``THERMOHYGRO`` is a two-field :cpp:class:`Marmot::Materials::TransportMaterial` for concrete, with the fields
``temperature`` :math:`T` and ``relative humidity`` :math:`h` (in this order; elements ``DC*TH`` of the
:ref:`transport elements <transportfiniteelement>`). It extends :ref:`HEATCONDUCTION <heatconduction>` (properties
0–10) and :ref:`MOISTUREDIFFUSION <moisturediffusion>` (properties 11–20) by couplings of the Bažant–Thonguthai type
(properties 21–31). With :math:`\theta = T - T_\mathrm{ref}`, :math:`T_a = T - T_0`:

.. math::

   w(h,T) &= w_\mathrm{iso}(h - \kappa\theta), \qquad
   \Delta H = H_\mathrm{dry}(T) - H_\mathrm{dry}(T_n) + c_w\,w(h,T)\,(T - T_n), \\
   \mathbf{q}_w &= -D_\mathrm{BN}(h)\,f(T)\,\big(\nabla h + \chi\,h\,\tfrac{\Theta_v}{T_a^2}\nabla T\big), \qquad
   f(T) = \big(\tfrac{T_a}{T_{a,\mathrm{ref}}}\big)^m e^{Q_D/R\,(1/T_{a,\mathrm{ref}} - 1/T_a)}, \\
   \mathbf{q}_T &= -(1 + \beta_{kh}h)\,\mathbf{k}(T)\,\nabla T + \ell\,\mathbf{q}_w, \\
   r_T &= \dot Q(t_e), \quad r_h = -\dot w_s(t_e), \qquad \dot t_e = e^{E_h/R\,(1/T_{a,\mathrm{ref}} - 1/T_a)},

with the hygrothermic coefficient :math:`\kappa = (\partial h/\partial T)_w`, the vapour pressure driving
:math:`\chi` (:math:`\chi = 1`: the moisture flux is driven by :math:`\nabla(h\,p_s(T))`, Clausius–Clapeyron
:math:`\mathrm{d}\ln p_s/\mathrm{d}T = \Theta_v/T_a^2`), the heat capacity of the evaporable water :math:`c_w`, the
latent heat carried by the vapour :math:`\ell = \psi L_v` (Künzel), the heat of hydration
:math:`Q = Q_\infty(1 - e^{-t_e/\tau_q})` and self-desiccation over the equivalent age :math:`t_e` (state variable
``equivalent age``).

Properties (21: no coupling, or 32): ``[rho, c0, k11, k22, k33, k12, k13, k23, beta_k, beta_c, T_ref, D1, alpha0, hc,
n, isotherm, K|C, k, Vm, dw_inf, tau_s, T_0, Q_D/R, m, kappa, chi, c_w, beta_kh, ell, E_h/R, Q_inf, tau_q]``. With all
couplings off, the material reproduces HEATCONDUCTION and MOISTUREDIFFUSION exactly. The tangent has all four blocks of
the capacity, flux and source derivatives. The latent heat leaving through a drying surface is applied with the load
``LatentHeatFilm``.

.. doxygenclass:: Marmot::Materials::ThermoHygro
   :project: Marmot
   :members:
