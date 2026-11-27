.. _automaticDifferentiation:

Automatic Differentiation
=========================

Marmot computes many derivatives by forward mode automatic differentiation (AD): the consistent tangent of a
material from its stress update, the Jacobian of a residual vector, the gradient and the Hessian of an energy
density. A function is evaluated with a number type that carries derivatives next to its value, and every
arithmetic operation and elementary function applies the chain rule to them. The derivatives are exact up to
round-off; there is no step size to choose as with finite differences.

The number types are Marmot's own and need no further library. Everything on this page lives in the namespace
``Marmot::Math::Differentiation::Automatic``, in three headers of ``MarmotMathCore``:

.. list-table::
   :header-rows: 1
   :align: left

   * - **Header**
     - **Content**
   * - ``MarmotAutomaticDifferentiationTypes.h``
     - number types, seeding, reading derivatives
   * - ``MarmotAutomaticDifferentiation.h``
     - derivatives of functions of scalars and of Eigen vectors
   * - ``MarmotAutomaticDifferentiationForFastor.h``
     - derivatives of functions of Fastor tensors

The code examples below assume

.. code-block:: cpp

   using namespace Marmot::Math::Differentiation::Automatic;
   using namespace Marmot::Math::Tensors::FastorTensors; // Tensor33d, Tensor33t< T >

Number types
------------

A dual number :math:`x + x'\,\varepsilon` extends a real number by a derivative part, with
:math:`\varepsilon^2 = 0`. Evaluating a function with it yields the value and the derivative at once,

.. math::
   f(x + \varepsilon) = f(x) + f'(x)\,\varepsilon .

A hyper-dual number nests this idea: both parts of a dual number are dual numbers again. With two levels there
are two units, :math:`\varepsilon_1^2 = \varepsilon_2^2 = 0` and :math:`\varepsilon_1\varepsilon_2 \neq 0`, and

.. math::
   f(x + \varepsilon_1 + \varepsilon_2) = f(x) + f'(x)\,(\varepsilon_1 + \varepsilon_2)
   + f''(x)\,\varepsilon_1\varepsilon_2 ,

so that the second derivative is exact as well. Each nesting level is one direction of differentiation. A
hyper-dual number therefore gives one directional derivative per evaluation of the function. The vector mode
numbers carry the partial derivatives with respect to all :math:`n` variables instead, and give a whole gradient
or Hessian in one evaluation.

.. list-table::
   :header-rows: 1
   :align: left

   * - **Type**
     - **Carries**
   * - ``HyperDual< order >``
     - value and all mixed derivatives along ``order`` directions
   * - ``Taylor< order >``
     - Taylor coefficients :math:`f^{(k)}/k!`, :math:`k = 0 \dots \mathrm{order}`
   * - ``JacobianDual< n >``
     - value and :math:`n` partial derivatives
   * - ``HessianDual< n >``
     - value, :math:`n` first and :math:`n(n+1)/2` second partial derivatives

- ``HyperDual< order >`` is a hyper-dual number with ``order`` nesting levels, :math:`2^\mathrm{order}` numbers
  in total; ``dual``, ``dual2nd`` and ``dual3rd`` are the orders 1, 2 and 3. Use it for a derivative with respect
  to one variable, for mixed and higher derivatives, and for functions whose arguments carry derivatives already.
- ``Taylor< order >`` is a truncated Taylor polynomial of a function of one variable. Use it for derivatives of
  high order, for which a hyper-dual number becomes expensive.
- ``JacobianDual< n >`` gives gradients and Jacobians, ``HessianDual< n >`` Hessians, in one evaluation.

``VectorXdual``, ``MatrixXdual`` and ``VectorXdual2nd`` are the Eigen vectors and matrices of ``dual`` and
``dual2nd``. All types except ``Taylor`` can be used as the scalar type of Eigen and Fastor objects, also in
expressions that mix them with ``double``.

.. warning::
   An expression of ``JacobianDual`` numbers is not a number yet: it is evaluated when it is assigned to a
   ``JacobianDual``, and until then it refers to the numbers it was built from. Never store such an expression in
   ``auto``, and never return one from a function or lambda with deduced return type — name the type instead.
   The other number types are plain values, for which ``auto`` is fine.

   .. code-block:: cpp

      const JacobianDual< 6 > f = a * b + c; // correct
      const auto              g = a * b + c; // wrong: g refers to temporaries

Differentiating by hand
-----------------------

The derivatives of a number are set and read by free functions; the numbers have no public members to rely on.

Hyper-dual numbers
..................

``seed< level >( x, 1.0 )`` makes ``x`` the variable that nesting level ``level`` differentiates with respect to,
``value( y )`` returns the value of a result and ``derivative< level >( y )`` its mixed derivative along the
levels 1 to ``level``.

.. code-block:: cpp

   dual x( 2.5 );
   seed< 1 >( x, 1.0 ); // x is the independent variable

   const dual y = x * x * exp( x );

   const double f     = value( y );           // x^2 exp( x )
   const double df_dx = derivative< 1 >( y ); // ( 2 x + x^2 ) exp( x )

With two levels, seeding two different variables gives their mixed derivative, and seeding one variable on both
levels gives its second derivative:

.. code-block:: cpp

   dual2nd a( 1.5 ), b( 0.5 );
   seed< 1 >( a, 1.0 ); // level 1 differentiates with respect to a
   seed< 2 >( b, 1.0 ); // level 2 differentiates with respect to b

   const dual2nd y = a * a * sin( b );

   const double dy_da    = derivative< 1 >( y ); // 2 a sin( b )
   const double d2y_dadb = derivative< 2 >( y ); // 2 a cos( b )

   dual2nd x( 1.5 );
   seed< 1 >( x, 1.0 );
   seed< 2 >( x, 1.0 ); // both levels differentiate with respect to x

   const dual2nd z = x * x * x;

   const double dz_dx   = derivative< 1 >( z ); // 3 x^2
   const double d2z_dx2 = derivative< 2 >( z ); // 6 x

A number constructed from a ``double`` is a constant. ``seed< level >( x, 0.0 )`` makes a seeded number a
constant again, which is how the functions below differentiate with respect to one tensor component after the
other.

Vector mode numbers
...................

``seedDirection( x, i )`` makes ``x`` the variable number ``i``; ``partialDerivative( f, i )`` and, for
``HessianDual``, ``secondPartialDerivative( f, i, j )`` read the derivatives of a result.

.. code-block:: cpp

   using scalar = JacobianDual< 2 >;

   scalar x( 1.5 ), y( 0.5 );
   seedDirection( x, 0 ); // x is variable 0
   seedDirection( y, 1 ); // y is variable 1

   const scalar f = x * x * sin( y ) + y; // name the type: not auto

   const double df_dx = partialDerivative( f, 0 ); // 2 x sin( y )
   const double df_dy = partialDerivative( f, 1 ); // x^2 cos( y ) + 1

.. code-block:: cpp

   using scalar = HessianDual< 2 >;

   scalar x( 1.5 ), y( 0.5 );
   seedDirection( x, 0 );
   seedDirection( y, 1 );

   const auto f = x * x * sin( y ); // a plain number type: auto is fine

   const double df_dx    = partialDerivative( f, 0 );          // 2 x sin( y )
   const double d2f_dxdy = secondPartialDerivative( f, 0, 1 ); // 2 x cos( y )

The three-dimensional co-rotational beam element (:ref:`beamelement`) differentiates its local map this way,
with ``JacobianDual< 12 >`` and ``HessianDual< 12 >``.

Taylor numbers
..............

``taylorDerivatives< order >( f, x )`` returns the value and the first ``order`` derivatives of a function of one
variable from one evaluation:

.. code-block:: cpp

   const auto [f, df, d2f, d3f] = taylorDerivatives< 3 >(
     []( const Taylor< 3 >& x ) { return exp( -x ) / ( 1. + x ); },
     0.5 );

A ``Taylor`` number can also be used directly. Its coefficients are read and written with ``operator[]``;
coefficient :math:`k` is the :math:`k`-th derivative divided by :math:`k!`:

.. code-block:: cpp

   Taylor< 2 > t( 0.5 ); // the constant 0.5
   t[1] = 1.0;           // the variable itself: dt/dt = 1

   const Taylor< 2 > g = t * t * t; // g[0] = t^3, g[1] = 3 t^2, g[2] = 6 t / 2!

Marmot uses Taylor numbers for the derivatives of high order of the compliance functions in the Post-Widder
formula of the Kelvin chain models (:doc:`b4`, :doc:`linearviscoelasticpowerlaw`).

Derivatives of functions
------------------------

For the common cases the seeding is already written. Each of the following functions takes the function to
differentiate and the point of evaluation; the number type of the function to differentiate is fixed by the
function chosen.

.. list-table::
   :header-rows: 1
   :align: left

   * - **Function**
     - **Returns**
     - **Number type**
     - **Evaluations**
   * - ``df_dx``
     - :math:`f'(x)`
     - ``dual``
     - 1
   * - ``taylorDerivatives< order >``
     - :math:`f(x),\ f'(x),\ \dots,\ f^{(\mathrm{order})}(x)`
     - ``Taylor< order >``
     - 1
   * - ``dF_dX< n >``
     - :math:`\boldsymbol{F},\ \partial\boldsymbol{F}/\partial\boldsymbol{X}`
     - ``JacobianDual< n >``
     - 1
   * - ``dF_dX``
     - :math:`\boldsymbol{F},\ \partial\boldsymbol{F}/\partial\boldsymbol{X}`
     - ``dual``
     - :math:`n`
   * - ``df_dT``
     - :math:`\partial f/\partial\boldsymbol{T}`
     - ``JacobianDual< nComponents >``
     - 1
   * - ``dF_dT``
     - :math:`\boldsymbol{F},\ \partial\boldsymbol{F}/\partial\boldsymbol{T}`
     - ``JacobianDual< nComponents >``
     - 1
   * - ``SecondOrder::d2f_dT2``
     - :math:`f,\ \partial f/\partial\boldsymbol{T},\ \partial^2 f/\partial\boldsymbol{T}^2`
     - ``HessianDual< dim * dim >``
     - 1
   * - ``SecondOrder::d2f_dTensor_dScalar< dim >``
     - :math:`\partial^2 f/\partial\boldsymbol{T}\partial s`
     - ``dual2nd``
     - :math:`\mathrm{dim}^2`
   * - ``ThirdOrder::d3f_dT3< dim >``
     - :math:`f,\ \partial f/\partial\boldsymbol{T},\ \partial^2 f/\partial\boldsymbol{T}^2,\
       \partial^3 f/\partial\boldsymbol{T}^3`
     - ``dual3rd``
     - 165

Lower case :math:`f` is a scalar function, :math:`x` and :math:`s` are scalars, :math:`\boldsymbol{X}` is an Eigen
vector with :math:`n` entries and :math:`\boldsymbol{T}` a Fastor tensor with ``nComponents`` components, e.g. 9 for
a :math:`3\times3` tensor. ``dF_dX< n >`` takes vectors of fixed size, ``dF_dX`` vectors of dynamic size. The three
last functions are limited to tensors of second order, of size ``dim`` :math:`\times` ``dim``; the 165 evaluations
of ``d3f_dT3`` are those for ``dim`` = 3. A derivative with respect to a tensor carries the indices of the tensor
last: ``dF_dT`` of two second order tensors returns :math:`\partial F_{ij} / \partial T_{kl}` as component
``( i, j, k, l )``.

Three more functions differentiate a function at a point that carries derivatives itself, so that the result
can be differentiated further, e.g. with respect to time or to the variables of an enclosing computation:
``df_dx`` at a ``dual`` with a function of ``dual2nd``, ``dF_dX_2nd`` at a ``VectorXdual`` with a function of
``VectorXdual2nd``, and ``df_dT< order >`` at a tensor of ``HyperDual< order >`` with a function of
``HyperDual< order + 1 >``. They evaluate the function once per variable, ``df_dT< order >``
once more for the value.

Scalars and vectors
...................

.. code-block:: cpp

   const double dfdx = df_dx( []( const dual& x ) { return x * x * exp( x ); }, 2.5 );

The Jacobian of a residual vector, as it is needed for a Newton iteration:

.. code-block:: cpp

   using scalar = JacobianDual< 2 >;
   using Vector = Eigen::Vector< scalar, 2 >;

   const Eigen::Vector2d X( 1.5, 0.5 );

   const auto [R, dR_dX] = dF_dX< 2 >(
     []( const Vector& X_ ) {
       Vector R_;
       R_( 0 ) = X_( 0 ) * X_( 0 ) + X_( 1 ) - 3.;
       R_( 1 ) = X_( 0 ) * X_( 1 ) - 1.;
       return R_;
     },
     X );

If the number of unknowns is known only at run time, write the function with ``VectorXdual`` and call
``dF_dX( F, X )`` with an ``Eigen::VectorXd``. It evaluates the function once per unknown; prefer ``dF_dX< n >``
whenever the size is fixed.

Tensors
.......

A function is best written once as a template of the number type. It can then be evaluated with ``double``
and differentiated to any order:

.. code-block:: cpp

   template < typename T >
   T strainEnergyDensity( const Tensor33t< T >& C, double K, double G )
   {
     const T J  = sqrt( Fastor::determinant( C ) );
     const T I1 = Fastor::trace( C );

     const T psi = K / 8. * ( J - 1. / J ) * ( J - 1. / J ) + G / 2. * ( I1 * pow( J, -2. / 3 ) - 3. );
     return psi;
   }

The gradient, with respect to a ``Tensor33d C``:

.. code-block:: cpp

   using scalar = JacobianDual< 9 >;

   const Tensor33d dPsi_dC = df_dT(
     [&]( const Tensor33t< scalar >& C_ ) { return strainEnergyDensity( C_, K, G ); },
     C );

Value, gradient and Hessian:

.. code-block:: cpp

   using scalar = HessianDual< 9 >;

   const auto [psi, dPsi_dC, d2Psi_dCdC] = SecondOrder::d2f_dT2(
     [&]( const Tensor33t< scalar >& C_ ) { return strainEnergyDensity( C_, K, G ); },
     C );

Up to the third derivative, which
:doc:`compressiblefinitestrainlinearviscoelasticity` needs in every stress update:

.. code-block:: cpp

   const auto [psi, dPsi_dC, d2Psi_dCdC, d3Psi_dCdCdC] = ThirdOrder::d3f_dT3< 3 >(
     [&]( const Tensor33t< dual3rd >& C_ ) { return strainEnergyDensity( C_, K, G ); },
     C,
     true );

``d3f_dT3`` evaluates the function once for every unordered triple of tensor components, since mixed partial
derivatives do not depend on the order of differentiation: 165 times for a :math:`3\times3` tensor. The last
argument, ``isSymmetric``, reduces this to 92. It may only be set if the tensor is symmetric and the function does
not change when its argument is transposed, as is the case for every function of the invariants of the right
Cauchy-Green tensor. ``df_dT< order >`` and ``d2f_dTensor_dScalar`` accept the flag under the same condition; it has
no effect in the functions that evaluate once.

``dF_dT`` differentiates a tensor valued function. In contrast to the other functions, which accept a lambda
directly, it needs the function as a named ``std::function``:

.. code-block:: cpp

   using scalar = JacobianDual< 9 >;

   std::function< Tensor33t< scalar >( const Tensor33t< scalar >& ) > kirchhoffStress =
     [&]( const Tensor33t< scalar >& F_ ) {
       Tensor33t< scalar > tau;
       // ... compute the Kirchhoff stress from the deformation gradient F_
       return tau;
     };

   const auto [tau, dTau_dF] = dF_dT( kirchhoffStress, F ); // dTau_dF( i, j, k, l )

Materials
---------

Two base classes turn a stress update written with AD numbers into a material with a consistent tangent. A
derived material implements ``computeStressAD``; the base class implements ``computeStress``, which evaluates
``computeStressAD`` once and returns the stress together with its derivative.

.. list-table::
   :header-rows: 1
   :align: left

   * - **Base class**
     - **Number type**
     - **Tangent**
     - **Computed with**
   * - ``MaterialHypoElasticAD``
     - ``ADScalar`` = ``JacobianDual< 6 >``
     - :math:`\partial \sig / \partial \eps`
     - ``dF_dX< 6 >``
   * - ``MaterialFiniteStrainAD``
     - ``ADScalar< 3 >`` = ``JacobianDual< 9 >``
     - :math:`\partial \boldsymbol{\tau} / \partial \boldsymbol{F}`
     - ``dF_dT``

Both classes live in ``Marmot::Materials``. The stress and the strain increment of ``MaterialHypoElasticAD`` are
vectors in :ref:`Voigt notation <voigtnotation>`.

Every quantity that depends on the strain increment (or on the deformation gradient) must be an ``ADScalar``.
Energies and state variables are plain numbers and are assigned from values:

.. code-block:: cpp

   void MyMaterial::computeStressAD( state3DAD&       state,
                                     const Vector6AD& dStrain,
                                     const timeInfo&  timeInfo ) const
   {
     // C is the elastic stiffness, a Marmot::Math::Matrix6d
     const Vector6AD dStress = C * dStrain;
     const ADScalar  dEnergy = dStrain.dot( state.stress + 0.5 * dStress );

     state.stress += dStress;

     state.elasticEnergyDensity += value( dEnergy );
     state.dissipation = 0.0;
   }

``value`` returns the value of a number; ``Marmot::Math::Tensors::makeReal`` returns that of a number, and
converts a whole Eigen matrix or Fastor tensor.

The materials :doc:`adlinearelastic`, :doc:`advonmises` and :doc:`adcompressibleneohooke` demonstrate the two
base classes; :doc:`finitestrainisotropicbiotviscoelasticity` and
:doc:`finitestrainorthotropicbiotviscoelasticity` obtain their tangents this way.

Pitfalls
--------

**Branches follow the values.** Comparisons compare the values of the numbers, and the derivative is that of the
branch taken. This is what a stress update needs — the tangent of an elastic step is the elastic one — but it
also means that the derivative at a kink is the derivative of one side.

**Points where the derivative does not exist.** The value is fine there, the derivative is not:

.. list-table::
   :header-rows: 1
   :align: left

   * - **Expression**
     - **Derivative at** :math:`x = 0`
   * - ``sqrt( x )``, ``cbrt( x )``
     - ``inf``
   * - ``sqrt( x * x )``
     - ``nan``
   * - ``abs( x )``
     - 1
   * - ``pow( x, 1.5 )``
     - 0

The second row is the norm of a tensor that vanishes, such as the equivalent stress of a hydrostatic stress
state. Branch on the value before such an expression is evaluated.

**Conversions to double drop the derivatives.** They are explicit for this reason (``value``, ``makeReal``,
``static_cast< double >``). Use them for state variables, energies and convergence checks, never for a quantity
on the way from the variables to the result.

**Newton iterations on AD numbers.** If a Newton iteration updates an AD number with a derivative taken as a
plain number,

.. code-block:: cpp

   x = x - g( x, p ) / dg_dx; // x, p: AD numbers; dg_dx: double

then the derivative :math:`\partial x / \partial p` that ``x`` carries after the update is the exact one,
:math:`-(\partial g / \partial p) / (\partial g / \partial x)`, but evaluated at the iterate before the update. It
lags one iteration behind the value: for :math:`g = x^3 + p\,x - 1` an iteration stopped at a residual of
:math:`10^{-14}` carries a derivative with an error of :math:`10^{-8}`. Either do one more update after
convergence, which makes the derivative as accurate as the value, or iterate with ``double`` and differentiate the
converged solution by the implicit function theorem.

**Seed independent variables only.** ``seed`` overwrites the derivatives a number carries along the level, and
``seedDirection`` expects a number without derivatives.

Elementary functions
--------------------

The hyper-dual numbers, ``JacobianDual`` and ``HessianDual`` support the arithmetic operators among themselves and
with ``double``, the comparison operators, output with ``operator<<`` and

.. code-block:: text

   sqrt  cbrt  exp  log  log10  pow
   sin  cos  tan  asin  acos  atan  atan2  sinh  cosh  tanh
   abs  fabs  min  max  isnan  isinf  isfinite

``pow`` accepts a number as base, as exponent or as both. ``Taylor`` supports the same functions except
``atan2``, ``min`` and ``max``. The functions are found by argument-dependent lookup: call them unqualified, as
``sqrt( x )``, not as ``std::sqrt( x )``, so that the same code works for ``double`` and for every AD number.
For ``min`` and ``max`` this takes a ``using std::max;`` in a function template that is also evaluated with
``double``.

Implementation
--------------

The interface is documented in :doc:`../codedocumentation/mathcore`, the two material base classes in
:doc:`../codedocumentation/mechanicscore` and :doc:`../codedocumentation/finitestrainmechanicscore`. Code outside
of ``MarmotAutomaticDifferentiationTypes.h`` uses only the names on this page, never the implementing classes or
their members.

The use of hyper-dual numbers for the derivatives of finite strain material models is described in
A. Dummer, M. Neuner, P. Gamnitzer, G. Hofstetter (2024). Robust and efficient implementation of finite strain
generalized continuum models for material failure: Analytical, numerical, and automatic differentiation with
hyper-dual numbers. *Computer Methods in Applied Mechanics and Engineering* 426:116987.
