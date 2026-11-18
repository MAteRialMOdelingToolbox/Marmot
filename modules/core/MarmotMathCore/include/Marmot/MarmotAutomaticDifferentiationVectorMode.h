/* ---------------------------------------------------------------------
 *                                       _
 *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_
 * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
 * | | | | | | (_| | |  | | | | | | (_) | |_
 * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
 *
 * Unit of Strength of Materials and Structural Analysis
 * University of Innsbruck,
 * 2020 - today
 *
 * festigkeitslehre@uibk.ac.at
 *
 * This file is part of the MAteRialMOdellingToolbox (marmot).
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of marmot.
 * ---------------------------------------------------------------------
 */

/** @file MarmotAutomaticDifferentiationVectorMode.h
 *  @brief Marmot's own vector-mode forward automatic differentiation with expression templates.
 *
 *  VectorMode::Dual< nDirections > carries a value and its partial derivatives with respect to nDirections independent
 *  variables, so that one evaluation of a function yields a full Jacobian row. It is tailored to the first-derivative
 *  drivers df_dT, dF_dT and dF_dX< n > (and hence the tangents of the AD materials), where nDirections is small and
 *  known at compile time (6 or 9):
 *
 *  - Arithmetic does not compute derivative arrays for intermediate results. It returns lightweight expression nodes,
 *    and the derivatives are computed in a single loop over the directions when an expression is assigned to a Dual.
 *  - Values are evaluated eagerly: when a node is created, it computes its value and the local partial derivatives
 *    with respect to its operands once. The loop over the directions then only combines the operands' derivatives
 *    linearly, i.e., elementary functions such as sqrt, exp, log or pow are never re-evaluated per direction.
 *  - Assignments may refer to the assigned number, e.g. `x = x * y + z`.
 *  - The partial derivatives are stored before the value, so that pairs of them are 16-byte aligned for SIMD.
 *
 *  Expression nodes keep their sub-expressions by value, but refer to the Dual numbers they are built from. An
 *  expression must therefore not outlive its full-expression: never store one in `auto`, and never return one from a
 *  function or lambda with a deduced return type. Name the number type (Dual, or JacobianDual) instead.
 *
 *  Use this header through the backend-neutral vocabulary of MarmotAutomaticDifferentiationTypes.h (JacobianDual,
 *  value, seedDirection, partialDerivative).
 */

#pragma once
#include <Eigen/Core>
#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <ostream>
#include <type_traits>
#include <utility>

namespace Marmot::Math::Differentiation::Automatic::VectorMode {

  /** @brief Common base of Dual and all expression nodes, which identifies them as vector-mode expressions */
  struct ExpressionBase {};

  /** @brief A vector-mode expression: a Dual or an expression node
   *  @tparam E type to check; references and cv-qualifiers are ignored
   *
   *  Every expression provides `nDirections`, `val()` and `dx( i )`.
   */
  template < typename E >
  concept Expression = std::is_base_of_v< ExpressionBase, std::remove_cvref_t< E > >;

  /** @brief An arithmetic constant (e.g. double or int), which carries no derivatives
   *  @tparam S type to check; references and cv-qualifiers are ignored
   */
  template < typename S >
  concept Constant = std::is_arithmetic_v< std::remove_cvref_t< S > >;

  template < std::size_t nDirections_ >
  class Dual;

  /** @brief True if T is a vector-mode Dual
   *  @tparam T type to check
   */
  template < typename T >
  inline constexpr bool isDual = false;

  /** @brief True for vector-mode Dual numbers
   *  @tparam nDirections_ number of directions
   */
  template < std::size_t nDirections_ >
  inline constexpr bool isDual< Dual< nDirections_ > > = true;

  /** @brief How an expression node stores an operand of type E (as deduced by a forwarding reference)
   *  @tparam E operand type
   *
   *  A Dual is referred to, since copying its derivatives would defeat the purpose of the expression templates; nodes
   *  are copied, since they are small and usually temporaries.
   */
  template < typename E >
  using Operand = std::
    conditional_t< isDual< std::remove_cvref_t< E > >, const std::remove_cvref_t< E >&, std::remove_cvref_t< E > >;

  /** @brief Number of directions of an expression
   *  @tparam E expression type; references and cv-qualifiers are ignored
   */
  template < Expression E >
  inline constexpr std::size_t directionsOf = std::remove_cvref_t< E >::nDirections;

  /** @brief Expression node of a unary operation f( a ), or of a binary operation with a constant operand
   *  @tparam A stored operand type (see Operand)
   *
   *  The derivatives are partial * a.dx( i ), with partial = f'( a ) evaluated once at construction.
   */
  template < typename A >
  class UnaryOperation : public ExpressionBase {
  public:
    /** @brief Number of directions */
    static constexpr std::size_t nDirections = directionsOf< A >;

    /** @brief Creates the node
     *  @tparam E operand type, as deduced by the creating function
     *  @param a operand
     *  @param value value of the operation
     *  @param partial derivative of the operation with respect to a
     */
    template < typename E >
    UnaryOperation( E&& a, double value, double partial )
      : a_( std::forward< E >( a ) ), val_( value ), partial_( partial )
    {
    }

    /** @brief Value of the expression
     *  @return value
     */
    double val() const { return val_; }

    /** @brief Derivative of the expression in one direction
     *  @param i direction
     *  @return partial derivative with respect to independent variable i
     */
    double dx( std::size_t i ) const { return partial_ * a_.dx( i ); }

  private:
    /** @brief Operand */
    A a_;
    /** @brief Value */
    double val_;
    /** @brief Derivative with respect to the operand */
    double partial_;
  };

  /** @brief Expression node of a binary operation f( a, b ) with two non-constant operands
   *  @tparam A stored type of the first operand (see Operand)
   *  @tparam B stored type of the second operand (see Operand)
   *
   *  The derivatives are partialA * a.dx( i ) + partialB * b.dx( i ), with the partial derivatives evaluated once at
   *  construction.
   */
  template < typename A, typename B >
  class BinaryOperation : public ExpressionBase {
  public:
    /** @brief Number of directions */
    static constexpr std::size_t nDirections = directionsOf< A >;
    static_assert( directionsOf< B > == nDirections, "the operands differ in the number of directions" );

    /** @brief Creates the node
     *  @tparam EA first operand type, as deduced by the creating function
     *  @tparam EB second operand type, as deduced by the creating function
     *  @param a first operand
     *  @param b second operand
     *  @param value value of the operation
     *  @param partialA derivative of the operation with respect to a
     *  @param partialB derivative of the operation with respect to b
     */
    template < typename EA, typename EB >
    BinaryOperation( EA&& a, EB&& b, double value, double partialA, double partialB )
      : a_( std::forward< EA >( a ) ),
        b_( std::forward< EB >( b ) ),
        val_( value ),
        partialA_( partialA ),
        partialB_( partialB )
    {
    }

    /** @brief Value of the expression
     *  @return value
     */
    double val() const { return val_; }

    /** @brief Derivative of the expression in one direction
     *  @param i direction
     *  @return partial derivative with respect to independent variable i
     */
    double dx( std::size_t i ) const { return partialA_ * a_.dx( i ) + partialB_ * b_.dx( i ); }

  private:
    /** @brief First operand */
    A a_;
    /** @brief Second operand */
    B b_;
    /** @brief Value */
    double val_;
    /** @brief Derivative with respect to the first operand */
    double partialA_;
    /** @brief Derivative with respect to the second operand */
    double partialB_;
  };

  /** @brief Expression node of a + b or a - b, without multiplications in the derivatives
   *  @tparam A stored type of the first operand (see Operand)
   *  @tparam B stored type of the second operand (see Operand)
   *  @tparam isDifference true for a - b, false for a + b
   */
  template < typename A, typename B, bool isDifference >
  class Sum : public ExpressionBase {
  public:
    /** @brief Number of directions */
    static constexpr std::size_t nDirections = directionsOf< A >;
    static_assert( directionsOf< B > == nDirections, "the operands differ in the number of directions" );

    /** @brief Creates the node
     *  @tparam EA first operand type, as deduced by the creating function
     *  @tparam EB second operand type, as deduced by the creating function
     *  @param a first operand
     *  @param b second operand
     *  @param value a + b, or a - b
     */
    template < typename EA, typename EB >
    Sum( EA&& a, EB&& b, double value ) : a_( std::forward< EA >( a ) ), b_( std::forward< EB >( b ) ), val_( value )
    {
    }

    /** @brief Value of the expression
     *  @return value
     */
    double val() const { return val_; }

    /** @brief Derivative of the expression in one direction
     *  @param i direction
     *  @return partial derivative with respect to independent variable i
     */
    double dx( std::size_t i ) const
    {
      if constexpr ( isDifference )
        return a_.dx( i ) - b_.dx( i );
      else
        return a_.dx( i ) + b_.dx( i );
    }

  private:
    /** @brief First operand */
    A a_;
    /** @brief Second operand */
    B b_;
    /** @brief Value */
    double val_;
  };

  /** @brief Expression node of a + c, a - c (sign +1) or -a, c - a (sign -1) with a constant c
   *  @tparam A stored operand type (see Operand)
   *  @tparam isNegated true if the derivatives are those of -a
   */
  template < typename A, bool isNegated >
  class Shift : public ExpressionBase {
  public:
    /** @brief Number of directions */
    static constexpr std::size_t nDirections = directionsOf< A >;

    /** @brief Creates the node
     *  @tparam E operand type, as deduced by the creating function
     *  @param a operand
     *  @param value value of the operation
     */
    template < typename E >
    Shift( E&& a, double value ) : a_( std::forward< E >( a ) ), val_( value )
    {
    }

    /** @brief Value of the expression
     *  @return value
     */
    double val() const { return val_; }

    /** @brief Derivative of the expression in one direction
     *  @param i direction
     *  @return partial derivative with respect to independent variable i
     */
    double dx( std::size_t i ) const
    {
      if constexpr ( isNegated )
        return -a_.dx( i );
      else
        return a_.dx( i );
    }

  private:
    /** @brief Operand */
    A a_;
    /** @brief Value */
    double val_;
  };

  /** @brief Vector-mode dual number: a value and its partial derivatives with respect to nDirections_ variables
   *  @tparam nDirections_ number of independent variables (directions)
   *
   *  Arithmetic on Dual numbers yields expression nodes, which are evaluated when assigned to a Dual (see the file
   *  documentation). A Dual converts implicitly from arithmetic constants (with zero derivatives) and from expressions.
   */
  template < std::size_t nDirections_ >
  class Dual : public ExpressionBase {
  public:
    /** @brief Number of directions */
    static constexpr std::size_t nDirections = nDirections_;

    /** @brief Creates a zero with zero derivatives */
    Dual() : dx_{}, val_( 0.0 ) {}

    /** @brief Creates a constant, with zero derivatives
     *  @tparam S arithmetic type
     *  @param value value
     */
    template < Constant S >
    Dual( S value ) : dx_{}, val_( static_cast< double >( value ) )
    {
    }

    /** @brief Evaluates an expression
     *  @tparam E expression type
     *  @param e expression
     */
    template < Expression E >
    Dual( const E& e )
    {
      assign( e );
    }

    /** @brief Evaluates an expression into this number; the expression may refer to this number
     *  @tparam E expression type
     *  @param e expression
     *  @return this number
     */
    template < Expression E >
    Dual& operator=( const E& e )
    {
      assign( e );
      return *this;
    }

    /** @brief Sets this number to a constant, with zero derivatives
     *  @tparam S arithmetic type
     *  @param value value
     *  @return this number
     */
    template < Constant S >
    Dual& operator=( S value )
    {
      val_ = static_cast< double >( value );
      dx_.fill( 0.0 );
      return *this;
    }

    /** @brief Value
     *  @return value
     */
    double val() const { return val_; }

    /** @brief Value
     *  @return reference to the value
     */
    double& val() { return val_; }

    /** @brief Partial derivative with respect to one independent variable
     *  @param i direction
     *  @return partial derivative
     */
    double dx( std::size_t i ) const { return dx_[i]; }

    /** @brief Partial derivative with respect to one independent variable
     *  @param i direction
     *  @return reference to the partial derivative
     */
    double& dx( std::size_t i ) { return dx_[i]; }

    /** @brief Adds an expression
     *  @tparam E expression type
     *  @param e expression
     *  @return this number
     */
    template < Expression E >
    Dual& operator+=( const E& e )
    {
      return *this = *this + e;
    }

    /** @brief Subtracts an expression
     *  @tparam E expression type
     *  @param e expression
     *  @return this number
     */
    template < Expression E >
    Dual& operator-=( const E& e )
    {
      return *this = *this - e;
    }

    /** @brief Multiplies by an expression
     *  @tparam E expression type
     *  @param e expression
     *  @return this number
     */
    template < Expression E >
    Dual& operator*=( const E& e )
    {
      return *this = *this * e;
    }

    /** @brief Divides by an expression
     *  @tparam E expression type
     *  @param e expression
     *  @return this number
     */
    template < Expression E >
    Dual& operator/=( const E& e )
    {
      return *this = *this / e;
    }

    /** @brief Adds a constant
     *  @tparam S arithmetic type
     *  @param c constant
     *  @return this number
     */
    template < Constant S >
    Dual& operator+=( S c )
    {
      val_ += static_cast< double >( c );
      return *this;
    }

    /** @brief Subtracts a constant
     *  @tparam S arithmetic type
     *  @param c constant
     *  @return this number
     */
    template < Constant S >
    Dual& operator-=( S c )
    {
      val_ -= static_cast< double >( c );
      return *this;
    }

    /** @brief Multiplies by a constant
     *  @tparam S arithmetic type
     *  @param c constant
     *  @return this number
     */
    template < Constant S >
    Dual& operator*=( S c )
    {
      const double c_ = static_cast< double >( c );
      val_ *= c_;
      for ( double& d : dx_ )
        d *= c_;
      return *this;
    }

    /** @brief Divides by a constant
     *  @tparam S arithmetic type
     *  @param c constant
     *  @return this number
     */
    template < Constant S >
    Dual& operator/=( S c )
    {
      const double c_ = static_cast< double >( c );
      val_ /= c_;
      for ( double& d : dx_ )
        d /= c_;
      return *this;
    }

  private:
    /** @brief Evaluates an expression into this number
     *  @tparam E expression type
     *  @param e expression
     *
     *  The derivatives are collected in a local array first: the compiler then knows that writing them cannot change
     *  the operands of e, which lets it keep the partial derivatives cached in e in registers and vectorize the loop.
     *  This also makes it safe for e to refer to this number.
     */
    template < typename E >
    void assign( const E& e )
    {
      static_assert( E::nDirections == nDirections, "the expression differs in the number of directions" );
      std::array< double, nDirections > dx;
      for ( std::size_t i = 0; i < nDirections; i++ )
        dx[i] = e.dx( i );
      val_ = e.val();
      dx_  = dx;
    }

    /** @brief Partial derivatives */
    std::array< double, nDirections > dx_;
    /** @brief Value */
    double val_;
  };

  // ---------------------------------------------------------------------------------------------------------------
  // node factories
  // ---------------------------------------------------------------------------------------------------------------

  /** @brief Creates a UnaryOperation node
   *  @tparam E operand type
   *  @param a operand
   *  @param value value of the operation
   *  @param partial derivative of the operation with respect to a
   *  @return expression node
   */
  template < Expression E >
  UnaryOperation< Operand< E > > unaryOperation( E&& a, double value, double partial )
  {
    return UnaryOperation< Operand< E > >( std::forward< E >( a ), value, partial );
  }

  /** @brief Creates a BinaryOperation node
   *  @tparam EA first operand type
   *  @tparam EB second operand type
   *  @param a first operand
   *  @param b second operand
   *  @param value value of the operation
   *  @param partialA derivative of the operation with respect to a
   *  @param partialB derivative of the operation with respect to b
   *  @return expression node
   */
  template < Expression EA, Expression EB >
  BinaryOperation< Operand< EA >, Operand< EB > > binaryOperation( EA&&   a,
                                                                   EB&&   b,
                                                                   double value,
                                                                   double partialA,
                                                                   double partialB )
  {
    return BinaryOperation< Operand< EA >, Operand< EB > >( std::forward< EA >( a ),
                                                            std::forward< EB >( b ),
                                                            value,
                                                            partialA,
                                                            partialB );
  }

  // ---------------------------------------------------------------------------------------------------------------
  // arithmetic operators
  // ---------------------------------------------------------------------------------------------------------------

  /** @brief Unary plus
   *  @tparam E expression type
   *  @param a operand
   *  @return a + 0
   */
  template < Expression E >
  Shift< Operand< E >, false > operator+( E&& a )
  {
    const double value = a.val();
    return Shift< Operand< E >, false >( std::forward< E >( a ), value );
  }

  /** @brief Negation
   *  @tparam E expression type
   *  @param a operand
   *  @return -a
   */
  template < Expression E >
  Shift< Operand< E >, true > operator-( E&& a )
  {
    const double value = -a.val();
    return Shift< Operand< E >, true >( std::forward< E >( a ), value );
  }

  /** @brief Sum of two expressions
   *  @tparam EA first operand type
   *  @tparam EB second operand type
   *  @param a first operand
   *  @param b second operand
   *  @return a + b
   */
  template < Expression EA, Expression EB >
  Sum< Operand< EA >, Operand< EB >, false > operator+( EA&& a, EB&& b )
  {
    const double value = a.val() + b.val();
    return Sum< Operand< EA >, Operand< EB >, false >( std::forward< EA >( a ), std::forward< EB >( b ), value );
  }

  /** @brief Sum of an expression and a constant
   *  @tparam E expression type
   *  @tparam S arithmetic type
   *  @param a expression
   *  @param c constant
   *  @return a + c
   */
  template < Expression E, Constant S >
  Shift< Operand< E >, false > operator+( E&& a, S c )
  {
    const double value = a.val() + static_cast< double >( c );
    return Shift< Operand< E >, false >( std::forward< E >( a ), value );
  }

  /** @brief Sum of a constant and an expression
   *  @tparam S arithmetic type
   *  @tparam E expression type
   *  @param c constant
   *  @param a expression
   *  @return c + a
   */
  template < Constant S, Expression E >
  Shift< Operand< E >, false > operator+( S c, E&& a )
  {
    return std::forward< E >( a ) + c;
  }

  /** @brief Difference of two expressions
   *  @tparam EA first operand type
   *  @tparam EB second operand type
   *  @param a first operand
   *  @param b second operand
   *  @return a - b
   */
  template < Expression EA, Expression EB >
  Sum< Operand< EA >, Operand< EB >, true > operator-( EA&& a, EB&& b )
  {
    const double value = a.val() - b.val();
    return Sum< Operand< EA >, Operand< EB >, true >( std::forward< EA >( a ), std::forward< EB >( b ), value );
  }

  /** @brief Difference of an expression and a constant
   *  @tparam E expression type
   *  @tparam S arithmetic type
   *  @param a expression
   *  @param c constant
   *  @return a - c
   */
  template < Expression E, Constant S >
  Shift< Operand< E >, false > operator-( E&& a, S c )
  {
    const double value = a.val() - static_cast< double >( c );
    return Shift< Operand< E >, false >( std::forward< E >( a ), value );
  }

  /** @brief Difference of a constant and an expression
   *  @tparam S arithmetic type
   *  @tparam E expression type
   *  @param c constant
   *  @param a expression
   *  @return c - a
   */
  template < Constant S, Expression E >
  Shift< Operand< E >, true > operator-( S c, E&& a )
  {
    const double value = static_cast< double >( c ) - a.val();
    return Shift< Operand< E >, true >( std::forward< E >( a ), value );
  }

  /** @brief Product of two expressions
   *  @tparam EA first operand type
   *  @tparam EB second operand type
   *  @param a first operand
   *  @param b second operand
   *  @return a * b
   */
  template < Expression EA, Expression EB >
  BinaryOperation< Operand< EA >, Operand< EB > > operator*( EA&& a, EB&& b )
  {
    const double a_ = a.val();
    const double b_ = b.val();
    return binaryOperation( std::forward< EA >( a ), std::forward< EB >( b ), a_ * b_, b_, a_ );
  }

  /** @brief Product of an expression and a constant
   *  @tparam E expression type
   *  @tparam S arithmetic type
   *  @param a expression
   *  @param c constant
   *  @return a * c
   */
  template < Expression E, Constant S >
  UnaryOperation< Operand< E > > operator*( E&& a, S c )
  {
    const double c_    = static_cast< double >( c );
    const double value = a.val() * c_;
    return unaryOperation( std::forward< E >( a ), value, c_ );
  }

  /** @brief Product of a constant and an expression
   *  @tparam S arithmetic type
   *  @tparam E expression type
   *  @param c constant
   *  @param a expression
   *  @return c * a
   */
  template < Constant S, Expression E >
  UnaryOperation< Operand< E > > operator*( S c, E&& a )
  {
    const double c_    = static_cast< double >( c );
    const double value = c_ * a.val();
    return unaryOperation( std::forward< E >( a ), value, c_ );
  }

  /** @brief Quotient of two expressions
   *  @tparam EA first operand type
   *  @tparam EB second operand type
   *  @param a dividend
   *  @param b divisor
   *  @return a / b
   */
  template < Expression EA, Expression EB >
  BinaryOperation< Operand< EA >, Operand< EB > > operator/( EA&& a, EB&& b )
  {
    const double inverse = 1.0 / b.val();
    const double value   = a.val() / b.val();
    return binaryOperation( std::forward< EA >( a ), std::forward< EB >( b ), value, inverse, -value * inverse );
  }

  /** @brief Quotient of an expression and a constant
   *  @tparam E expression type
   *  @tparam S arithmetic type
   *  @param a dividend
   *  @param c divisor
   *  @return a / c
   */
  template < Expression E, Constant S >
  UnaryOperation< Operand< E > > operator/( E&& a, S c )
  {
    const double c_    = static_cast< double >( c );
    const double value = a.val() / c_;
    return unaryOperation( std::forward< E >( a ), value, 1.0 / c_ );
  }

  /** @brief Quotient of a constant and an expression
   *  @tparam S arithmetic type
   *  @tparam E expression type
   *  @param c dividend
   *  @param a divisor
   *  @return c / a
   */
  template < Constant S, Expression E >
  UnaryOperation< Operand< E > > operator/( S c, E&& a )
  {
    const double value = static_cast< double >( c ) / a.val();
    return unaryOperation( std::forward< E >( a ), value, -value / a.val() );
  }

  // ---------------------------------------------------------------------------------------------------------------
  // comparisons (of the values)
  // ---------------------------------------------------------------------------------------------------------------

  /** @brief Value of an expression or a constant, for the comparison operators
   *  @tparam T expression or arithmetic type
   *  @param x expression or constant
   *  @return value of x
   */
  template < typename T >
  double valueOf( const T& x )
  {
    if constexpr ( Expression< T > )
      return x.val();
    else
      return static_cast< double >( x );
  }

  /** @brief An operand pair of a comparison: two expressions, or an expression and a constant
   *  @tparam A first operand type
   *  @tparam B second operand type
   */
  template < typename A, typename B >
  concept Comparable = ( Expression< A > && (Expression< B > || Constant< B >)) || ( Constant< A > && Expression< B > );

  /** @brief Compares the values
   *  @tparam A first operand type
   *  @tparam B second operand type
   *  @param a first operand
   *  @param b second operand
   *  @return value of a == value of b
   */
  template < typename A, typename B >
    requires Comparable< A, B >
  bool operator==( const A& a, const B& b )
  {
    return valueOf( a ) == valueOf( b );
  }

  /** @brief Compares the values
   *  @tparam A first operand type
   *  @tparam B second operand type
   *  @param a first operand
   *  @param b second operand
   *  @return value of a != value of b
   */
  template < typename A, typename B >
    requires Comparable< A, B >
  bool operator!=( const A& a, const B& b )
  {
    return valueOf( a ) != valueOf( b );
  }

  /** @brief Compares the values
   *  @tparam A first operand type
   *  @tparam B second operand type
   *  @param a first operand
   *  @param b second operand
   *  @return value of a < value of b
   */
  template < typename A, typename B >
    requires Comparable< A, B >
  bool operator<( const A& a, const B& b )
  {
    return valueOf( a ) < valueOf( b );
  }

  /** @brief Compares the values
   *  @tparam A first operand type
   *  @tparam B second operand type
   *  @param a first operand
   *  @param b second operand
   *  @return value of a <= value of b
   */
  template < typename A, typename B >
    requires Comparable< A, B >
  bool operator<=( const A& a, const B& b )
  {
    return valueOf( a ) <= valueOf( b );
  }

  /** @brief Compares the values
   *  @tparam A first operand type
   *  @tparam B second operand type
   *  @param a first operand
   *  @param b second operand
   *  @return value of a > value of b
   */
  template < typename A, typename B >
    requires Comparable< A, B >
  bool operator>( const A& a, const B& b )
  {
    return valueOf( a ) > valueOf( b );
  }

  /** @brief Compares the values
   *  @tparam A first operand type
   *  @tparam B second operand type
   *  @param a first operand
   *  @param b second operand
   *  @return value of a >= value of b
   */
  template < typename A, typename B >
    requires Comparable< A, B >
  bool operator>=( const A& a, const B& b )
  {
    return valueOf( a ) >= valueOf( b );
  }

  // ---------------------------------------------------------------------------------------------------------------
  // elementary functions
  // ---------------------------------------------------------------------------------------------------------------

  /** @brief Square root
   *  @tparam E expression type
   *  @param a argument
   *  @return sqrt( a )
   */
  template < Expression E >
  UnaryOperation< Operand< E > > sqrt( E&& a )
  {
    const double value = std::sqrt( a.val() );
    return unaryOperation( std::forward< E >( a ), value, 0.5 / value );
  }

  /** @brief Cube root
   *  @tparam E expression type
   *  @param a argument
   *  @return cbrt( a )
   */
  template < Expression E >
  UnaryOperation< Operand< E > > cbrt( E&& a )
  {
    const double value = std::cbrt( a.val() );
    return unaryOperation( std::forward< E >( a ), value, 1.0 / ( 3.0 * value * value ) );
  }

  /** @brief Exponential function
   *  @tparam E expression type
   *  @param a argument
   *  @return exp( a )
   */
  template < Expression E >
  UnaryOperation< Operand< E > > exp( E&& a )
  {
    const double value = std::exp( a.val() );
    return unaryOperation( std::forward< E >( a ), value, value );
  }

  /** @brief Natural logarithm
   *  @tparam E expression type
   *  @param a argument
   *  @return log( a )
   */
  template < Expression E >
  UnaryOperation< Operand< E > > log( E&& a )
  {
    const double x = a.val();
    return unaryOperation( std::forward< E >( a ), std::log( x ), 1.0 / x );
  }

  /** @brief Decimal logarithm
   *  @tparam E expression type
   *  @param a argument
   *  @return log10( a )
   */
  template < Expression E >
  UnaryOperation< Operand< E > > log10( E&& a )
  {
    const double x = a.val();
    return unaryOperation( std::forward< E >( a ), std::log10( x ), 1.0 / ( x * std::numbers::ln10 ) );
  }

  /** @brief Power with a constant exponent
   *  @tparam E expression type
   *  @tparam S arithmetic type
   *  @param a base
   *  @param exponent exponent
   *  @return pow( a, exponent )
   */
  template < Expression E, Constant S >
  UnaryOperation< Operand< E > > pow( E&& a, S exponent )
  {
    const double x     = a.val();
    const double n     = static_cast< double >( exponent );
    const double value = std::pow( x, n );
    // n x^(n-1), as n x^n / x unless x = 0
    const double partial = n == 0.0 ? 0.0 : x != 0.0 ? n * value / x : n * std::pow( x, n - 1.0 );
    return unaryOperation( std::forward< E >( a ), value, partial );
  }

  /** @brief Power of a constant base
   *  @tparam S arithmetic type
   *  @tparam E expression type
   *  @param base base
   *  @param a exponent
   *  @return pow( base, a )
   */
  template < Constant S, Expression E >
  UnaryOperation< Operand< E > > pow( S base, E&& a )
  {
    const double b     = static_cast< double >( base );
    const double value = std::pow( b, a.val() );
    return unaryOperation( std::forward< E >( a ), value, value == 0.0 ? 0.0 : value * std::log( b ) );
  }

  /** @brief Power of two expressions
   *  @tparam EA base type
   *  @tparam EB exponent type
   *  @param a base
   *  @param b exponent
   *  @return pow( a, b )
   */
  template < Expression EA, Expression EB >
  BinaryOperation< Operand< EA >, Operand< EB > > pow( EA&& a, EB&& b )
  {
    const double x        = a.val();
    const double n        = b.val();
    const double value    = std::pow( x, n );
    const double partialA = n == 0.0 ? 0.0 : x != 0.0 ? n * value / x : n * std::pow( x, n - 1.0 );
    const double partialB = value == 0.0 ? 0.0 : value * std::log( x );
    return binaryOperation( std::forward< EA >( a ), std::forward< EB >( b ), value, partialA, partialB );
  }

  /** @brief Sine
   *  @tparam E expression type
   *  @param a argument
   *  @return sin( a )
   */
  template < Expression E >
  UnaryOperation< Operand< E > > sin( E&& a )
  {
    const double x = a.val();
    return unaryOperation( std::forward< E >( a ), std::sin( x ), std::cos( x ) );
  }

  /** @brief Cosine
   *  @tparam E expression type
   *  @param a argument
   *  @return cos( a )
   */
  template < Expression E >
  UnaryOperation< Operand< E > > cos( E&& a )
  {
    const double x = a.val();
    return unaryOperation( std::forward< E >( a ), std::cos( x ), -std::sin( x ) );
  }

  /** @brief Tangent
   *  @tparam E expression type
   *  @param a argument
   *  @return tan( a )
   */
  template < Expression E >
  UnaryOperation< Operand< E > > tan( E&& a )
  {
    const double value = std::tan( a.val() );
    return unaryOperation( std::forward< E >( a ), value, 1.0 + value * value );
  }

  /** @brief Arc sine
   *  @tparam E expression type
   *  @param a argument
   *  @return asin( a )
   */
  template < Expression E >
  UnaryOperation< Operand< E > > asin( E&& a )
  {
    const double x = a.val();
    return unaryOperation( std::forward< E >( a ), std::asin( x ), 1.0 / std::sqrt( 1.0 - x * x ) );
  }

  /** @brief Arc cosine
   *  @tparam E expression type
   *  @param a argument
   *  @return acos( a )
   */
  template < Expression E >
  UnaryOperation< Operand< E > > acos( E&& a )
  {
    const double x = a.val();
    return unaryOperation( std::forward< E >( a ), std::acos( x ), -1.0 / std::sqrt( 1.0 - x * x ) );
  }

  /** @brief Arc tangent
   *  @tparam E expression type
   *  @param a argument
   *  @return atan( a )
   */
  template < Expression E >
  UnaryOperation< Operand< E > > atan( E&& a )
  {
    const double x = a.val();
    return unaryOperation( std::forward< E >( a ), std::atan( x ), 1.0 / ( 1.0 + x * x ) );
  }

  /** @brief Arc tangent of y / x, using the signs of both to determine the quadrant
   *  @tparam EY type of y
   *  @tparam EX type of x
   *  @param y numerator
   *  @param x denominator
   *  @return atan2( y, x )
   */
  template < Expression EY, Expression EX >
  BinaryOperation< Operand< EY >, Operand< EX > > atan2( EY&& y, EX&& x )
  {
    const double y_ = y.val();
    const double x_ = x.val();
    const double r2 = x_ * x_ + y_ * y_;
    return binaryOperation( std::forward< EY >( y ), std::forward< EX >( x ), std::atan2( y_, x_ ), x_ / r2, -y_ / r2 );
  }

  /** @brief Hyperbolic sine
   *  @tparam E expression type
   *  @param a argument
   *  @return sinh( a )
   */
  template < Expression E >
  UnaryOperation< Operand< E > > sinh( E&& a )
  {
    const double x = a.val();
    return unaryOperation( std::forward< E >( a ), std::sinh( x ), std::cosh( x ) );
  }

  /** @brief Hyperbolic cosine
   *  @tparam E expression type
   *  @param a argument
   *  @return cosh( a )
   */
  template < Expression E >
  UnaryOperation< Operand< E > > cosh( E&& a )
  {
    const double x = a.val();
    return unaryOperation( std::forward< E >( a ), std::cosh( x ), std::sinh( x ) );
  }

  /** @brief Hyperbolic tangent
   *  @tparam E expression type
   *  @param a argument
   *  @return tanh( a )
   */
  template < Expression E >
  UnaryOperation< Operand< E > > tanh( E&& a )
  {
    const double value = std::tanh( a.val() );
    return unaryOperation( std::forward< E >( a ), value, 1.0 - value * value );
  }

  /** @brief Absolute value; its derivative at 0 is taken as 1
   *  @tparam E expression type
   *  @param a argument
   *  @return |a|
   */
  template < Expression E >
  UnaryOperation< Operand< E > > abs( E&& a )
  {
    const double x = a.val();
    return unaryOperation( std::forward< E >( a ), std::abs( x ), x < 0.0 ? -1.0 : 1.0 );
  }

  /** @brief Absolute value; its derivative at 0 is taken as 1
   *  @tparam E expression type
   *  @param a argument
   *  @return |a|
   */
  template < Expression E >
  UnaryOperation< Operand< E > > fabs( E&& a )
  {
    return abs( std::forward< E >( a ) );
  }

  /** @brief Larger of two numbers; for equal values, the derivatives of a
   *  @tparam nDirections number of directions
   *  @param a first number
   *  @param b second number
   *  @return a copy of the larger number
   */
  template < std::size_t nDirections >
  Dual< nDirections > max( const Dual< nDirections >& a, const Dual< nDirections >& b )
  {
    return a.val() < b.val() ? b : a;
  }

  /** @brief Smaller of two numbers; for equal values, the derivatives of a
   *  @tparam nDirections number of directions
   *  @param a first number
   *  @param b second number
   *  @return a copy of the smaller number
   */
  template < std::size_t nDirections >
  Dual< nDirections > min( const Dual< nDirections >& a, const Dual< nDirections >& b )
  {
    return b.val() < a.val() ? b : a;
  }

  /** @brief True if the value is not a number
   *  @tparam E expression type
   *  @param a expression
   *  @return std::isnan of the value
   */
  template < Expression E >
  bool isnan( const E& a )
  {
    return std::isnan( a.val() );
  }

  /** @brief True if the value is infinite
   *  @tparam E expression type
   *  @param a expression
   *  @return std::isinf of the value
   */
  template < Expression E >
  bool isinf( const E& a )
  {
    return std::isinf( a.val() );
  }

  /** @brief True if the value is finite
   *  @tparam E expression type
   *  @param a expression
   *  @return std::isfinite of the value
   */
  template < Expression E >
  bool isfinite( const E& a )
  {
    return std::isfinite( a.val() );
  }

  /** @brief Writes the value and the partial derivatives
   *  @tparam E expression type
   *  @param out stream
   *  @param a expression
   *  @return out
   */
  template < Expression E >
  std::ostream& operator<<( std::ostream& out, const E& a )
  {
    out << a.val() << " [";
    for ( std::size_t i = 0; i < directionsOf< E >; i++ )
      out << ( i ? ", " : "" ) << a.dx( i );
    return out << "]";
  }

} // namespace Marmot::Math::Differentiation::Automatic::VectorMode

namespace Eigen {
  /** @brief Eigen support for vector-mode dual numbers
   *  @tparam nDirections number of directions
   */
  template < std::size_t nDirections >
  struct NumTraits< Marmot::Math::Differentiation::Automatic::VectorMode::Dual< nDirections > > : NumTraits< double > {
    /** @brief Real type */
    using Real = Marmot::Math::Differentiation::Automatic::VectorMode::Dual< nDirections >;
    /** @brief Non-integer type */
    using NonInteger = Real;
    /** @brief Literal type */
    using Literal = double;
    /** @brief Nested type */
    using Nested = Real;
    /** @brief Properties */
    enum {
      IsComplex             = 0,
      IsInteger             = 0,
      IsSigned              = 1,
      RequireInitialization = 1,
      ReadCost              = 1,
      AddCost               = 3,
      MulCost               = 3
    };
  };

  /** @brief Lets Eigen mix vector-mode dual numbers with double
   *  @tparam nDirections number of directions
   *  @tparam BinaryOp Eigen binary functor
   */
  template < std::size_t nDirections, typename BinaryOp >
  struct ScalarBinaryOpTraits< Marmot::Math::Differentiation::Automatic::VectorMode::Dual< nDirections >,
                               double,
                               BinaryOp > {
    /** @brief Result type */
    using ReturnType = Marmot::Math::Differentiation::Automatic::VectorMode::Dual< nDirections >;
  };

  /** @brief Lets Eigen mix double with vector-mode dual numbers
   *  @tparam nDirections number of directions
   *  @tparam BinaryOp Eigen binary functor
   */
  template < std::size_t nDirections, typename BinaryOp >
  struct ScalarBinaryOpTraits< double,
                               Marmot::Math::Differentiation::Automatic::VectorMode::Dual< nDirections >,
                               BinaryOp > {
    /** @brief Result type */
    using ReturnType = Marmot::Math::Differentiation::Automatic::VectorMode::Dual< nDirections >;
  };
} // namespace Eigen
