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
#pragma once
#include "Marmot/BeamFiberSection.h"
#include "Marmot/MarmotElement.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotFiniteElement.h"
#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <string>
#include <vector>

namespace Marmot::Discretizations::Elements {

  /**
   * @class Marmot::Discretizations::Elements::BeamElement
   * @brief Euler-Bernoulli beam (frame) element in two or three spatial dimensions with 2 or 3 nodes, geometrically
   * linear, with a fiber section (BeamFiberSection) for any hypoelastic material of Marmot.
   *
   * @tparam nDim   number of spatial dimensions: 2 (nodal dofs \f$u_x, u_y, \theta_z\f$) or 3 (nodal dofs
   *                \f$u_x, u_y, u_z, \theta_x, \theta_y, \theta_z\f$)
   * @tparam nNodes 2 (linear) or 3 (quadratic; node order end, end, mid as Bar3; the mid node must lie at the middle
   *                of the straight beam)
   *
   * Node fields: @b displacement and @b rotation (2D: one rotation about the out-of-plane axis; 3D: the rotation
   * vector). The dofs are ordered node by node, displacement before rotation.
   *
   * @par Kinematics
   * In the local frame \f$( \boldsymbol{e}_1, \boldsymbol{e}_2, \boldsymbol{e}_3 )\f$ (\f$\boldsymbol{e}_1\f$ along
   * the axis from node 1 to node 2, section coordinates \f$y, z\f$) the axial displacement \f$u\f$ and the twist
   * \f$\varphi\f$ are interpolated by the Lagrange polynomials of the nodes (linear or quadratic), the transverse
   * displacements \f$v, w\f$ by the Hermite polynomials of the nodal values and the nodal rotations
   * \f$\theta_z = v'\f$ and \f$\theta_y = - w'\f$ (cubic for 2 nodes, quintic for 3 nodes; Euler-Bernoulli: no shear
   * deformation, no shear locking). The generalized strains are \f$\varepsilon_0 = u'\f$, \f$\kappa_y = -w''\f$,
   * \f$\kappa_z = v''\f$ and \f$\chi = \varphi'\f$ (2D: \f$\varepsilon_0\f$ and \f$\kappa = \kappa_z\f$). For a linear
   * elastic material, the 2-node element reproduces the exact nodal values of a beam loaded at its nodes, the 3-node
   * element also the exact deflection of a beam under a uniform line load (quartic).
   *
   * Integration along the axis (full, exact for elastic beams): 3 Gauss points (2 nodes), 4 Gauss points (3 nodes).
   *
   * Geometrically linear: rotations must remain small. The section (BeamFiberSection) is independent of the
   * kinematics, so a co-rotational variant only has to provide the generalized strains of its co-rotated frame and
   * the corresponding (geometric) tangent.
   *
   * @par Section and material
   * Every integration point along the axis has its own BeamFiberSection; every fiber has
   * its own material instance (created by name through the hypoelastic material factory from the
   * Marmot::Discretizations::MaterialSection) and state. 2D fibers are in uniaxial stress
   * (Marmot::Materials::MaterialHypoElastic::computeUniaxialStress), 3D fibers additionally carry the torsional shear
   * strains (Marmot::Materials::MaterialHypoElastic::computeBeamStress).
   *
   * Element properties: the section is given by its integration points (any shape; the element knows no shapes):
   * - 2D: \f$[ n, y_1, A_1, \dots, y_n, A_n ]\f$, \f$y\f$ in the plane, normal to the axis (\f$\boldsymbol{e}_2 =
   *   \boldsymbol{e}_z \times \boldsymbol{e}_1\f$).
   * - 3D: \f$[ v_x, v_y, v_z, J, n, y_1, z_1, A_1, \dots, y_n, z_n, A_n ]\f$ with the orientation vector
   *   \f$\boldsymbol{v}\f$, whose part normal to the axis is the local \f$y\f$ axis \f$\boldsymbol{e}_2\f$
   *   (\f$\boldsymbol{e}_3 = \boldsymbol{e}_1 \times \boldsymbol{e}_2\f$), and the torsion constant \f$J\f$.
   *
   * The section coordinates are measured from the beam axis (the line through the nodes): a section whose centroid is
   * off the axis is an eccentric beam. The materials see the characteristic length \f$L/3\f$.
   *
   * Quadrature point states: 2D: @b normal force, @b bending moment, @b axial strain, @b curvature; 3D: @b normal
   * force, @b bending moments (\f$M_y, M_z\f$), @b torque, @b axial strain, @b curvatures (\f$\kappa_y, \kappa_z\f$),
   * @b twist rate; both: @b section forces and @b section strains (all components), @b elastic energy, @b dissipation
   * (integrated over the quadrature point's volume). The state of fiber \f$k\f$ is accessible as "fiber k stress"
   * (\f$\sigma_{11}, \sigma_{12}, \sigma_{13}\f$) or "fiber k <material state name>".
   *
   * Not implemented: distributed (surface) loads, inertia.
   */
  template < int nDim, int nNodes >
  class BeamElement : public Marmot::Discretizations::Elements::Element {

    static_assert( nDim == 2 || nDim == 3, "BeamElement: nDim must be 2 or 3" );
    static_assert( nNodes == 2 || nNodes == 3, "BeamElement: nNodes must be 2 or 3" );

  public:
    using Section = BeamFiberSection< nDim >;

    static constexpr int nHermite       = 2 * nNodes; ///< transverse dofs (value and slope per node)
    static constexpr int nRotations     = nDim == 2 ? 1 : 3;
    static constexpr int nDofPerNode    = nDim + nRotations;
    static constexpr int sizeLoadVector = nNodes * nDofPerNode;
    static constexpr int nGeneralized   = Section::nGeneralized;
    static constexpr int nHeader        = nDim == 2 ? 1 : 5; ///< properties before the section points
    static constexpr int pointStride    = nDim == 2 ? 2 : 3; ///< properties per section point

    using VectorDim  = Eigen::Matrix< double, nDim, 1 >;
    using MatrixDim  = Eigen::Matrix< double, nDim, nDim >;
    using RhsSized   = Eigen::Matrix< double, sizeLoadVector, 1 >;
    using KSized     = Eigen::Matrix< double, sizeLoadVector, sizeLoadVector >;
    using BSized     = Eigen::Matrix< double, nGeneralized, sizeLoadVector >;
    using VectorGen  = typename Section::VectorGen;
    using MatrixGen  = typename Section::MatrixGen;
    using InterpSize = Eigen::Matrix< double, nDim, sizeLoadVector >;

    /// A quadrature point along the axis, with its own fiber section
    struct QuadraturePoint {
      const double xi;        ///< parent coordinate in [-1, 1]
      const double weight;    ///< Gauss weight
      BSized       B;         ///< generalized strains of the global dofs
      double       dV = 0.0;  ///< length x weight / 2
      double*      stateVars; ///< the state of this point
      Section      section;   ///< the fibers with their materials
      QuadraturePoint( double xi, double weight ) : xi( xi ), weight( weight ), stateVars( nullptr ) {}
    };

    /// Element label
    const int elLabel;

    /// Element properties, see the class documentation
    Eigen::Map< const Eigen::VectorXd > elementProperties;
    /// Reference nodal coordinates, [dim, node]
    Eigen::Map< const Eigen::Matrix< double, nDim, nNodes > > coordinates;

    /// Length of the element
    double length = 0.0;
    /// Rows: the local axes e1, e2(, e3) in global coordinates
    MatrixDim R = MatrixDim::Identity();
    /// Global -> local dof transformation
    KSized T = KSized::Identity();
    /// The quadrature points along the axis
    std::vector< QuadraturePoint > qps;

    /// The material section, kept to create the fiber materials once the fibers exist
    std::string   materialName;
    const double* materialProperties  = nullptr;
    int           nMaterialProperties = 0;

    /// The names of the per-point states and their lengths, see the class documentation
    static const std::vector< std::pair< std::string, int > >& qpStateEntries()
    {
      static const std::vector< std::pair< std::string, int > >
        entries = nDim == 2 ? std::vector< std::pair< std::string, int > >{ { "normal force", 1 },
                                                                            { "bending moment", 1 },
                                                                            { "axial strain", 1 },
                                                                            { "curvature", 1 },
                                                                            { "elastic energy", 1 },
                                                                            { "dissipation", 1 } }
                            : std::vector< std::pair< std::string, int > >{ { "normal force", 1 },
                                                                            { "bending moments", 2 },
                                                                            { "torque", 1 },
                                                                            { "axial strain", 1 },
                                                                            { "curvatures", 2 },
                                                                            { "twist rate", 1 },
                                                                            { "elastic energy", 1 },
                                                                            { "dissipation", 1 } };
      return entries;
    }

    /// offsets in the state of a point
    static constexpr int idxSectionForces  = 0;
    static constexpr int idxSectionStrains = nGeneralized;
    static constexpr int idxElasticEnergy  = 2 * nGeneralized;
    static constexpr int idxDissipation    = 2 * nGeneralized + 1;
    static constexpr int nQpStateVarsOwn   = 2 * nGeneralized + 2;

    /// The parametric coordinates s in [0, 1] of the nodes (end, end, mid)
    static constexpr std::array< double, 3 > nodeS = { 0.0, 1.0, 0.5 };

    /// The coefficients of the Hermite polynomials in the monomials s^k: column 2a for the value at node a, column 2a+1
    /// for the slope dv/ds at node a
    Eigen::Matrix< double, nHermite, nHermite > hermiteCoefficients;

    explicit BeamElement( int elementID )
      : elLabel( elementID ), elementProperties( nullptr, 0 ), coordinates( nullptr )
    {
      using namespace Marmot::Discretizations::Elements::Quadrature;
      if constexpr ( nNodes == 2 ) {
        for ( const auto& qpInfo : Marmot::Discretizations::Elements::Quadrature::Spatial1D::gaussPointList3 )
          qps.emplace_back( qpInfo.xi( 0 ), qpInfo.weight );
      }
      else {
        // 4 point Gauss rule
        const double a  = std::sqrt( 3. / 7 - 2. / 7 * std::sqrt( 6. / 5 ) ),
                     b  = std::sqrt( 3. / 7 + 2. / 7 * std::sqrt( 6. / 5 ) );
        const double wa = ( 18 + std::sqrt( 30. ) ) / 36, wb = ( 18 - std::sqrt( 30. ) ) / 36;
        for ( const auto& [xi, w] :
              std::array< std::pair< double, double >, 4 >{ { { -b, wb }, { -a, wa }, { a, wa }, { b, wb } } } )
          qps.emplace_back( xi, w );
      }

      Eigen::Matrix< double, nHermite, nHermite > M = Eigen::Matrix< double, nHermite, nHermite >::Zero();
      for ( int a = 0; a < nNodes; a++ )
        for ( int k = 0; k < nHermite; k++ ) {
          M( 2 * a, k )     = std::pow( nodeS[a], k );
          M( 2 * a + 1, k ) = k == 0 ? 0.0 : k * std::pow( nodeS[a], k - 1 );
        }
      hermiteCoefficients = M.inverse();
    }

    /// The Hermite polynomials (row) or their second derivative w.r.t. s at s, ordered as the columns of
    /// hermiteCoefficients
    Eigen::Matrix< double, 1, nHermite > hermite( double s, int derivative ) const
    {
      Eigen::Matrix< double, 1, nHermite > monomials = Eigen::Matrix< double, 1, nHermite >::Zero();
      for ( int k = derivative; k < nHermite; k++ )
        monomials( k ) = ( derivative == 0 ? 1.0 : k * ( k - 1 ) ) * std::pow( s, k - derivative );
      return monomials * hermiteCoefficients;
    }

    /// The Lagrange shape functions and their derivatives w.r.t. s at s (node order end, end, mid)
    static std::pair< Eigen::Matrix< double, 1, nNodes >, Eigen::Matrix< double, 1, nNodes > > lagrange( double s )
    {
      const double xi = 2 * s - 1;
      if constexpr ( nNodes == 2 )
        return { Marmot::Discretizations::Elements::Spatial1D::Bar2::N( xi ),
                 2 * Marmot::Discretizations::Elements::Spatial1D::Bar2::dNdXi( xi ) };
      else
        return { Marmot::Discretizations::Elements::Spatial1D::Bar3::N( xi ),
                 2 * Marmot::Discretizations::Elements::Spatial1D::Bar3::dNdXi( xi ) };
    }

    int nQpStateVars() const { return nQpStateVarsOwn + qps[0].section.getNumberOfRequiredStateVars(); }

    int getNumberOfRequiredStateVars() override { return static_cast< int >( qps.size() ) * nQpStateVars(); }

    std::vector< std::vector< std::string > > getNodeFields() override
    {
      return std::vector< std::vector< std::string > >( nNodes, { "displacement", "rotation" } );
    }

    std::vector< int > getDofIndicesPermutationPattern() override
    {
      std::vector< int > pattern( sizeLoadVector );
      for ( int i = 0; i < sizeLoadVector; i++ )
        pattern[i] = i;
      return pattern;
    }

    int getNNodes() override { return nNodes; }

    int getNSpatialDimensions() override { return nDim; }

    int getNDofPerElement() override { return sizeLoadVector; }

    std::string getElementShape() override { return nNodes == 2 ? "bar2" : "bar3"; }

    void assignStateVars( double* stateVars, int nStateVars ) override
    {
      if ( nStateVars < getNumberOfRequiredStateVars() )
        throw std::invalid_argument( Marmot::Utilities::MakeString()
                                     << __PRETTY_FUNCTION__ << ": too few state variables" );
      for ( size_t i = 0; i < qps.size(); i++ )
        qps[i].stateVars = stateVars + i * nQpStateVars();
    }

    void assignProperty( const ElementProperties& property ) override
    {
      const int n = property.nElementProperties >= nHeader
                      ? static_cast< int >( property.elementProperties[nHeader - 1] )
                      : -1;
      if ( n < 1 || property.nElementProperties != nHeader + pointStride * n )
        throw std::invalid_argument(
          Marmot::Utilities::MakeString()
          << __PRETTY_FUNCTION__ << ": element " << elLabel << ": the properties of a beam in " << nDim << "D are "
          << ( nDim == 2 ? "[n, y_1, A_1, ...]" : "[vx, vy, vz, J, n, y_1, z_1, A_1, ...]" ) );
      new ( &elementProperties )
        Eigen::Map< const Eigen::VectorXd >( property.elementProperties, property.nElementProperties );
    }

    void assignProperty( const Marmot::Discretizations::MaterialSection& section ) override
    {
      materialName        = section.materialName;
      materialProperties  = section.materialProperties;
      nMaterialProperties = section.nMaterialProperties;
      createMaterials();
    }

    /// one material per fiber, as soon as both the fibers (initializeYourself) and the material section are known
    void createMaterials()
    {
      if ( materialName.empty() || qps[0].section.fibers.empty() )
        return;
      for ( auto& qp : qps )
        qp.section.createMaterials( materialName,
                                    materialProperties,
                                    nMaterialProperties,
                                    elLabel,
                                    length / static_cast< double >( qps.size() ) );
    }

    void assignNodeCoordinates( const double* coords ) override
    {
      new ( &coordinates ) Eigen::Map< const Eigen::Matrix< double, nDim, nNodes > >( coords );
    }

    /// The number of section points
    int nSectionPoints() const { return static_cast< int >( elementProperties[nHeader - 1] ); }

    /// The cross section area (the sum of the areas of the section points)
    double area() const { return qps[0].section.integrals()[0]; }

    /// The local displacements [u, v(, w)] at the parent coordinate xi in terms of the global dofs
    InterpSize displacementInterpolation( double xi ) const
    {
      const double s                               = 0.5 * ( xi + 1 );
      const auto [N, dNds]                         = lagrange( s );
      const Eigen::Matrix< double, 1, nHermite > H = hermite( s, 0 );

      InterpSize N_ = InterpSize::Zero();
      for ( int a = 0; a < nNodes; a++ ) {
        const int    i  = a * nDofPerNode;
        const double Hv = H( 2 * a ), Ht = length * H( 2 * a + 1 ); // the slope dof is dv/dx = dv/ds / L
        N_( 0, i ) = N( a );
        if constexpr ( nDim == 2 ) {
          N_( 1, i + 1 ) = Hv, N_( 1, i + 2 ) = Ht;
        }
        else {
          N_( 1, i + 1 ) = Hv, N_( 1, i + 5 ) = Ht;  // v with theta_z = v'
          N_( 2, i + 2 ) = Hv, N_( 2, i + 4 ) = -Ht; // w with theta_y = -w'
        }
      }
      return N_ * T;
    }

    /// The generalized strains at the parent coordinate xi in terms of the global dofs
    BSized strainOperator( double xi ) const
    {
      const double s                                 = 0.5 * ( xi + 1 );
      const double L                                 = length;
      const auto [N, dNds]                           = lagrange( s );
      const Eigen::Matrix< double, 1, nHermite > d2H = hermite( s, 2 );

      BSized B = BSized::Zero();
      for ( int a = 0; a < nNodes; a++ ) {
        const int    i  = a * nDofPerNode;
        const double dN = dNds( a ) / L;                                       // d/dx
        const double dv = d2H( 2 * a ) / ( L * L ), dt = d2H( 2 * a + 1 ) / L; // d2/dx2 of the value and slope parts
        B( 0, i ) = dN;
        if constexpr ( nDim == 2 ) {
          B( 1, i + 1 ) = dv, B( 1, i + 2 ) = dt; // kappa = v''
        }
        else {
          B( 1, i + 2 ) = -dv, B( 1, i + 4 ) = dt; // kappa_y = -w''
          B( 2, i + 1 ) = dv, B( 2, i + 5 ) = dt;  // kappa_z = v''
          B( 3, i + 3 ) = dN;                      // chi = phi'
        }
      }
      return B * T;
    }

    void initializeYourself() override
    {
      const VectorDim axis = coordinates.col( 1 ) - coordinates.col( 0 );
      length               = axis.norm();
      if ( length <= 0 )
        throw std::invalid_argument( Marmot::Utilities::MakeString()
                                     << __PRETTY_FUNCTION__ << ": element " << elLabel << " has zero length" );
      const VectorDim e1 = axis / length;
      if constexpr ( nNodes == 3 ) {
        if ( ( coordinates.col( 2 ) - 0.5 * ( coordinates.col( 0 ) + coordinates.col( 1 ) ) ).norm() > 1e-6 * length )
          throw std::invalid_argument( Marmot::Utilities::MakeString()
                                       << __PRETTY_FUNCTION__ << ": element " << elLabel
                                       << ": the mid node must lie at the middle of the beam" );
      }
      if constexpr ( nDim == 2 ) {
        R.row( 0 ) = e1.transpose();
        R.row( 1 ) << -e1( 1 ), e1( 0 );
      }
      else {
        const Eigen::Vector3d v  = elementProperties.template segment< 3 >( 0 );
        const Eigen::Vector3d vn = v - v.dot( e1 ) * e1;
        if ( v.norm() == 0 || vn.norm() < 1e-8 * v.norm() )
          throw std::invalid_argument( Marmot::Utilities::MakeString()
                                       << __PRETTY_FUNCTION__ << ": element " << elLabel
                                       << ": the orientation vector must not be parallel to the axis" );
        const Eigen::Vector3d e2 = vn.normalized();
        R.row( 0 )               = e1.transpose();
        R.row( 1 )               = e2.transpose();
        R.row( 2 )               = e1.cross( e2 ).transpose();
      }

      T.setZero();
      for ( int a = 0; a < nNodes; a++ ) {
        T.template block< nDim, nDim >( a * nDofPerNode, a * nDofPerNode ) = R;
        if constexpr ( nDim == 2 )
          T( a * nDofPerNode + 2, a * nDofPerNode + 2 ) = 1;
        else
          T.template block< 3, 3 >( a * nDofPerNode + 3, a * nDofPerNode + 3 ) = R;
      }

      const double J = nDim == 3 ? elementProperties[3] : 0.0;
      for ( auto& qp : qps ) {
        qp.B  = strainOperator( qp.xi );
        qp.dV = 0.5 * length * qp.weight;
        try {
          qp.section.setFibers( elementProperties.data() + nHeader, nSectionPoints(), J );
        }
        catch ( const std::invalid_argument& e ) {
          throw std::invalid_argument( Marmot::Utilities::MakeString() << "element " << elLabel << ": " << e.what() );
        }
      }
      createMaterials();
    }

    void setInitialConditions( StateTypes state, const double* ) override
    {
      switch ( state ) {
      case Marmot::Discretizations::Elements::Element::MarmotMaterialInitialization: {
        for ( auto& qp : qps )
          qp.section.initializeStateVars( qp.stateVars + nQpStateVarsOwn );
        break;
      }
      default:
        throw std::invalid_argument( Marmot::Utilities::MakeString()
                                     << __PRETTY_FUNCTION__ << ": invalid initial condition" );
      }
    }

    void computeKernels( const double* QTotal_, const double* dQ_, double* Pe_, double* Ke_, double time, double dT )
      override
    {
      Eigen::Map< const RhsSized > QTotal( QTotal_ );
      Eigen::Map< const RhsSized > dQ( dQ_ );
      Eigen::Map< RhsSized >       Pe( Pe_ );
      Eigen::Map< KSized >         Ke( Ke_ );

      for ( auto& qp : qps ) {
        VectorGen s;
        MatrixGen D;
        double    energy, dissipation;
        qp.section.computeSection( qp.stateVars + nQpStateVarsOwn, s, D, energy, dissipation, qp.B * dQ, time, dT );

        Eigen::Map< VectorGen >( qp.stateVars + idxSectionForces )  = s;
        Eigen::Map< VectorGen >( qp.stateVars + idxSectionStrains ) = qp.B * QTotal;
        qp.stateVars[idxElasticEnergy]                              = energy * qp.dV;
        qp.stateVars[idxDissipation]                                = dissipation * qp.dV;

        Pe += qp.B.transpose() * s * qp.dV;
        Ke += qp.B.transpose() * D * qp.B * qp.dV;
      }
    }

    void computeKernelsExplicit( const double* QTotal, const double* dQ, double* Pe, double time, double dT ) override
    {
      KSized K = KSized::Zero();
      computeKernels( QTotal, dQ, Pe, K.data(), time, dT );
    }

    void computeDistributedLoad( DistributedLoadTypes,
                                 double*,
                                 double*,
                                 int,
                                 const double*,
                                 const double*,
                                 double,
                                 double ) override
    {
      throw std::invalid_argument( Marmot::Utilities::MakeString()
                                   << __PRETTY_FUNCTION__ << ": distributed loads are not supported by beam elements" );
    }

    /// Consistent nodal forces and moments of a body force (per volume), i.e., a line load b A along the axis
    void computeBodyForce( double* P_, double*, const double* load, const double*, double, double ) override
    {
      Eigen::Map< RhsSized >        P( P_ );
      Eigen::Map< const VectorDim > b( load );
      const VectorDim               q = R * b * area(); // local line load
      for ( const auto& qp : qps )
        P += displacementInterpolation( qp.xi ).transpose() * q * qp.dV;
    }

    void computeInternalEnergy( double& internalEnergy ) override
    {
      for ( const auto& qp : qps )
        internalEnergy += qp.stateVars[idxElasticEnergy];
    }

    Marmot::Utilities::StateView getStateView( const std::string& stateName, int qpNumber ) override
    {
      const auto& qp = qps[qpNumber];
      if ( stateName == "section forces" )
        return { qp.stateVars + idxSectionForces, nGeneralized };
      if ( stateName == "section strains" )
        return { qp.stateVars + idxSectionStrains, nGeneralized };
      int offset = 0;
      for ( const auto& [name, len] : qpStateEntries() ) {
        if ( name == stateName )
          return { qp.stateVars + offset, len };
        offset += len;
      }
      // "fiber <k> <state>"
      if ( stateName.rfind( "fiber ", 0 ) == 0 ) {
        const size_t space = stateName.find( ' ', 6 );
        if ( space != std::string::npos )
          return qp.section.getFiberStateView( qp.stateVars + nQpStateVarsOwn,
                                               std::stoul( stateName.substr( 6, space - 6 ) ),
                                               stateName.substr( space + 1 ) );
      }
      throw std::invalid_argument( Marmot::Utilities::MakeString()
                                   << __PRETTY_FUNCTION__ << ": unknown state " << stateName );
    }

    std::vector< double > getCoordinatesAtCenter() override
    {
      const VectorDim x = 0.5 * ( coordinates.col( 0 ) + coordinates.col( 1 ) );
      return std::vector< double >( x.data(), x.data() + nDim );
    }

    std::vector< std::vector< double > > getCoordinatesAtQuadraturePoints() override
    {
      std::vector< std::vector< double > > result;
      for ( const auto& qp : qps ) {
        const double    s = 0.5 * ( qp.xi + 1 );
        const VectorDim x = ( 1 - s ) * coordinates.col( 0 ) + s * coordinates.col( 1 );
        result.emplace_back( x.data(), x.data() + nDim );
      }
      return result;
    }

    int getNumberOfQuadraturePoints() override { return static_cast< int >( qps.size() ); }
  };

} // namespace Marmot::Discretizations::Elements
