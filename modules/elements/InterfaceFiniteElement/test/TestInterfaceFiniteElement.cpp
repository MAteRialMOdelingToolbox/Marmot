#include "Marmot/InterfaceFiniteElement.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotTesting.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

using namespace Marmot;
using namespace Marmot::Testing;
using namespace Marmot::Elements;

namespace {

  /**
   * Records the characteristic element length seen at every computeStress call, then delegates.
   * The element holds ONE material for all quadrature points, so the length of the current point
   * has to be installed before each evaluation; this material makes that observable.
   */
  class LengthRecordingInterfaceMaterial : public MarmotInterfaceMaterialHypoElastic {
  public:
    LengthRecordingInterfaceMaterial()
      : MarmotInterfaceMaterialHypoElastic( "LINEARELASTIC", materialProperties.data(), 3, 0 )
    {
    }

    mutable std::vector< double > seenLengths;

    void computeStress( State&               state,
                        Tangents&            tangents,
                        const Deformation&   deformation,
                        const TimeIncrement& timeIncrement ) override
    {
      seenLengths.push_back( characteristicElementLength );
      MarmotInterfaceMaterialHypoElastic::computeStress( state, tangents, deformation, timeIncrement );
    }

  private:
    inline static const std::array< double, 3 > materialProperties = { 4000.0, 0.3, 0.01 };
  };

  class FailingInterfaceMaterial : public MarmotInterfaceMaterialHypoElastic {
  public:
    FailingInterfaceMaterial() : MarmotInterfaceMaterialHypoElastic( "LINEARELASTIC", materialProperties.data(), 3, 0 )
    {
    }

    void computeStress( State&, Tangents&, const Deformation&, const TimeIncrement& ) override
    {
      throw Marmot::StressUpdateFailed( "Deliberate interface-material update failure." );
    }

  private:
    inline static const std::array< double, 3 > materialProperties = { 1.0, 0.0, 1.0 };
  };

  template < typename DerivedA, typename DerivedB >
  void assertMatrixNear( const Eigen::MatrixBase< DerivedA >& actual,
                         const Eigen::MatrixBase< DerivedB >& expected,
                         double                               tol,
                         const std::string&                   message )
  {
    throwExceptionOnFailure( actual.rows() == expected.rows() && actual.cols() == expected.cols(),
                             message + ": matrix shape mismatch." );

    const double err = ( actual - expected ).template lpNorm< Eigen::Infinity >();
    throwExceptionOnFailure( err < tol, message + ": max error = " + std::to_string( err ) );
  }

  std::unique_ptr< InterfaceFiniteElement< 3, 8 > > makeSingleInputFileInterfaceElement()
  {
    constexpr int nDim   = 3;
    constexpr int nNodes = 8;

    const int  elId    = 3;
    const auto intType = FiniteElement::Quadrature::IntegrationTypes::FullIntegration;
    const auto secType = InterfaceFiniteElement< nDim, nNodes >::SectionType::Interface;

    auto element = std::make_unique< InterfaceFiniteElement< nDim, nNodes > >( elId, intType, secType );

    static std::array< double, nDim* nNodes > coordinates = {
      -0.500000, -0.500000, 0.088163,  0.500000,  -0.500000, -0.088163,
      0.500000,  0.500000,  -0.088163, -0.500000, 0.500000,  0.088163,

      -0.500000, -0.500000, 0.188163,  0.500000,  -0.500000, 0.011837,
      0.500000,  0.500000,  0.011837,  -0.500000, 0.500000,  0.188163,
    };
    element->assignNodeCoordinates( coordinates.data() );

    static std::array< double, 1 > elPropsVec = { 1.0 };
    ElementProperties              elProps( elPropsVec.data(), static_cast< int >( elPropsVec.size() ) );
    element->assignProperty( elProps );

    static std::array< double, 8 > materialProperties = { 4000.0, 0.30, 0.01, 1.06e-3, 0.334, 12, 1.0e-2, 1 };
    element->assignMaterial( "LINEARVISCOELASTICPOWERLAW", materialProperties.data(), materialProperties.size() );
    return element;
  }

  std::unique_ptr< InterfaceFiniteElement< 2, 4 > > makeTwoDimensionalInterfaceElement()
  {
    constexpr int nDim   = 2;
    constexpr int nNodes = 4;

    const int  elId    = 2;
    const auto intType = FiniteElement::Quadrature::IntegrationTypes::FullIntegration;
    const auto secType = InterfaceFiniteElement< nDim, nNodes >::SectionType::Interface;

    auto element = std::make_unique< InterfaceFiniteElement< nDim, nNodes > >( elId, intType, secType );

    static std::array< double, nDim* nNodes > coordinates = {
      0.0,
      0.0,
      1.0,
      0.0,
      0.0,
      0.0,
      1.0,
      0.0,
    };
    element->assignNodeCoordinates( coordinates.data() );

    static std::array< double, 1 > elPropsVec = { 1.0 };
    ElementProperties              elProps( elPropsVec.data(), static_cast< int >( elPropsVec.size() ) );
    element->assignProperty( elProps );

    static std::array< double, 3 > materialProperties = { 1000.0, 0.25, 0.1 };
    element->assignMaterial( "LINEARELASTIC", materialProperties.data(), materialProperties.size() );

    return element;
  }

  template < int nDim, int nNodes >
  void initializeStateAndMaterial( InterfaceFiniteElement< nDim, nNodes >& element, std::vector< double >& stateVars )
  {
    stateVars.assign( element.getNumberOfRequiredStateVars(), 0.0 );
    element.assignStateVars( stateVars.data(), static_cast< int >( stateVars.size() ) );
    element.initializeYourself();
    element.setInitialConditions( MarmotElement::MarmotMaterialInitialization, nullptr );
  }

  void TestMaterialInitializationResetsMaterialState()
  {
    constexpr int nDim   = 3;
    constexpr int nNodes = 8;

    const auto intType = FiniteElement::Quadrature::IntegrationTypes::FullIntegration;
    auto       element = std::make_unique< InterfaceFiniteElement< nDim, nNodes > >( 3, intType );

    static std::array< double, nDim* nNodes > coordinates = {
      -0.5, -0.5, 0.0, 0.5, -0.5, 0.0, 0.5, 0.5, 0.0, -0.5, 0.5, 0.0,
      -0.5, -0.5, 0.1, 0.5, -0.5, 0.1, 0.5, 0.5, 0.1, -0.5, 0.5, 0.1,
    };
    element->assignNodeCoordinates( coordinates.data() );

    static std::array< double, 1 > elementProperties = { 1.0 };
    ElementProperties              properties( elementProperties.data(), elementProperties.size() );
    element->assignProperty( properties );

    static std::array< double, 8 > materialProperties = { 4000.0, 0.30, 0.01, 1.06e-3, 0.334, 12, 1.0e-2, 1 };
    element->assignMaterial( "LINEARVISCOELASTICPOWERLAW", materialProperties.data(), materialProperties.size() );

    std::vector< double > stateVars( element->getNumberOfRequiredStateVars(), 0.0 );
    element->assignStateVars( stateVars.data(), stateVars.size() );

    for ( auto& qp : element->qps ) {
      qp.managedStateVars->materialStateVars.setOnes();
    }

    element->setInitialConditions( MarmotElement::MarmotMaterialInitialization, nullptr );

    for ( const auto& qp : element->qps ) {
      throwExceptionOnFailure( qp.managedStateVars->materialStateVars.isZero(),
                               "MarmotMaterialInitialization did not reset interface-material state." );
    }
  }

  void TestUnsupportedInertiaThrows()
  {
    constexpr int nDim         = 3;
    constexpr int nNodes       = 8;
    constexpr int nElementDofs = nDim * nNodes;

    const auto intType = FiniteElement::Quadrature::IntegrationTypes::FullIntegration;
    auto       element = std::make_unique< InterfaceFiniteElement< nDim, nNodes > >( 3, intType );

    std::array< double, nElementDofs * nElementDofs > consistentInertia{};
    std::array< double, nElementDofs >                lumpedInertia{};

    bool consistentInertiaRejected = false;
    try {
      element->computeConsistentInertia( consistentInertia.data() );
    }
    catch ( const std::runtime_error& ) {
      consistentInertiaRejected = true;
    }

    bool lumpedInertiaRejected = false;
    try {
      element->computeLumpedInertia( lumpedInertia.data() );
    }
    catch ( const std::runtime_error& ) {
      lumpedInertiaRejected = true;
    }

    throwExceptionOnFailure( consistentInertiaRejected,
                             "InterfaceFiniteElement did not reject consistent inertia computation." );
    throwExceptionOnFailure( lumpedInertiaRejected,
                             "InterfaceFiniteElement did not reject lumped inertia computation." );
  }

  void TestStressUpdateFailureRequestsSmallerTimeStep()
  {
    auto element = makeSingleInputFileInterfaceElement();

    std::vector< double > stateVars( element->getNumberOfRequiredStateVars(), 0.0 );
    element->assignStateVars( stateVars.data(), stateVars.size() );
    element->initializeYourself();

    element->material = std::make_unique< FailingInterfaceMaterial >();

    constexpr int                                     nElementDofs = 3 * 8;
    std::array< double, nElementDofs >                QTotal{};
    std::array< double, nElementDofs >                dQ{};
    std::array< double, nElementDofs >                Pe{};
    std::array< double, nElementDofs * nElementDofs > Ke{};
    const std::array< double, 2 >                     time               = { 0.0, 0.0 };
    bool                                              stressUpdateFailed = false;
    try {
      element->computeKernels( QTotal.data(), dQ.data(), Pe.data(), Ke.data(), time[0], 1.0 );
    }
    catch ( const Marmot::StressUpdateFailed& ) {
      stressUpdateFailed = true;
    }

    throwExceptionOnFailure( stressUpdateFailed, "InterfaceFiniteElement did not propagate StressUpdateFailed." );
  }

  template < typename QuadraturePointType >
  double integrationWeight( const QuadraturePointType& qp )
  {
    return qp.weight * qp.sqrtDetG;
  }

  struct ExpectedResponse {
    std::vector< Eigen::VectorXd > qpForce;
    std::vector< Eigen::VectorXd > qpSurfaceStress;
    std::vector< Eigen::VectorXd > qpResidual;
    Eigen::VectorXd                Pe;
  };

  ExpectedResponse computeExpectedResponseFromGaussPoints( InterfaceFiniteElement< 3, 8 >& element,
                                                           const Eigen::VectorXd&          dU,
                                                           const double*                   time,
                                                           double                          dT )
  {
    constexpr int nDim      = 3;
    constexpr int nNodes    = 8;
    constexpr int totalNDof = nDim * nNodes;
    constexpr int halfNDof  = totalNDof / 2;

    ExpectedResponse expected;
    expected.Pe = Eigen::VectorXd::Zero( totalNDof );

    for ( int q = 0; q < element.getNumberOfQuadraturePoints(); ++q ) {
      auto&       qp    = element.qps[q];
      const auto& Nside = qp.NmatSide;
      const auto& Bside = qp.BmatSide;
      const auto& Njump = qp.NmatJump;
      const auto& Bavg  = qp.BmatAverage;

      const int nTensor = static_cast< int >( Bside.rows() );

      Eigen::VectorXd dUGp( 2 * nDim );
      dUGp.segment( 0, nDim )    = Nside * dU.segment( halfNDof, halfNDof );
      dUGp.segment( nDim, nDim ) = Nside * dU.segment( 0, halfNDof );

      Eigen::VectorXd dSurfaceStrainGp( 2 * nTensor );
      dSurfaceStrainGp.segment( 0, nTensor )       = Bside * dU.segment( halfNDof, halfNDof );
      dSurfaceStrainGp.segment( nTensor, nTensor ) = Bside * dU.segment( 0, halfNDof );

      auto force         = qp.managedStateVars->force;
      auto surfaceStress = qp.managedStateVars->surfaceStress;

      Eigen::MatrixXd Qij = Eigen::MatrixXd::Zero( nDim, nDim );
      Eigen::MatrixXd Z   = Eigen::MatrixXd::Zero( nTensor, nTensor );
      Eigen::MatrixXd H   = Eigen::MatrixXd::Zero( nDim, nTensor );
      Eigen::MatrixXd Y   = Eigen::MatrixXd::Zero( nTensor, nTensor );

      MarmotInterfaceMaterialHypoElastic::State         state{ force.data(),
                                                       surfaceStress.data(),
                                                       qp.managedStateVars->materialStateVars.data() };
      MarmotInterfaceMaterialHypoElastic::Tangents      tangents{ Qij.data(), Z.data(), H.data(), Y.data() };
      MarmotInterfaceMaterialHypoElastic::Deformation   deformation{ dUGp.data(),
                                                                   dSurfaceStrainGp.data(),
                                                                   qp.normal.data() };
      MarmotInterfaceMaterialHypoElastic::TimeIncrement timeIncrement{ time[0], dT };
      element.material->computeStress( state, tangents, deformation, timeIncrement );

      const double J0xW = integrationWeight( qp );

      throwExceptionOnFailure( std::isfinite( J0xW ) && J0xW > 0.0,
                               "Invalid computed integration weight at qp " + std::to_string( q ) );
      throwExceptionOnFailure( force.allFinite(),
                               "Direct material force contains nan or inf at qp " + std::to_string( q ) );
      throwExceptionOnFailure( surfaceStress.allFinite(),
                               "Direct material surface stress contains nan or inf at qp " + std::to_string( q ) );
      throwExceptionOnFailure( Qij.allFinite(),
                               "Direct material Q tangent contains nan or inf at qp " + std::to_string( q ) );
      throwExceptionOnFailure( Z.allFinite(),
                               "Direct material Z tangent contains nan or inf at qp " + std::to_string( q ) );
      throwExceptionOnFailure( H.allFinite(),
                               "Direct material H tangent contains nan or inf at qp " + std::to_string( q ) );
      throwExceptionOnFailure( Y.allFinite(),
                               "Direct material Y tangent contains nan or inf at qp " + std::to_string( q ) );

      const Eigen::VectorXd Pqp = Njump.transpose() * force * J0xW + Bavg.transpose() * surfaceStress * J0xW;

      expected.qpForce.emplace_back( force );
      expected.qpSurfaceStress.emplace_back( surfaceStress );
      expected.qpResidual.emplace_back( Pqp );
      expected.Pe += Pqp;
    }

    return expected;
  }

  Eigen::VectorXd makeTopSideXOpeningIncrement()
  {
    constexpr int nDim      = 3;
    constexpr int nNodes    = 8;
    constexpr int totalNDof = nDim * nNodes;
    constexpr int halfNDof  = totalNDof / 2;

    Eigen::VectorXd dU = Eigen::VectorXd::Zero( totalNDof );

    dU( halfNDof + 3 * 1 + 0 ) = 1e-2;
    dU( halfNDof + 3 * 2 + 0 ) = 1e-2;

    return dU;
  }

} // namespace

void TestSingleInputFileElementGeometryMatrices()
{
  std::cout << "\n--- TestSingleInputFileElementGeometryMatrices ---\n";

  constexpr int nDim      = 3;
  constexpr int nNodes    = 8;
  constexpr int totalNDof = nDim * nNodes;
  constexpr int halfNDof  = totalNDof / 2;

  auto element = makeSingleInputFileInterfaceElement();
  element->initializeYourself();

  throwExceptionOnFailure( element->getNumberOfQuadraturePoints() == 4,
                           "IQuad4 full integration should have 4 Gauss points." );
  throwExceptionOnFailure( element->getNDofPerElement() == totalNDof, "Unexpected number of interface element dofs." );

  const double tol = 1e-12;

  for ( int q = 0; q < element->getNumberOfQuadraturePoints(); ++q ) {
    const auto& qp = element->qps[q];

    throwExceptionOnFailure( qp.NmatSide.rows() == nDim && qp.NmatSide.cols() == halfNDof,
                             "Unexpected NmatSide dimensions." );
    throwExceptionOnFailure( qp.NmatJump.rows() == nDim && qp.NmatJump.cols() == totalNDof,
                             "Unexpected NmatJump dimensions." );
    throwExceptionOnFailure( qp.BmatAverage.cols() == totalNDof, "Unexpected BmatAverage column count." );
    throwExceptionOnFailure( std::abs( qp.normal.norm() - 1.0 ) < tol, "Interface normal is not normalized." );

    const double J0xW = integrationWeight( qp );

    throwExceptionOnFailure( std::isfinite( J0xW ) && J0xW > 0.0,
                             "Computed Gauss-point integration weight must be positive and finite." );
    throwExceptionOnFailure( std::isfinite( qp.J0xW ) && std::abs( qp.J0xW - J0xW ) < 1e-12,
                             "Stored qp.J0xW differs from weight * sqrtDetG. Check element property lifetime." );

    const double xi  = qp.xi( 0 );
    const double eta = qp.xi( 1 );
    const double N1  = 0.25 * ( 1.0 - xi ) * ( 1.0 - eta );
    const double N2  = 0.25 * ( 1.0 + xi ) * ( 1.0 - eta );
    const double N3  = 0.25 * ( 1.0 + xi ) * ( 1.0 + eta );
    const double N4  = 0.25 * ( 1.0 - xi ) * ( 1.0 + eta );

    Eigen::Matrix< double, nDim, halfNDof > NsideExpected;
    NsideExpected << N1, 0.0, 0.0, N2, 0.0, 0.0, N3, 0.0, 0.0, N4, 0.0, 0.0, 0.0, N1, 0.0, 0.0, N2, 0.0, 0.0, N3, 0.0,
      0.0, N4, 0.0, 0.0, 0.0, N1, 0.0, 0.0, N2, 0.0, 0.0, N3, 0.0, 0.0, N4;

    assertMatrixNear( qp.NmatSide, NsideExpected, tol, "NmatSide differs from the IQuad4 shape-function matrix." );

    Eigen::Matrix< double, nDim, totalNDof > NjumpExpected;
    NjumpExpected.setZero();
    NjumpExpected.template block< nDim, halfNDof >( 0, 0 )        = -NsideExpected;
    NjumpExpected.template block< nDim, halfNDof >( 0, halfNDof ) = NsideExpected;

    assertMatrixNear( qp.NmatJump, NjumpExpected, tol, "NmatJump differs from top-minus-bottom convention." );

    Eigen::Matrix< double, halfNDof, 1 > constantSideDisplacement;
    constantSideDisplacement << 0.123, -0.456, 0.789, 0.123, -0.456, 0.789, 0.123, -0.456, 0.789, 0.123, -0.456, 0.789;

    throwExceptionOnFailure( ( qp.BmatSide * constantSideDisplacement ).template lpNorm< Eigen::Infinity >() < tol,
                             "BmatSide must annihilate constant side displacements." );
  }
}

void TestSurfaceGradientOperatorProjectsOnlyGradientDirection()
{
  std::cout << "\n--- TestSurfaceGradientOperatorProjectsOnlyGradientDirection ---\n";

  const double tol = 1e-12;

  // BSurfaceMatrix must build grad_s u = grad(u) . T: the GRADIENT DIRECTION is projected onto the
  // tangent plane, the DISPLACEMENT COMPONENT is not. That keeps the in-plane derivative of the normal
  // displacement, which MarmotInterfaceMaterialHypoElastic needs for the thin-layer gradient
  // grad(u) = (1/h) [u] (x) n + <grad_s u>. The same formula must hold for every nDim, so that a
  // plane-strain model and its 3D extrusion see the same surface strain.
  auto checkOperator = [&]( const auto& element, int nDim, int nInterfaceNodes, const std::string& label ) {
    for ( const auto& qp : element.qps ) {

      bool differsFromFullyProjected = false;

      for ( int A = 0; A < nInterfaceNodes; ++A ) {
        for ( int k = 0; k < nDim; ++k ) {

          double expected = 0.0;
          for ( int j = 0; j < nDim; ++j )
            expected += qp.gradN( j, A ) * qp.tangentProjection( j, k );

          for ( int i = 0; i < nDim; ++i ) {
            for ( int m = 0; m < nDim; ++m ) {

              const double actual = qp.BmatSide( i * nDim + k, A * nDim + m );
              const double wanted = ( i == m ) ? expected : 0.0;

              throwExceptionOnFailure( std::abs( actual - wanted ) < tol,
                                       label +
                                         ": BmatSide must project the gradient direction only. "
                                         "Entry (" +
                                         std::to_string( i * nDim + k ) + "," + std::to_string( A * nDim + m ) +
                                         ") is " + std::to_string( actual ) + ", expected " + std::to_string( wanted ) +
                                         "." );

              double fullyProjected = 0.0;
              for ( int j = 0; j < nDim; ++j )
                fullyProjected += qp.tangentProjection( i, m ) * qp.gradN( j, A ) * qp.tangentProjection( j, k );

              differsFromFullyProjected = differsFromFullyProjected || std::abs( fullyProjected - wanted ) > tol;
            }
          }
        }
      }

      // Guard the regression: on this geometry the fully projected operator T(i,m) gradN(j,A) T(j,k) is a
      // genuinely different matrix, so reintroducing it would fail the assertions above rather than pass
      // vacuously.
      throwExceptionOnFailure( differsFromFullyProjected,
                               label + ": test geometry must distinguish the projected-component operator." );
    }
  };

  {
    auto element = makeSingleInputFileInterfaceElement();
    element->initializeYourself();
    checkOperator( *element, 3, 4, "3D interface" );
  }

  {
    // 45 degree interface: the tangent is not axis aligned, so a purely normal nodal displacement has a
    // nonzero in-plane derivative and the two candidate operators disagree.
    constexpr int nDim   = 2;
    constexpr int nNodes = 4;

    auto element = std::make_unique<
      InterfaceFiniteElement< nDim, nNodes > >( 1,
                                                FiniteElement::Quadrature::IntegrationTypes::FullIntegration,
                                                InterfaceFiniteElement< nDim, nNodes >::SectionType::Interface );

    static std::array< double, nDim* nNodes > nodeCoordsVec = { 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0 };
    element->assignNodeCoordinates( nodeCoordsVec.data() );

    static std::array< double, 1 > elPropsVec = { 1.0 };
    ElementProperties              elProps( elPropsVec.data(), static_cast< int >( elPropsVec.size() ) );
    element->assignProperty( elProps );

    element->initializeYourself();
    checkOperator( *element, 2, 2, "2D angled interface" );
  }
}

void TestSingleInputFileElementMaterialResponseIsFinite()
{
  std::cout << "\n--- TestSingleInputFileElementMaterialResponseIsFinite ---\n";

  auto element = makeSingleInputFileInterfaceElement();

  std::vector< double > stateVars;
  initializeStateAndMaterial( *element, stateVars );

  const Eigen::VectorXd dU       = makeTopSideXOpeningIncrement();
  double                time     = 0.0;
  double                dT       = 0.1;
  const auto            expected = computeExpectedResponseFromGaussPoints( *element, dU, &time, dT );

  throwExceptionOnFailure( expected.Pe.allFinite(), "Material-only expected residual contains nan or inf." );
  throwExceptionOnFailure( expected.Pe.template lpNorm< Eigen::Infinity >() > 0.0,
                           "The top-side x opening should produce a nonzero material residual." );
}

void TestSingleInputFileElementGaussPointStiffnessAndResidual()
{
  std::cout << "\n--- TestSingleInputFileElementGaussPointStiffnessAndResidual ---\n";

  constexpr int nDim      = 3;
  constexpr int nNodes    = 8;
  constexpr int totalNDof = nDim * nNodes;

  auto expectedElement = makeSingleInputFileInterfaceElement();
  auto testedElement   = makeSingleInputFileInterfaceElement();

  std::vector< double > expectedStateVars;
  std::vector< double > testedStateVars;

  initializeStateAndMaterial( *expectedElement, expectedStateVars );
  initializeStateAndMaterial( *testedElement, testedStateVars );

  const Eigen::VectorXd dUEigen = makeTopSideXOpeningIncrement();

  std::vector< double > U( totalNDof, 0.0 );
  std::vector< double > dU( dUEigen.data(), dUEigen.data() + dUEigen.size() );
  std::vector< double > Pe( totalNDof, 0.0 );
  std::vector< double > Ke( totalNDof * totalNDof, 0.0 );

  double timeForExpected = 0.0;
  double timeForElement  = 0.0;
  double dT              = 0.1;

  const auto expected = computeExpectedResponseFromGaussPoints( *expectedElement, dUEigen, &timeForExpected, dT );

  testedElement->computeKernels( U.data(), dU.data(), Pe.data(), Ke.data(), timeForElement, dT );

  Eigen::Map< Eigen::Matrix< double, totalNDof, 1 > >                          PeActual( Pe.data() );
  Eigen::Map< Eigen::Matrix< double, totalNDof, totalNDof, Eigen::RowMajor > > KeActual( Ke.data() );

  const double tolResidual = 1e-9;

  throwExceptionOnFailure( PeActual.allFinite(), "Element residual contains nan or inf." );
  throwExceptionOnFailure( KeActual.allFinite(), "Element stiffness contains nan or inf." );
  throwExceptionOnFailure( PeActual.template lpNorm< Eigen::Infinity >() > 0.0,
                           "The chosen top-side x opening should produce a nonzero residual." );
  throwExceptionOnFailure( KeActual.template lpNorm< Eigen::Infinity >() > 0.0,
                           "The linear elastic interface should produce a nonzero stiffness." );

  for ( int q = 0; q < testedElement->getNumberOfQuadraturePoints(); ++q ) {
    const auto& qp = testedElement->qps[q];

    throwExceptionOnFailure( expected.qpResidual[q].allFinite(),
                             "Expected Gauss-point residual contains nan or inf at qp " + std::to_string( q ) );
    throwExceptionOnFailure( expected.qpResidual[q].template lpNorm< Eigen::Infinity >() > 0.0,
                             "The chosen increment should produce a nonzero Gauss-point residual at qp " +
                               std::to_string( q ) );

    assertMatrixNear( qp.managedStateVars->force,
                      expected.qpForce[q],
                      tolResidual,
                      "Gauss-point force differs from direct material evaluation at qp " + std::to_string( q ) );
    assertMatrixNear( qp.managedStateVars->surfaceStress,
                      expected.qpSurfaceStress[q],
                      tolResidual,
                      "Gauss-point surface stress differs from direct material evaluation at qp " +
                        std::to_string( q ) );
  }

  assertMatrixNear( PeActual,
                    expected.Pe,
                    tolResidual,
                    "Assembled residual differs from sum of Gauss-point residuals." );

  auto computePeForIncrement = [&]( const Eigen::VectorXd& dUIncrement ) {
    auto element = makeSingleInputFileInterfaceElement();

    std::vector< double > stateVars;
    initializeStateAndMaterial( *element, stateVars );

    std::vector< double > Ulocal( totalNDof, 0.0 );
    std::vector< double > dUlocal( dUIncrement.data(), dUIncrement.data() + dUIncrement.size() );
    std::vector< double > Pelocal( totalNDof, 0.0 );
    std::vector< double > Kelocal( totalNDof * totalNDof, 0.0 );

    double timeLocal = 0.0;
    double dTLocal   = 0.1;

    element->computeKernels( Ulocal.data(), dUlocal.data(), Pelocal.data(), Kelocal.data(), timeLocal, dTLocal );

    return Eigen::Map< Eigen::Matrix< double, totalNDof, 1 > >( Pelocal.data() ).eval();
  };

  Eigen::Matrix< double, totalNDof, totalNDof > KeFiniteDifference;
  const double                                  eps = 1e-7;

  for ( int j = 0; j < totalNDof; ++j ) {
    Eigen::VectorXd dUPlus  = dUEigen;
    Eigen::VectorXd dUMinus = dUEigen;

    dUPlus( j ) += eps;
    dUMinus( j ) -= eps;

    KeFiniteDifference.col( j ) = ( computePeForIncrement( dUPlus ) - computePeForIncrement( dUMinus ) ) /
                                  ( 2.0 * eps );
  }

  const double err = ( KeActual - KeFiniteDifference ).template lpNorm< Eigen::Infinity >();

  const double fdNorm = KeFiniteDifference.template lpNorm< Eigen::Infinity >();

  const double relErr = err / std::max( 1.0, fdNorm );

  Eigen::Index maxRow         = 0;
  Eigen::Index maxCol         = 0;
  const double maxAbsMismatch = ( KeActual - KeFiniteDifference ).cwiseAbs().maxCoeff( &maxRow, &maxCol );

  std::cout << "max |Ke - KeFD| = " << maxAbsMismatch << " at (" << maxRow << ", " << maxCol << ")\n";
  std::cout << "KeActual(" << maxRow << "," << maxCol << ") = " << KeActual( maxRow, maxCol ) << "\n";
  std::cout << "KeFiniteDifference(" << maxRow << "," << maxCol << ") = " << KeFiniteDifference( maxRow, maxCol )
            << "\n";

  std::cout << "KeActual row " << maxRow << ": " << KeActual.row( maxRow ) << "\n";
  std::cout << "KeFD row " << maxRow << ": " << KeFiniteDifference.row( maxRow ) << "\n";
  std::cout << "KeActual - KeFD row " << maxRow << ": " << ( KeActual - KeFiniteDifference ).row( maxRow ) << "\n";

  std::cout << "stiffness finite-difference check: "
            << "relative error for Ke = +dPe/ddU is " << relErr << "\n";

  throwExceptionOnFailure( relErr < 1e-5, "Assembled stiffness is not consistent with Ke = +dPe/ddU." );
}

void TestSingleInputFileElementRigidTranslationGivesZeroResidual()
{
  std::cout << "\n--- TestSingleInputFileElementRigidTranslationGivesZeroResidual ---\n";

  constexpr int nDim      = 3;
  constexpr int nNodes    = 8;
  constexpr int totalNDof = nDim * nNodes;

  auto element = makeSingleInputFileInterfaceElement();

  std::vector< double > stateVars;
  initializeStateAndMaterial( *element, stateVars );

  std::vector< double > U( totalNDof, 0.0 );
  std::vector< double > dU( totalNDof, 0.0 );
  std::vector< double > Pe( totalNDof, 0.0 );
  std::vector< double > Ke( totalNDof * totalNDof, 0.0 );

  for ( int a = 0; a < nNodes; ++a ) {
    dU[3 * a + 0] = 0.01;
    dU[3 * a + 1] = -0.02;
    dU[3 * a + 2] = 0.03;
  }

  double time = 0.0;
  double dT   = 0.1;

  element->computeKernels( U.data(), dU.data(), Pe.data(), Ke.data(), time, dT );

  Eigen::Map< Eigen::Matrix< double, totalNDof, 1 > > PeActual( Pe.data() );

  throwExceptionOnFailure( PeActual.allFinite(), "Rigid-translation residual contains nan or inf." );
  throwExceptionOnFailure( PeActual.template lpNorm< Eigen::Infinity >() < 1e-9,
                           "Equal rigid translation on both sides should produce zero interface residual." );
}

void TestAssignStateVarsPreservesHistoryAcrossIncrements()
{
  std::cout << "\n--- TestAssignStateVarsPreservesHistoryAcrossIncrements ---\n";

  constexpr int nDim      = 3;
  constexpr int nNodes    = 8;
  constexpr int totalNDof = nDim * nNodes;

  auto element = makeSingleInputFileInterfaceElement();

  std::vector< double > stateVars;
  initializeStateAndMaterial( *element, stateVars ); // proper one-time init: step1/inc1 equivalent

  std::vector< double > U( totalNDof, 0.0 );
  std::vector< double > dU( totalNDof, 0.0 );
  std::vector< double > Pe( totalNDof, 0.0 );
  std::vector< double > Ke( totalNDof * totalNDof, 0.0 );

  const Eigen::VectorXd dUEigen = makeTopSideXOpeningIncrement();
  for ( int i = 0; i < totalNDof; ++i )
    dU[i] = dUEigen[i];

  double time = 0.0;
  double dT   = 1.0;

  // simulate increment 1: some nonzero deformation drives material state away from zero
  element->computeKernels( U.data(), dU.data(), Pe.data(), Ke.data(), time, dT );

  const Eigen::VectorXd stateAfterInc1 = element->qps[0].managedStateVars->materialStateVars;

  throwExceptionOnFailure( stateAfterInc1.template lpNorm< Eigen::Infinity >() > 0.0,
                           "Precondition failed: material state should be nonzero after a real increment." );

  // simulate what Abaqus/UEL does on every subsequent call: re-point state vars at the
  // same persistent buffer, WITHOUT requesting MarmotMaterialInitialization again
  element->assignStateVars( stateVars.data(), static_cast< int >( stateVars.size() ) );

  const Eigen::VectorXd stateAfterReassign = element->qps[0].managedStateVars->materialStateVars;

  throwExceptionOnFailure( stateAfterReassign.isApprox( stateAfterInc1 ),
                           "assignStateVars wiped accumulated material history — "
                           "this breaks history-dependent materials (e.g. creep) across increments!" );
}

void TestTwoDimensionalInterfaceElementComputesWithEmbeddedMaterial()
{
  std::cout << "\n--- TestTwoDimensionalInterfaceElementComputesWithEmbeddedMaterial ---\n";

  constexpr int nDim      = 2;
  constexpr int nNodes    = 4;
  constexpr int totalNDof = nDim * nNodes;
  constexpr int halfNDof  = totalNDof / 2;
  constexpr int nTensor   = nDim * nDim;

  auto element = makeTwoDimensionalInterfaceElement();

  std::vector< double > stateVars;
  initializeStateAndMaterial( *element, stateVars );

  Eigen::Matrix< double, totalNDof, 1 >         U;
  Eigen::Matrix< double, totalNDof, 1 >         dU;
  Eigen::Matrix< double, totalNDof, 1 >         Pe;
  Eigen::Matrix< double, totalNDof, totalNDof > Ke;

  U.setZero();
  dU.setZero();
  Pe.setZero();
  Ke.setZero();

  dU << 0.0, -1.0e-5, 0.0, 2.0e-5, 1.0e-5, 1.2e-4, -2.0e-5, 8.0e-5;

  Eigen::Matrix< double, totalNDof, 1 >         expectedPe;
  Eigen::Matrix< double, totalNDof, totalNDof > expectedKe;
  expectedPe.setZero();
  expectedKe.setZero();

  std::vector< Eigen::Matrix< double, nDim, 1 > >    expectedForce;
  std::vector< Eigen::Matrix< double, nTensor, 1 > > expectedSurfaceStress;
  const double                                       time = 0.0;
  const double                                       dT   = 1.0;

  for ( auto& qp : element->qps ) {
    const auto& Nside = qp.NmatSide;
    const auto& Bside = qp.BmatSide;
    const auto& Njump = qp.NmatJump;
    const auto& Bavg  = qp.BmatAverage;

    Eigen::Matrix< double, 2 * nDim, 1 > dUGp;
    dUGp.template segment< nDim >( 0 )    = Nside * dU.template segment< halfNDof >( halfNDof );
    dUGp.template segment< nDim >( nDim ) = Nside * dU.template segment< halfNDof >( 0 );

    Eigen::Matrix< double, 2 * nTensor, 1 > dSurfaceStrainGp;
    dSurfaceStrainGp.template segment< nTensor >( 0 )       = Bside * dU.template segment< halfNDof >( halfNDof );
    dSurfaceStrainGp.template segment< nTensor >( nTensor ) = Bside * dU.template segment< halfNDof >( 0 );

    Eigen::Vector3d                                force3d = Eigen::Vector3d::Zero();
    Eigen::Matrix< double, 9, 1 >                  surfaceStress3d;
    Eigen::Matrix< double, 6, 1 >                  dU3d;
    Eigen::Matrix< double, 18, 1 >                 dSurfaceStrain3d;
    Eigen::Vector3d                                normal3d = Eigen::Vector3d::Zero();
    Eigen::Matrix< double, 3, 3, Eigen::RowMajor > Q3d;
    Eigen::Matrix< double, 9, 9, Eigen::RowMajor > Z3d;
    Eigen::Matrix< double, 3, 9, Eigen::RowMajor > H3d;
    Eigen::Matrix< double, 9, 9, Eigen::RowMajor > Y3d;

    surfaceStress3d.setZero();
    dU3d.setZero();
    dSurfaceStrain3d.setZero();
    Q3d.setZero();
    Z3d.setZero();
    H3d.setZero();
    Y3d.setZero();

    for ( int i = 0; i < nDim; ++i ) {
      normal3d( i ) = qp.normal( i );
      dU3d( i )     = dUGp( i );
      dU3d( 3 + i ) = dUGp( nDim + i );

      for ( int j = 0; j < nDim; ++j ) {
        const int index2d = i * nDim + j;
        const int index3d = i * 3 + j;

        dSurfaceStrain3d( index3d )     = dSurfaceStrainGp( index2d );
        dSurfaceStrain3d( 9 + index3d ) = dSurfaceStrainGp( nTensor + index2d );
      }
    }

    std::vector< double > materialStateVars( qp.managedStateVars->materialStateVars.size(), 0.0 );

    MarmotInterfaceMaterialHypoElastic::State       materialState{ force3d.data(),
                                                             surfaceStress3d.data(),
                                                             materialStateVars.data() };
    MarmotInterfaceMaterialHypoElastic::Tangents    materialTangents{ Q3d.data(), Z3d.data(), H3d.data(), Y3d.data() };
    MarmotInterfaceMaterialHypoElastic::Deformation materialDeformation{ dU3d.data(),
                                                                         dSurfaceStrain3d.data(),
                                                                         normal3d.data() };
    MarmotInterfaceMaterialHypoElastic::TimeIncrement materialTimeIncrement{ time, dT };
    element->material->computeStress( materialState, materialTangents, materialDeformation, materialTimeIncrement );

    Eigen::Matrix< double, nDim, 1 >                           force2d;
    Eigen::Matrix< double, nTensor, 1 >                        surfaceStress2d;
    Eigen::Matrix< double, nDim, nDim, Eigen::RowMajor >       Q2d;
    Eigen::Matrix< double, nTensor, nTensor, Eigen::RowMajor > Z2d;
    Eigen::Matrix< double, nDim, nTensor, Eigen::RowMajor >    H2d;
    Eigen::Matrix< double, nTensor, nTensor, Eigen::RowMajor > Y2d;

    force2d.setZero();
    surfaceStress2d.setZero();
    Q2d.setZero();
    Z2d.setZero();
    H2d.setZero();
    Y2d.setZero();

    for ( int i = 0; i < nDim; ++i ) {
      force2d( i ) = force3d( i );

      for ( int j = 0; j < nDim; ++j ) {
        const int index2d = i * nDim + j;
        const int index3d = i * 3 + j;

        surfaceStress2d( index2d ) = surfaceStress3d( index3d );
        Q2d( i, j )                = Q3d( i, j );

        for ( int k = 0; k < nDim; ++k ) {
          const int tensorCol2d = j * nDim + k;
          const int tensorCol3d = j * 3 + k;

          H2d( i, tensorCol2d ) = H3d( i, tensorCol3d );

          for ( int l = 0; l < nDim; ++l ) {
            const int tensorRow2d  = i * nDim + j;
            const int tensorRow3d  = i * 3 + j;
            const int tensorCol2d4 = k * nDim + l;
            const int tensorCol3d4 = k * 3 + l;

            Z2d( tensorRow2d, tensorCol2d4 ) = Z3d( tensorRow3d, tensorCol3d4 );
            Y2d( tensorRow2d, tensorCol2d4 ) = Y3d( tensorRow3d, tensorCol3d4 );
          }
        }
      }
    }

    expectedForce.emplace_back( force2d );
    expectedSurfaceStress.emplace_back( surfaceStress2d );

    expectedPe += Njump.transpose() * force2d * qp.J0xW;
    expectedPe += Bavg.transpose() * surfaceStress2d * qp.J0xW;

    expectedKe += ( Njump.transpose() * Q2d * Njump + Bavg.transpose() * Z2d * Bavg + Bavg.transpose() * Y2d * Bavg +
                    Njump.transpose() * H2d * Bavg + Bavg.transpose() * H2d.transpose() * Njump ) *
                  qp.J0xW;
  }

  element->computeKernels( U.data(), dU.data(), Pe.data(), Ke.data(), time, dT );

  assertMatrixNear( Pe, expectedPe, 1e-10, "2D interface element residual differs from embedded material response." );
  assertMatrixNear( Ke, expectedKe, 1e-10, "2D interface element tangent differs from embedded material response." );

  for ( int q = 0; q < element->getNumberOfQuadraturePoints(); ++q ) {
    assertMatrixNear( element->qps[q].managedStateVars->force,
                      expectedForce[q],
                      1e-10,
                      "2D interface force differs at qp " + std::to_string( q ) );
    assertMatrixNear( element->qps[q].managedStateVars->surfaceStress,
                      expectedSurfaceStress[q],
                      1e-10,
                      "2D interface surface stress differs at qp " + std::to_string( q ) );
  }
}

void TestAngledInterfaceKinematics()
{
  std::cout << "\n--- TestAngledInterfaceKinematics ---\n";

  constexpr int nDim   = 2;
  constexpr int nNodes = 4;

  const int  elId    = 1;
  const auto intType = FiniteElement::Quadrature::IntegrationTypes::FullIntegration;
  const auto secType = InterfaceFiniteElement< nDim, nNodes >::SectionType::Interface;

  auto element = std::make_unique< InterfaceFiniteElement< nDim, nNodes > >( elId, intType, secType );

  static std::array< double, nDim* nNodes > nodeCoordsVec = { 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0 };

  element->assignNodeCoordinates( nodeCoordsVec.data() );

  static std::array< double, 1 > elPropsVec = { 1.0 };
  ElementProperties              elProps( elPropsVec.data(), static_cast< int >( elPropsVec.size() ) );
  element->assignProperty( elProps );

  element->initializeYourself();

  for ( int q = 0; q < element->getNumberOfQuadraturePoints(); ++q ) {
    const auto& qp = element->qps[q];

    std::cout << "qp " << q << " xi = " << qp.xi.transpose() << " normal = " << qp.normal.transpose()
              << " norm = " << qp.normal.norm() << " sqrtDetG = " << qp.sqrtDetG << "\n";

    throwExceptionOnFailure( std::abs( qp.normal.norm() - 1.0 ) < 1e-12, "Interface normal is not normalized." );
    throwExceptionOnFailure( std::isfinite( qp.sqrtDetG ) && qp.sqrtDetG > 0.0,
                             "Angled interface sqrtDetG must be positive and finite." );
  }
}

void TestSharedMaterialSeesEachQuadraturePointCharacteristicLength()
{
  std::cout << "\n--- TestSharedMaterialSeesEachQuadraturePointCharacteristicLength ---\n";

  constexpr int nDim      = 3;
  constexpr int nNodes    = 8;
  constexpr int totalNDof = nDim * nNodes;

  // An IRREGULAR quad. A parallelogram has a constant surface Jacobian, so its quadrature points
  // all share one characteristic length and a shared material could not reveal a missing per-point
  // update. This trapezoid makes sqrtDetG vary across the points.
  auto element = std::make_unique<
    InterfaceFiniteElement< nDim, nNodes > >( 7,
                                              FiniteElement::Quadrature::IntegrationTypes::FullIntegration,
                                              InterfaceFiniteElement< nDim, nNodes >::SectionType::Interface );

  static std::array< double, nDim* nNodes > coordinates = {
    0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0,

    0.0, 0.0, 0.1, 3.0, 0.0, 0.1, 1.0, 1.0, 0.1, 0.0, 1.0, 0.1,
  };
  element->assignNodeCoordinates( coordinates.data() );

  static std::array< double, 1 > elPropsVec = { 1.0 };
  ElementProperties              elProps( elPropsVec.data(), static_cast< int >( elPropsVec.size() ) );
  element->assignProperty( elProps );

  static std::array< double, 3 > materialProperties = { 4000.0, 0.3, 0.01 };
  element->assignMaterial( "LINEARELASTIC", materialProperties.data(), materialProperties.size() );

  std::vector< double > stateVars( element->getNumberOfRequiredStateVars(), 0.0 );
  element->assignStateVars( stateVars.data(), static_cast< int >( stateVars.size() ) );
  element->initializeYourself();
  element->setInitialConditions( MarmotElement::MarmotMaterialInitialization, nullptr );

  std::vector< double > expectedLengths;
  for ( const auto& qp : element->qps )
    expectedLengths.push_back( qp.characteristicLength );

  const double spread = *std::max_element( expectedLengths.begin(), expectedLengths.end() ) -
                        *std::min_element( expectedLengths.begin(), expectedLengths.end() );

  std::cout << "characteristic-length spread across quadrature points = " << spread << "\n";

  throwExceptionOnFailure( spread > 1e-6,
                           "Test geometry must have quadrature points of differing characteristic length, "
                           "otherwise a shared material could not reveal a missing per-point update. spread = " +
                             std::to_string( spread ) );

  // swap in a material that records the length it is given at each evaluation; same base material and
  // therefore the same state-variable layout
  auto  recording   = std::make_unique< LengthRecordingInterfaceMaterial >();
  auto* recorder    = recording.get();
  element->material = std::move( recording );

  std::vector< double > U( totalNDof, 0.0 );
  std::vector< double > dU( totalNDof, 0.0 );
  std::vector< double > Pe( totalNDof, 0.0 );
  std::vector< double > Ke( totalNDof * totalNDof, 0.0 );

  for ( int a = 4; a < nNodes; ++a )
    dU[nDim * a + 0] = 1.0e-4;

  double time = 0.0;
  double dT   = 1.0;

  element->computeKernels( U.data(), dU.data(), Pe.data(), Ke.data(), time, dT );

  throwExceptionOnFailure( recorder->seenLengths.size() == expectedLengths.size(),
                           "Material was not evaluated once per quadrature point: " +
                             std::to_string( recorder->seenLengths.size() ) + " calls for " +
                             std::to_string( expectedLengths.size() ) + " points." );

  for ( size_t q = 0; q < expectedLengths.size(); ++q )
    throwExceptionOnFailure( std::abs( recorder->seenLengths[q] - expectedLengths[q] ) < 1e-12,
                             "Shared material saw the wrong characteristic length at quadrature point " +
                               std::to_string( q ) + ": got " + std::to_string( recorder->seenLengths[q] ) +
                               ", expected " + std::to_string( expectedLengths[q] ) + "." );
}

void TestElementShapeKeywordsAreEnSightNames()
{
  std::cout << "\n--- TestElementShapeKeywordsAreEnSightNames ---\n";

  // The result-file keyword is an EnSight cell name spanning the element's nodes, not the
  // computational interface shape: EnSight and ParaView reject "iquad4"/"iline2". The element
  // delegates to the geometry class, so the names live in one map there.
  {
    auto element = makeSingleInputFileInterfaceElement();
    throwExceptionOnFailure( element->getElementShape() == "hexa8",
                             "3D interface element must report the EnSight keyword \"hexa8\", got \"" +
                               element->getElementShape() + "\"." );
  }

  {
    constexpr int nDim   = 2;
    constexpr int nNodes = 4;

    auto element = std::make_unique<
      InterfaceFiniteElement< nDim, nNodes > >( 1,
                                                FiniteElement::Quadrature::IntegrationTypes::FullIntegration,
                                                InterfaceFiniteElement< nDim, nNodes >::SectionType::Interface );

    throwExceptionOnFailure( element->getElementShape() == "bar2",
                             "2D interface element must report the EnSight keyword \"bar2\", got \"" +
                               element->getElementShape() + "\"." );
  }
}

int main()
{
  auto tests = std::vector< std::function< void() > >{ TestMaterialInitializationResetsMaterialState,
                                                       TestUnsupportedInertiaThrows,
                                                       TestStressUpdateFailureRequestsSmallerTimeStep,
                                                       TestSingleInputFileElementGeometryMatrices,
                                                       TestSurfaceGradientOperatorProjectsOnlyGradientDirection,
                                                       TestSingleInputFileElementMaterialResponseIsFinite,
                                                       TestSingleInputFileElementGaussPointStiffnessAndResidual,
                                                       TestSingleInputFileElementRigidTranslationGivesZeroResidual,
                                                       TestAssignStateVarsPreservesHistoryAcrossIncrements,
                                                       TestTwoDimensionalInterfaceElementComputesWithEmbeddedMaterial,
                                                       TestAngledInterfaceKinematics,
                                                       TestSharedMaterialSeesEachQuadraturePointCharacteristicLength,
                                                       TestElementShapeKeywordsAreEnSightNames };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
