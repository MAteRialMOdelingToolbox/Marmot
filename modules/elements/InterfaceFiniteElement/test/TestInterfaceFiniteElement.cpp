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

  class FailingInterfaceMaterial : public MarmotInterfaceMaterialHypoElastic {
  public:
    FailingInterfaceMaterial() : MarmotInterfaceMaterialHypoElastic( nullptr, 0, 0 ) {}

    void computeStress( State&, Tangents&, const Deformation&, const TimeIncrement& ) override
    {
      throw Marmot::StressUpdateFailed( "Deliberate interface-material update failure." );
    }
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

    static std::array< double, 3 > materialProperties = { 4e5, 0.3, 0.01 };
    element->assignMaterial( "LINEARELASTICINTERFACE",
                             materialProperties.data(),
                             static_cast< int >( materialProperties.size() ) );

    return element;
  }

  void initializeStateAndMaterial( InterfaceFiniteElement< 3, 8 >& element, std::vector< double >& stateVars )
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

    static std::array< double, 8 > materialProperties = { 1e4, 0.3, 0.1, 1e-2, 1e-8, 1, 1e-2, 1.0 };
    element->assignMaterial( "LINEARVISCOELASTICINTERFACE", materialProperties.data(), materialProperties.size() );

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

    for ( auto& qp : element->qps ) {
      qp.material = std::make_unique< FailingInterfaceMaterial >();
    }

    constexpr int                                     nElementDofs = 3 * 8;
    std::array< double, nElementDofs >                QTotal{};
    std::array< double, nElementDofs >                dQ{};
    std::array< double, nElementDofs >                Pe{};
    std::array< double, nElementDofs * nElementDofs > Ke{};
    const std::array< double, 2 >                     time   = { 0.0, 0.0 };
    double                                            pNewDT = 1.0;

    element->computeYourself( QTotal.data(), dQ.data(), Pe.data(), Ke.data(), time.data(), 1.0, pNewDT );

    throwExceptionOnFailure( pNewDT == 0.5, "InterfaceFiniteElement did not request a smaller time step." );
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
                                                           double                          dT,
                                                           double&                         pNewDT )
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
      qp.material->computeStress( state, tangents, deformation, timeIncrement );

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

      const Eigen::VectorXd Pqp = -Njump.transpose() * force * J0xW - Bavg.transpose() * surfaceStress * J0xW;

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

void TestSingleInputFileElementMaterialResponseIsFinite()
{
  std::cout << "\n--- TestSingleInputFileElementMaterialResponseIsFinite ---\n";

  auto element = makeSingleInputFileInterfaceElement();

  std::vector< double > stateVars;
  initializeStateAndMaterial( *element, stateVars );

  const Eigen::VectorXd dU     = makeTopSideXOpeningIncrement();
  double                time   = 0.0;
  double                dT     = 0.1;
  double                pNewDT = 1.0;

  const auto expected = computeExpectedResponseFromGaussPoints( *element, dU, &time, dT, pNewDT );

  throwExceptionOnFailure( expected.Pe.allFinite(), "Material-only expected residual contains nan or inf." );
  throwExceptionOnFailure( expected.Pe.template lpNorm< Eigen::Infinity >() > 0.0,
                           "The top-side x opening should produce a nonzero material residual." );
}

void TestSingleInputFileElementLinearElasticGaussPointStiffnessAndResidual()
{
  std::cout << "\n--- TestSingleInputFileElementLinearElasticGaussPointStiffnessAndResidual ---\n";

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
  double pNewDTExpected  = 1.0;
  double pNewDTElement   = 1.0;

  const auto expected = computeExpectedResponseFromGaussPoints( *expectedElement,
                                                                dUEigen,
                                                                &timeForExpected,
                                                                dT,
                                                                pNewDTExpected );

  testedElement->computeYourself( U.data(), dU.data(), Pe.data(), Ke.data(), &timeForElement, dT, pNewDTElement );

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

    double timeLocal   = 0.0;
    double dTLocal     = 0.1;
    double pNewDTLocal = 1.0;

    element->computeYourself( Ulocal.data(),
                              dUlocal.data(),
                              Pelocal.data(),
                              Kelocal.data(),
                              &timeLocal,
                              dTLocal,
                              pNewDTLocal );

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

  const double err = ( KeActual + KeFiniteDifference ).template lpNorm< Eigen::Infinity >();

  const double fdNorm = KeFiniteDifference.template lpNorm< Eigen::Infinity >();

  const double relErr = err / std::max( 1.0, fdNorm );

  Eigen::Index maxRow         = 0;
  Eigen::Index maxCol         = 0;
  const double maxAbsMismatch = ( KeActual + KeFiniteDifference ).cwiseAbs().maxCoeff( &maxRow, &maxCol );

  std::cout << "max |Ke + KeFD| = " << maxAbsMismatch << " at (" << maxRow << ", " << maxCol << ")\n";
  std::cout << "KeActual(" << maxRow << "," << maxCol << ") = " << KeActual( maxRow, maxCol ) << "\n";
  std::cout << "KeFiniteDifference(" << maxRow << "," << maxCol << ") = " << KeFiniteDifference( maxRow, maxCol )
            << "\n";

  std::cout << "KeActual row " << maxRow << ": " << KeActual.row( maxRow ) << "\n";
  std::cout << "KeFD row " << maxRow << ": " << KeFiniteDifference.row( maxRow ) << "\n";
  std::cout << "KeActual + KeFD row " << maxRow << ": " << ( KeActual + KeFiniteDifference ).row( maxRow ) << "\n";

  std::cout << "stiffness finite-difference check: "
            << "relative error for Ke = -dPe/ddU is " << relErr << "\n";

  throwExceptionOnFailure( relErr < 1e-5, "Assembled stiffness is not consistent with Ke = -dPe/ddU." );
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

  double time   = 0.0;
  double dT     = 0.1;
  double pNewDT = 1.0;

  element->computeYourself( U.data(), dU.data(), Pe.data(), Ke.data(), &time, dT, pNewDT );

  Eigen::Map< Eigen::Matrix< double, totalNDof, 1 > > PeActual( Pe.data() );

  throwExceptionOnFailure( PeActual.allFinite(), "Rigid-translation residual contains nan or inf." );
  throwExceptionOnFailure( PeActual.template lpNorm< Eigen::Infinity >() < 1e-9,
                           "Equal rigid translation on both sides should produce zero interface residual." );
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

int main()
{
  auto tests = std::vector<
    std::function< void() > >{ TestMaterialInitializationResetsMaterialState,
                               TestUnsupportedInertiaThrows,
                               TestStressUpdateFailureRequestsSmallerTimeStep,
                               TestSingleInputFileElementGeometryMatrices,
                               TestSingleInputFileElementMaterialResponseIsFinite,
                               TestSingleInputFileElementLinearElasticGaussPointStiffnessAndResidual,
                               TestSingleInputFileElementRigidTranslationGivesZeroResidual,
                               TestAngledInterfaceKinematics };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
