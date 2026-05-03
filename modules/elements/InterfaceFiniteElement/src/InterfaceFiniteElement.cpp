#include "Marmot/InterfaceFiniteElement.h"

#include <cstddef>
#include <stdexcept>

namespace Marmot::Elements {

  template < int nDim, int nNodes >
  InterfaceFiniteElement< nDim, nNodes >::InterfaceFiniteElement(
    int                                         elementID,
    FiniteElement::Quadrature::IntegrationTypes integrationType,
    SectionType                                 sectionType )
    : ParentGeometryElement(),
      elementProperties( Eigen::Map< const Eigen::VectorXd >( nullptr, 0 ) ),
      elLabel( elementID ),
      sectionType( sectionType )
  {
    const auto qpInfos = FiniteElement::Quadrature::getGaussPointInfo( this->shape, integrationType );

    for ( const auto& qpInfo : qpInfos ) {
      QuadraturePoint qp( qpInfo.xi, qpInfo.weight );
      qps.push_back( std::move( qp ) );
    }
  }

  template < int nDim, int nNodes >
  int InterfaceFiniteElement< nDim, nNodes >::getNumberOfRequiredStateVars()
  {
    return qps[0].getNumberOfRequiredStateVars() * qps.size();
  }

  template < int nDim, int nNodes >
  std::vector< std::vector< std::string > > InterfaceFiniteElement< nDim, nNodes >::getNodeFields()
  {
    using namespace std;

    static vector< vector< string > > nodeFields;

    if ( nodeFields.empty() ) {
      for ( int i = 0; i < nNodes; i++ ) {
        nodeFields.push_back( vector< string >() );
        nodeFields[i].push_back( "displacement" );
      }
    }

    return nodeFields;
  }

  template < int nDim, int nNodes >
  std::vector< int > InterfaceFiniteElement< nDim, nNodes >::getDofIndicesPermutationPattern()
  {
    static std::vector< int > permutationPattern;

    if ( permutationPattern.empty() ) {
      for ( int i = 0; i < nNodes * nDim; i++ )
        permutationPattern.push_back( i );
    }

    return permutationPattern;
  }

  template < int nDim, int nNodes >
  void InterfaceFiniteElement< nDim, nNodes >::assignStateVars( double* stateVars, int nStateVars )
  {
    const int nQpStateVars = nStateVars / qps.size();

    for ( size_t i = 0; i < qps.size(); i++ ) {
      auto&   qp          = qps[i];
      double* qpStateVars = stateVars + ( i * nQpStateVars );
      qp.assignStateVars( qpStateVars, nQpStateVars );
    }
  }

  template < int nDim, int nNodes >
  void InterfaceFiniteElement< nDim, nNodes >::assignProperty( const ElementProperties& elementPropertiesInfo )
  {
    new ( &elementProperties ) Eigen::Map< const Eigen::VectorXd >( elementPropertiesInfo.elementProperties,
                                                                    elementPropertiesInfo.nElementProperties );
  }

  template < int nDim, int nNodes >
  void InterfaceFiniteElement< nDim, nNodes >::assignProperty( const MarmotMaterialSection& section )
  {
    for ( auto& qp : qps ) {
      qp.material = std::unique_ptr< Material >( dynamic_cast< Material* >(
        MarmotLibrary::MarmotMaterialHypoElasticInterfaceFactory::createMaterial( section.materialName,
                                                                                  section.materialProperties,
                                                                                  section.nMaterialProperties,
                                                                                  elLabel ) ) );

      if ( !qp.material ) {
        throw std::invalid_argument(
          MakeString() << __PRETTY_FUNCTION__
                       << ": invalid material assigned; cannot cast to MarmotMaterialHypoElasticInterface!" );
      }
    }
  }

  template < int nDim, int nNodes >
  void InterfaceFiniteElement< nDim, nNodes >::assignMaterial( const std::string& materialName,
                                                               const double*      materialProperties,
                                                               int                nMaterialProperties )
  {
    for ( auto& qp : qps ) {
      qp.material = std::unique_ptr< Material >( dynamic_cast< Material* >(
        MarmotLibrary::MarmotMaterialHypoElasticInterfaceFactory::createMaterial( materialName,
                                                                                  materialProperties,
                                                                                  nMaterialProperties,
                                                                                  elLabel ) ) );

      if ( !qp.material ) {
        throw std::invalid_argument(
          MakeString() << __PRETTY_FUNCTION__
                       << ": invalid material assigned; cannot cast to MarmotMaterialHypoElasticInterface!" );
      }
    }
  }

  template < int nDim, int nNodes >
  void InterfaceFiniteElement< nDim, nNodes >::assignNodeCoordinates( const double* coordinates )
  {
    ParentGeometryElement::assignNodeCoordinates( coordinates );
  }

  template < int nDim, int nNodes >
  void InterfaceFiniteElement< nDim, nNodes >::initializeYourself()
  {
    const double thickness = elementProperties.size() > 0 ? elementProperties[0] : 1.0;

    for ( QuadraturePoint& qp : qps ) {
      const bool fullyProjectedB = ( nDim == 3 );
      const auto geom            = this->evaluateAt( qp.xi, 0, fullyProjectedB );

      qp.N                 = geom.N;
      qp.dNdXi             = geom.dNdXi;
      qp.J                 = geom.J;
      qp.G                 = geom.G;
      qp.sqrtDetG          = geom.sqrtDetG;
      qp.detJ              = geom.sqrtDetG;
      qp.gradN             = geom.gradN;
      qp.normal            = geom.n;
      qp.normalProjection  = geom.normalProjection;
      qp.tangentProjection = geom.tangentProjection;

      qp.NmatSide    = geom.NmatSide;
      qp.BmatSide    = geom.BmatSide;
      qp.NmatJump    = geom.NmatJump;
      qp.BmatAverage = geom.BmatAverage;

      qp.J0xW = qp.weight * qp.sqrtDetG * thickness;

      if ( qp.material ) {
        if constexpr ( nDim == 3 ) {
          qp.material->setCharacteristicElementLength( std::sqrt( qp.sqrtDetG ) );
        }
        else if constexpr ( nDim == 2 ) {
          qp.material->setCharacteristicElementLength( qp.sqrtDetG );
        }
      }
    }
  }

  template < int nDim, int nNodes >
  void InterfaceFiniteElement< nDim, nNodes >::computeYourself( const double* QTotal_,
                                                                const double* dQ_,
                                                                double*       Pe_,
                                                                double*       Ke_,
                                                                const double* time,
                                                                double        dT,
                                                                double&       pNewDT )
  {
    (void) QTotal_;

    Eigen::Map< const RhsSized > dQ( dQ_ );
    Eigen::Map< KeSizedMatrix >  Ke( Ke_ );
    Eigen::Map< RhsSized >       Pe( Pe_ );

    constexpr int halfSize = nNodes * nDim / 2;

    for ( QuadraturePoint& qp : qps ) {
      const auto& Nside = qp.NmatSide;
      const auto& Bside = qp.BmatSide;
      const auto& Njump = qp.NmatJump;
      const auto& Bavg  = qp.BmatAverage;

      const auto dQBottom = dQ.template segment< halfSize >( 0 );
      const auto dQTop    = dQ.template segment< halfSize >( halfSize );

      InterfaceDisplSized dU_GPs;
      dU_GPs.template segment< nDim >( 0 )    = Nside * dQTop;
      dU_GPs.template segment< nDim >( nDim ) = Nside * dQBottom;

      InterfaceSurfaceGradSized dSurface_strain_GPs;
      dSurface_strain_GPs.template segment< nTensor >( 0 )       = Bside * dQTop;
      dSurface_strain_GPs.template segment< nTensor >( nTensor ) = Bside * dQBottom;

      ForceSized         force          = qp.managedStateVars->force;
      SurfaceStressSized surface_stress = qp.managedStateVars->surfaceStress;

      QMatrixSized Q_ij;
      ZMatrixSized Z_ijkl;
      HMatrixSized H_ijk;
      YMatrixSized Y_ijkl;

      Q_ij.setZero();
      Z_ijkl.setZero();
      H_ijk.setZero();
      Y_ijkl.setZero();

      qp.material->computeStress( force.data(),
                                  surface_stress.data(),
                                  Q_ij.data(),
                                  Z_ijkl.data(),
                                  H_ijk.data(),
                                  Y_ijkl.data(),
                                  dU_GPs.data(),
                                  dSurface_strain_GPs.data(),
                                  qp.normal.data(),
                                  time,
                                  dT,
                                  pNewDT );

      if ( pNewDT < 1.0 )
        return;

      qp.managedStateVars->force         = force;
      qp.managedStateVars->surfaceStress = surface_stress;
      qp.managedStateVars->displacement += dU_GPs;
      qp.managedStateVars->surfaceStrain += dSurface_strain_GPs;

      Pe -= Njump.transpose() * force * qp.J0xW;
      Pe -= Bavg.transpose() * surface_stress * qp.J0xW;

      Ke += ( Njump.transpose() * Q_ij * Njump
              + Bavg.transpose() * Z_ijkl * Bavg
              + Bavg.transpose() * Y_ijkl * Bavg
              + Njump.transpose() * H_ijk * Bavg
              + Bavg.transpose() * H_ijk.transpose() * Njump )
            * qp.J0xW;
    }
  }

  template < int nDim, int nNodes >
  void InterfaceFiniteElement< nDim, nNodes >::setInitialConditions( StateTypes state, const double* values )
  {
    switch ( state ) {
    case MarmotElement::MarmotMaterialInitialization: {
      /* Interface material state is managed through qp.managedStateVars in v26.05 style. */
      break;
    }

    case MarmotElement::MarmotMaterialStateVars: {
      throw std::invalid_argument( "Please use initializeStateVars directly on material" );
    }

    default:
      throw std::invalid_argument( MakeString()
                                   << __PRETTY_FUNCTION__ << ": invalid initial condition for InterfaceFiniteElement" );
    }
  }

  template < int nDim, int nNodes >
  void InterfaceFiniteElement< nDim, nNodes >::computeDistributedLoad( MarmotElement::DistributedLoadTypes loadType,
                                                                       double*                             P,
                                                                       double*                             K,
                                                                       const int                           elementFace,
                                                                       const double*                       load,
                                                                       const double*                       QTotal,
                                                                       const double*                       time,
                                                                       double                              dT )
  {
    throw std::invalid_argument(
      MakeString() << __PRETTY_FUNCTION__ << ": distributed loads are not implemented for InterfaceFiniteElement." );
  }

  template < int nDim, int nNodes >
  void InterfaceFiniteElement< nDim, nNodes >::computeBodyForce( double*       P,
                                                                 double*       K,
                                                                 const double* load,
                                                                 const double* QTotal,
                                                                 const double* time,
                                                                 double        dT )
  {
    throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__
                                              << ": body forces are not implemented for InterfaceFiniteElement." );
  }

  template < int nDim, int nNodes >
  void InterfaceFiniteElement< nDim, nNodes >::computeConsistentInertia( double* M )
  {
    Eigen::Map< KeSizedMatrix > Me( M );
    Me.setZero();
  }

  template < int nDim, int nNodes >
  void InterfaceFiniteElement< nDim, nNodes >::computeLumpedInertia( double* M )
  {
    Eigen::Map< RhsSized > Me( M );
    Me.setZero();
  }

  template < int nDim, int nNodes >
  std::vector< double > InterfaceFiniteElement< nDim, nNodes >::getCoordinatesAtCenter()
  {
    std::vector< double > coords( nDim );

    Eigen::Map< VectorDim > coordsMap( coords.data() );

    const auto centerXi = XiSized::Zero();
    const auto Ncenter  = this->N( centerXi );
    const auto Nmat     = this->NMatrix( Ncenter );

    const auto xSide = this->getSideCoordinates( 0 );

    coordsMap = Nmat * xSide;

    return coords;
  }

  template < int nDim, int nNodes >
  std::vector< std::vector< double > > InterfaceFiniteElement< nDim, nNodes >::getCoordinatesAtQuadraturePoints()
  {
    std::vector< std::vector< double > > listedCoords;

    for ( const auto& qp : qps ) {
      std::vector< double > coords( nDim );

      Eigen::Map< VectorDim > coordsMap( coords.data() );

      const auto xSide = this->getSideCoordinates( 0 );
      coordsMap        = qp.NmatSide * xSide;

      listedCoords.push_back( coords );
    }

    return listedCoords;
  }

  template < int nDim, int nNodes >
  int InterfaceFiniteElement< nDim, nNodes >::getNumberOfQuadraturePoints()
  {
    return static_cast< int >( qps.size() );
  }

  template class InterfaceFiniteElement< 2, 4 >;
  template class InterfaceFiniteElement< 3, 8 >;

} // namespace Marmot::Elements