#pragma once

#include "Marmot/Marmot.h"
#include "Marmot/MarmotConstants.h"
#include "Marmot/MarmotElement.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotGeometryInterfaceElement.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMaterialHypoElasticInterface.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotStateVarVectorManager.h"
#include "Marmot/MarmotTypedefs.h"

#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

using namespace Marmot;
using namespace Eigen;

namespace Marmot::Elements {

  template < int nDim, int nNodes >
  class InterfaceFiniteElement
    : public MarmotElement,
      public MarmotGeometryInterfaceElement< nDim, nNodes > {

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    enum SectionType {
      Interface,
    };

    static constexpr int nDofPerNodeU    = nDim;
    static constexpr int nInterfaceNodes = nNodes / 2;

    static constexpr int sizeLoadVector = nNodes * nDim;
    static constexpr int nCoordinates   = nNodes * nDim;

    static constexpr int nSideDofs = nInterfaceNodes * nDim;
    static constexpr int nTensor   = nDim * nDim;

    using ParentGeometryElement = MarmotGeometryInterfaceElement< nDim, nNodes >;

    using XiSized              = typename ParentGeometryElement::XiSized;
    using NSized               = typename ParentGeometryElement::NSized;
    using dNdXiSized           = typename ParentGeometryElement::dNdXiSized;
    using SurfaceJacobianSized = typename ParentGeometryElement::SurfaceJacobianSized;
    using MetricSized          = typename ParentGeometryElement::MetricSized;
    using GradSized            = typename ParentGeometryElement::GradSized;

    using VectorDim = typename ParentGeometryElement::VectorDim;
    using TensorDim = typename ParentGeometryElement::TensorDim;

    using NMatrixSized     = typename ParentGeometryElement::NMatrixSized;
    using NJumpMatrixSized = typename ParentGeometryElement::NJumpMatrixSized;

    using BSurfaceSized    = typename ParentGeometryElement::BSurfaceSized;
    using BAvgSurfaceSized = typename ParentGeometryElement::BAvgSurfaceSized;

    using RhsSized      = Matrix< double, sizeLoadVector, 1 >;
    using KeSizedMatrix = Matrix< double, sizeLoadVector, sizeLoadVector >;

    using ForceSized                = Matrix< double, nDim, 1 >;
    using SurfaceStressSized        = Matrix< double, nTensor, 1 >;
    using InterfaceDisplSized       = Matrix< double, 2 * nDim, 1 >;
    using InterfaceSurfaceGradSized = Matrix< double, 2 * nTensor, 1 >;

    /*
     * Important:
     * These material tangent matrices are passed to a C/C++ material function
     * through raw double* pointers, while the Python implementation expects
     * C-order / row-major tensor flattening.
     *
     * Eigen is column-major by default, so these must be RowMajor to avoid
     * transposed/scrambled Q, Z, H, Y when the material writes into data().
     */
    using QMatrixSized = Matrix< double, nDim, nDim, RowMajor >;
    using ZMatrixSized = Matrix< double, nTensor, nTensor, RowMajor >;
    using HMatrixSized = Matrix< double, nDim, nTensor, RowMajor >;
    using YMatrixSized = Matrix< double, nTensor, nTensor, RowMajor >;

    using Material = MarmotMaterialHypoElasticInterface;

    Map< const VectorXd > elementProperties;

    const int         elLabel;
    const SectionType sectionType;

    struct QuadraturePoint {
      EIGEN_MAKE_ALIGNED_OPERATOR_NEW

      const XiSized xi;
      const double  weight;

      double detJ;
      double sqrtDetG;
      double J0xW;

      NSized               N;
      dNdXiSized           dNdXi;
      SurfaceJacobianSized J;
      MetricSized          G;
      GradSized            gradN;

      VectorDim normal;
      TensorDim normalProjection;
      TensorDim tangentProjection;

      /*
       * One-side operators:
       *
       * NmatSide:
       *   maps one side's nodal displacement vector to u at the interface qp.
       *
       * BmatSide:
       *   maps one side's nodal displacement vector to the surface-gradient
       *   quantity passed to the material.
       */
      NMatrixSized  NmatSide;
      BSurfaceSized BmatSide;

      /*
       * Whole-element operators:
       *
       * NmatJump:
       *   jump u = u_top - u_bottom
       *   NmatJump = [ -Nside , +Nside ]
       *
       * BmatAverage:
       *   grad_s u_avg = 0.5 * (grad_s u_bottom + grad_s u_top)
       *   BmatAverage = 0.5 * [ Bside , Bside ]
       */
      NJumpMatrixSized NmatJump;
      BAvgSurfaceSized BmatAverage;

      class QPStateVarManager : public MarmotStateVarVectorManager {

        /*
         * Python-compatible persistent state layout:
         *
         *   force
         *   surface stress
         *   displacement
         *   surface strain
         *   material state variables
         *
         * The displacement and surface strain entries store the accumulated
         * top/bottom quantities:
         *
         *   displacement   = [u_top, u_bottom]
         *   surface strain = [grad_s u_top, grad_s u_bottom]
         */
        inline const static auto layout = makeLayout( {
          { .name = "force", .length = nDim },
          { .name = "surface stress", .length = nDim * nDim },
          { .name = "displacement", .length = 2 * nDim },
          { .name = "surface strain", .length = 2 * nDim * nDim },
          { .name = "begin of material state", .length = 0 },
        } );

      public:
        Eigen::Map< ForceSized >                force;
        Eigen::Map< SurfaceStressSized >        surfaceStress;
        Eigen::Map< InterfaceDisplSized >       displacement;
        Eigen::Map< InterfaceSurfaceGradSized > surfaceStrain;

        Eigen::Map< Eigen::VectorXd > materialStateVars;

        static int getNumberOfRequiredStateVarsQuadraturePointOnly()
        {
          return layout.nRequiredStateVars;
        }

        QPStateVarManager( double* theStateVarVector, int nStateVars )
          : MarmotStateVarVectorManager( theStateVarVector, layout ),
            force( &find( "force" ) ),
            surfaceStress( &find( "surface stress" ) ),
            displacement( &find( "displacement" ) ),
            surfaceStrain( &find( "surface strain" ) ),
            materialStateVars( &find( "begin of material state" ),
                               nStateVars - getNumberOfRequiredStateVarsQuadraturePointOnly() )
        {
        }
      };

      std::unique_ptr< QPStateVarManager > managedStateVars;
      std::unique_ptr< Material >          material;

      int getNumberOfRequiredStateVarsQuadraturePointOnly()
      {
        return QPStateVarManager::getNumberOfRequiredStateVarsQuadraturePointOnly();
      }

      int getNumberOfRequiredStateVars()
      {
        return getNumberOfRequiredStateVarsQuadraturePointOnly()
             + material->getNumberOfRequiredStateVars();
      }

      void assignStateVars( double* stateVars, int nStateVars )
      {
        managedStateVars = std::make_unique< QPStateVarManager >( stateVars, nStateVars );

        material->assignStateVars( managedStateVars->materialStateVars.data(),
                                   managedStateVars->materialStateVars.size() );
      }

      QuadraturePoint( XiSized xi, double weight )
        : xi( xi ),
          weight( weight ),
          detJ( 0.0 ),
          sqrtDetG( 0.0 ),
          J0xW( 0.0 ),
          N( NSized::Zero() ),
          dNdXi( dNdXiSized::Zero() ),
          J( SurfaceJacobianSized::Zero() ),
          G( MetricSized::Zero() ),
          gradN( GradSized::Zero() ),
          normal( VectorDim::Zero() ),
          normalProjection( TensorDim::Zero() ),
          tangentProjection( TensorDim::Zero() ),
          NmatSide( NMatrixSized::Zero() ),
          BmatSide( BSurfaceSized::Zero() ),
          NmatJump( NJumpMatrixSized::Zero() ),
          BmatAverage( BAvgSurfaceSized::Zero() )
      {
      }
    };

    std::vector< QuadraturePoint, Eigen::aligned_allocator< QuadraturePoint > > qps;

    InterfaceFiniteElement( int                                         elementID,
                            FiniteElement::Quadrature::IntegrationTypes integrationType,
                            SectionType                                 sectionType = SectionType::Interface );

    int getNumberOfRequiredStateVars();

    std::vector< std::vector< std::string > > getNodeFields();

    std::vector< int > getDofIndicesPermutationPattern();

    int getNNodes() { return nNodes; }

    int getNSpatialDimensions() { return nDim; }

    int getNDofPerElement() { return sizeLoadVector; }

    /*
     * Important:
     * ParentGeometryElement::getElementShape() returns the computational
     * interface shape, e.g. "iquad4".
     *
     * EnSight/ParaView do not understand "iquad4" as a geometry keyword.
     * Therefore, for output/visualization we map interface elements to valid
     * EnSight element names.
     *
     * This mirrors the Python EdelweissFE provider:
     *
     *   IQuad4 -> ensightType = "hexa8"
     *
     * The computational element remains an interface element. This only affects
     * the geometry keyword written to result files.
     */
    std::string getElementShape()
    {
      if constexpr ( nDim == 3 && nNodes == 8 ) {
        return "hexa8";
      }
      else if constexpr ( nDim == 2 && nNodes == 4 ) {
        return "line2";
      }
      else {
        return ParentGeometryElement::getElementShape();
      }
    }

    void assignStateVars( double* stateVars, int nStateVars );

    void assignProperty( const ElementProperties& marmotElementProperty );

    void assignProperty( const MarmotMaterialSection& marmotElementProperty );

    void assignMaterial( const std::string& materialName,
                         const double*      materialProperties,
                         int                nMaterialProperties ) override;

    void assignNodeCoordinates( const double* coordinates );

    void initializeYourself();

    void setInitialConditions( StateTypes state, const double* values );

    void computeDistributedLoad( MarmotElement::DistributedLoadTypes loadType,
                                 double*                             P,
                                 double*                             K,
                                 const int                           elementFace,
                                 const double*                       load,
                                 const double*                       QTotal,
                                 const double*                       time,
                                 double                              dT );

    void computeBodyForce( double*       P,
                           double*       K,
                           const double* load,
                           const double* QTotal,
                           const double* time,
                           double        dT );

    void computeYourself( const double* QTotal,
                          const double* dQ,
                          double*       Pe,
                          double*       Ke,
                          const double* time,
                          double        dT,
                          double&       pNewdT );

    void computeConsistentInertia( double* M );

    void computeLumpedInertia( double* M );

    StateView getStateView( const std::string& stateName, int qpNumber )
    {
      const auto& qp = qps[qpNumber];

      if ( qp.managedStateVars->contains( stateName ) ) {
        return qp.managedStateVars->getStateView( stateName );
      }

      if ( stateName == "sdv" ) {
        std::cout << __PRETTY_FUNCTION__
                  << " on 'sdv' is discouraged and deprecated, please use precise state name";
        return { qp.managedStateVars->materialStateVars.data(),
                 static_cast< int >( qp.managedStateVars->materialStateVars.size() ) };
      }

      return qp.material->getStateView( stateName );
    }

    std::vector< double > getCoordinatesAtCenter();

    std::vector< std::vector< double > > getCoordinatesAtQuadraturePoints();

    int getNumberOfQuadraturePoints();
  };

} // namespace Marmot::Elements