#include "Marmot/CompressibleHyperelasticity.h"
#include "Marmot/MarmotDeformationMeasures.h"
#include "Marmot/MarmotEnergyDensityFunctions.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include "Marmot/MarmotStressMeasures.h"
#include <Fastor/expressions/linalg_ops/unary_trans_op.h>
#include <Fastor/tensor/Tensor.h>
#include <Fastor/tensor_algebra/einsum_explicit.h>
#include <Fastor/tensor_algebra/indicial.h>
#include <cmath>
#include <stdexcept>
#include <string>
#include <variant>

namespace Marmot::Materials {

  using namespace Marmot;
  using namespace Fastor;
  using namespace Marmot::Math::Tensors::FastorIndices;
  using namespace Marmot::Math::Tensors::FastorTensors;

  namespace {
    namespace EDF = ContinuumMechanics::Elasticity::EnergyDensityFunctions;

    [[noreturn]] void throwInvalidProperties( const std::string& message )
    {
      throw std::invalid_argument( "CompressibleHyperelasticity: " + message );
    }

    /* A base selector is an integer stored as a double; anything else is rejected rather than truncated. */
    int baseSelector( const double value, const std::string& name )
    {
      if ( !( value >= 0. && value <= 1000. ) || std::abs( value - std::round( value ) ) > 1e-9 )
        throwInvalidProperties( "the " + name + " base must be a non-negative integer, got " + std::to_string( value ) +
                                "." );
      return static_cast< int >( std::round( value ) );
    }

    int nIsochoricParameters( const int isochoricBase )
    {
      switch ( isochoricBase ) {
      case CompressibleHyperelasticity::NeoHooke: return 1;
      case CompressibleHyperelasticity::Yeoh: return 3;
      case CompressibleHyperelasticity::MooneyRivlin: return 2;
      case CompressibleHyperelasticity::ArrudaBoyce8Chain: return 2;
      default:
        throwInvalidProperties( "unknown isochoric base " + std::to_string( isochoricBase ) +
                                " (0: NeoHooke, 1: Yeoh, 2: MooneyRivlin, 3: ArrudaBoyce8Chain)." );
      }
    }

    int nVolumetricParameters( const int volumetricBase )
    {
      switch ( volumetricBase ) {
      case CompressibleHyperelasticity::SussmanBathe:
      case CompressibleHyperelasticity::Penalty:
      case CompressibleHyperelasticity::SimoTaylor:
      case CompressibleHyperelasticity::LiuHofstetterMang:
      case CompressibleHyperelasticity::DollSchweizerhofExponential:
      case CompressibleHyperelasticity::DollSchweizerhofLogarithmic: return 1;
      case CompressibleHyperelasticity::Ogden: return 2;
      case CompressibleHyperelasticity::DollSchweizerhof: return 3;
      default:
        throwInvalidProperties( "unknown volumetric base " + std::to_string( volumetricBase ) +
                                " (1: SussmanBathe, 2: Penalty, 3: SimoTaylor, 4: Ogden, 5: LiuHofstetterMang, "
                                "6: DollSchweizerhof, 7: DollSchweizerhofExponential, "
                                "8: DollSchweizerhofLogarithmic)." );
      }
    }
  } // namespace

  CompressibleHyperelasticity::CompressibleHyperelasticity( const double* materialProperties,
                                                            int           nMaterialProperties,
                                                            int           materialLabel )
    : Marmot::Materials::MaterialFiniteStrain( materialProperties, nMaterialProperties, materialLabel )
  {
    if ( nMaterialProperties < 2 )
      throwInvalidProperties( "the first two material properties select the isochoric and the volumetric base." );

    const int isochoricBase  = baseSelector( materialProperties[0], "isochoric" );
    const int volumetricBase = baseSelector( materialProperties[1], "volumetric" );
    const int nIsochoric     = nIsochoricParameters( isochoricBase );
    const int nVolumetric    = nVolumetricParameters( volumetricBase );
    indexDensity             = 2 + nIsochoric + nVolumetric;

    if ( nMaterialProperties < indexDensity )
      throwInvalidProperties( "expected at least " + std::to_string( indexDensity ) +
                              " material properties (2 bases, " + std::to_string( nIsochoric ) + " isochoric and " +
                              std::to_string( nVolumetric ) + " volumetric parameters), got " +
                              std::to_string( nMaterialProperties ) + "." );

    if ( nMaterialProperties > indexDensity + 1 )
      throwInvalidProperties( "expected at most " + std::to_string( indexDensity + 1 ) +
                              " material properties (2 bases, " + std::to_string( nIsochoric ) + " isochoric and " +
                              std::to_string( nVolumetric ) + " volumetric parameters, optional density), got " +
                              std::to_string( nMaterialProperties ) + "." );

    const double* iso = materialProperties + 2;
    const double* vol = iso + nIsochoric;

    switch ( isochoricBase ) {
    case NeoHooke: isochoricPart = EDF::Isochoric::NeoHooke{ iso[0] }; break;
    case Yeoh: isochoricPart = EDF::Isochoric::Yeoh{ iso[0], iso[1], iso[2] }; break;
    case MooneyRivlin: isochoricPart = EDF::Isochoric::MooneyRivlin{ iso[0], iso[1] }; break;
    default:
      if ( !( iso[1] > 1.0 ) )
        throwInvalidProperties( "the Arruda-Boyce locking stretch lambdaL must be > 1, got " +
                                std::to_string( iso[1] ) + "." );
      isochoricPart = EDF::Isochoric::ArrudaBoyce8Chain{ iso[0], iso[1] };
    }

    switch ( volumetricBase ) {
    case SussmanBathe: volumetricPart = EDF::Volumetric::SussmanBathe{ vol[0] }; break;
    case Penalty: volumetricPart = EDF::Volumetric::Penalty{ vol[0] }; break;
    case SimoTaylor: volumetricPart = EDF::Volumetric::SimoTaylor{ vol[0] }; break;
    case Ogden:
      if ( vol[1] == 0.0 || !std::isfinite( vol[1] ) )
        throwInvalidProperties( "the Ogden exponent theta must be nonzero, got " + std::to_string( vol[1] ) + "." );
      volumetricPart = EDF::Volumetric::Ogden{ vol[0], vol[1] };
      break;
    case LiuHofstetterMang: volumetricPart = EDF::Volumetric::LiuHofstetterMang{ vol[0] }; break;
    case DollSchweizerhof:
      if ( !( vol[1] > 0.0 ) || !( vol[2] > 1.0 ) )
        throwInvalidProperties(
          "the Doll-Schweizerhof shape parameters must fulfil alpha > 0 and beta > 1, got alpha = " +
          std::to_string( vol[1] ) + ", beta = " + std::to_string( vol[2] ) + "." );
      volumetricPart = EDF::Volumetric::DollSchweizerhof{ vol[0], vol[1], vol[2] };
      break;
    case DollSchweizerhofExponential: volumetricPart = EDF::Volumetric::DollSchweizerhofExponential{ vol[0] }; break;
    default: volumetricPart = EDF::Volumetric::DollSchweizerhofLogarithmic{ vol[0] };
    }

    stateLayout.finalize();
  }

  double CompressibleHyperelasticity::getDensity( const double* stateVars ) const
  {
    if ( this->nMaterialProperties <= indexDensity )
      throw std::runtime_error( std::string( Marmot::Utilities::MakeString()
                                             << __PRETTY_FUNCTION__ << ": No density given! nMaterialProperties < "
                                             << indexDensity + 1 << "." ) );
    return this->materialProperties[indexDensity];
  }

  void CompressibleHyperelasticity::computeStress( ConstitutiveResponse< 3 >& response,
                                                   AlgorithmicModuli< 3 >&    tangents,
                                                   const Deformation< 3 >&    deformation,
                                                   const TimeIncrement&       timeIncrement ) const
  {
    const auto& F_ = deformation.F;

    using namespace ContinuumMechanics;
    // compute Cauchy-Green deformation
    const auto [C, dC_dF] = Marmot::ContinuumMechanics::Kinematics::DeformationMeasures::FirstOrder::rightCauchyGreen(
      F_ );

    // compute energy density, first and second partial derivatives wrt Cauchy Green deformation
    const auto [psi_, dPsi_dC, d2Psi_dCdC] = std::
      visit( [&C]( const auto& isochoric,
                   const auto&
                     volumetric ) { return EDF::SecondOrder::DecoupledPotential( C, isochoric, volumetric ); },
             isochoricPart,
             volumetricPart );

    // compute Kirchhoff stress
    Tensor33d PK2 = 2. * dPsi_dC;

    const auto [tau,
                dTau_dPK2,
                dTau_dF] = Marmot::ContinuumMechanics::StressMeasures::FirstOrder::KirchhoffStressFromPK2( PK2, F_ );
    response.tau         = tau;
    response.elasticEnergyDensity = psi_;
    response.dissipation          = 0.0;

    // compute tangent operator
    tangents.dTau_dF = 2.0 * einsum< ijKL, KLMN >( einsum< ijKL, IJKL >( dTau_dPK2, d2Psi_dCdC ), dC_dF ) + dTau_dF;
  }
} // namespace Marmot::Materials
