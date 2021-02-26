#pragma once
#include "Eigen/Core"
#include "Fastor/Fastor.h"
#include "Marmot/MarmotTensor.h"

namespace Marmot {

  namespace FastorStandardTensors {

    using Tensor3d    = Fastor::Tensor< double, 3 >;
    using Tensor33d   = Fastor::Tensor< double, 3, 3 >;
    using Tensor333d  = Fastor::Tensor< double, 3, 3, 3 >;
    using Tensor3333d = Fastor::Tensor< double, 3, 3, 3, 3 >;

    using TensorMap3d    = Fastor::TensorMap< double, 3 >;
    using TensorMap33d   = Fastor::TensorMap< double, 3, 3 >;
    using TensorMap333d  = Fastor::TensorMap< double, 3, 3, 3 >;
    using TensorMap3333d = Fastor::TensorMap< double, 3, 3, 3, 3 >;

    namespace Spatial3D {

      inline const Tensor33d I = Tensor33d( ( Eigen::Matrix3d() << Eigen::Matrix3d::Identity() ).finished().data(),
                                            Fastor::ColumnMajor );

      inline const Tensor333d LeviCivita = Tensor333d( Marmot::ContinuumMechanics::CommonTensors::LeviCivita3D.data(),
                                                       Fastor::ColumnMajor );

      inline const Tensor3333d IHyd = Tensor3333d( Marmot::ContinuumMechanics::CommonTensors::I2xI2.data(),
                                                   Fastor::ColumnMajor );

      inline const Tensor3333d ISymm = Tensor3333d( Marmot::ContinuumMechanics::CommonTensors::Isym.data(),
                                                    Fastor::ColumnMajor );

      inline const Tensor3333d ISkew = Tensor3333d( Marmot::ContinuumMechanics::CommonTensors::Iskew.data(),
                                                    Fastor::ColumnMajor );

      inline const Tensor3333d I4 = Tensor3333d( Marmot::ContinuumMechanics::CommonTensors::IFourthOrder.data(),
                                                 Fastor::ColumnMajor );

      inline const Tensor3333d
        ITranspose = Tensor3333d( Marmot::ContinuumMechanics::CommonTensors::IFourthOrderTranspose.data(),
                                  Fastor::ColumnMajor );

      inline const Tensor3333d Deviatoric = I4 - 1. / 3 * IHyd;

      inline const Tensor3333d DeviatoricTranspose = Fastor::transpose( DeviatoricTranspose );

    } // namespace Spatial3D

  } // namespace FastorStandardTensors

  namespace ContinuumMechanics::Micropolar::GeneralizedInvariants {

    FastorStandardTensors::Tensor3333d d2J2_dStress_dStress( double a1, double a2 );

  }

  namespace FastorIndices {
    enum { i_, j_, k_, l_, m_, n_, A_, B_, I_, J_, K_, L_, M_, N_, P_ };

    using i    = Fastor::Index< i_ >;
    using ij   = Fastor::Index< i_, j_ >;
    using ik   = Fastor::Index< i_, k_ >;
    using ijk  = Fastor::Index< i_, j_, k_ >;
    using ijkl = Fastor::Index< i_, j_, k_, l_ >;
    using ijkB = Fastor::Index< i_, j_, k_, B_ >;
    using ijmn = Fastor::Index< i_, j_, m_, n_ >;
    using ijnk = Fastor::Index< i_, j_, n_, k_ >;
    using ijB  = Fastor::Index< i_, j_, B_ >;
    using in   = Fastor::Index< i_, n_ >;
    using inkB = Fastor::Index< i_, n_, k_, B_ >;
    using inB  = Fastor::Index< i_, n_, B_ >;
    using iA   = Fastor::Index< i_, A_ >;
    using iAkB = Fastor::Index< i_, A_, k_, B_ >;
    using iB   = Fastor::Index< i_, B_ >;
    using ji   = Fastor::Index< j_, i_ >;
    using jl   = Fastor::Index< j_, l_ >;
    using jin  = Fastor::Index< j_, i_, n_ >;
    using jkl  = Fastor::Index< j_, k_, l_ >;
    using jA   = Fastor::Index< j_, A_ >;
    using kl   = Fastor::Index< k_, l_ >;
    using kA   = Fastor::Index< k_, A_ >;
    using kB   = Fastor::Index< k_, B_ >;
    using lB   = Fastor::Index< l_, B_ >;
    using mnkB = Fastor::Index< m_, n_, k_, B_ >;
    using nB   = Fastor::Index< n_, B_ >;
    using A    = Fastor::Index< A_ >;
    using Ai   = Fastor::Index< A_, i_ >;
    using B    = Fastor::Index< B_ >;

    using iL   = Fastor::Index< i_, L_ >;
    using jk   = Fastor::Index< j_, k_ >;
    using iIKL = Fastor::Index< i_, I_, K_, L_ >;
    using ijKL = Fastor::Index< i_, j_, K_, L_ >;
    using iJKL = Fastor::Index< i_, J_, K_, L_ >;
    using iIjJ = Fastor::Index< i_, I_, j_, J_ >;
    using ij   = Fastor::Index< i_, j_ >;
    using jJ   = Fastor::Index< j_, J_ >;
    using IJ   = Fastor::Index< I_, J_ >;
    using IJKL = Fastor::Index< I_, J_, K_, L_ >;
    using IJML = Fastor::Index< I_, J_, M_, L_ >;
    using iI   = Fastor::Index< i_, I_ >;
    using iJ   = Fastor::Index< i_, J_ >;
    using iL   = Fastor::Index< i_, L_ >;
    using Jk   = Fastor::Index< J_, k_ >;
    using JI   = Fastor::Index< J_, I_ >;
    using JL   = Fastor::Index< J_, L_ >;
    using KI   = Fastor::Index< K_, I_ >;
    using KJ   = Fastor::Index< K_, J_ >;
    using KLMN = Fastor::Index< K_, L_, M_, N_ >;
    using KLm  = Fastor::Index< K_, L_, m_ >;
    using ijKJ = Fastor::Index< i_, j_, K_, J_ >;
    using Kk   = Fastor::Index< K_, k_ >;
    using IK   = Fastor::Index< I_, K_ >;
    using iIkK = Fastor::Index< i_, I_, k_, K_ >;
    using KL   = Fastor::Index< K_, L_ >;
    using MK   = Fastor::Index< M_, K_ >;

    using Mi   = Fastor::Index< M_, i_ >;
    using MJKL = Fastor::Index< M_, J_, K_, L_ >;
    using iK   = Fastor::Index< i_, K_ >;
    using iKjL = Fastor::Index< i_, K_, j_, L_ >;
    using im   = Fastor::Index< i_, m_ >;
    using jL   = Fastor::Index< j_, L_ >;
    using jLm  = Fastor::Index< j_, L_, m_ >;
    using KL   = Fastor::Index< K_, L_ >;
    using KLMP = Fastor::Index< K_, L_, M_, P_ >;
    using ML   = Fastor::Index< M_, L_ >;
    using Pm   = Fastor::Index< P_, m_ >;
    using Nm   = Fastor::Index< N_, m_ >;
    using MPm  = Fastor::Index< M_, P_, m_ >;
    using KLNM = Fastor::Index< K_, L_, N_, M_ >;
    using i    = Fastor::Index< i_ >;
    using j    = Fastor::Index< j_ >;
    using jk   = Fastor::Index< j_, k_ >;
    using jL   = Fastor::Index< j_, L_ >;
    using ij   = Fastor::Index< i_, j_ >;
    using ijk  = Fastor::Index< i_, j_, k_ >;
    using ijL  = Fastor::Index< i_, j_, L_ >;
    using ijLm = Fastor::Index< i_, j_, L_, m_ >;
    using ijm  = Fastor::Index< i_, j_, m_ >;
    using ik   = Fastor::Index< i_, k_ >;
    using iL   = Fastor::Index< i_, L_ >;
    using im   = Fastor::Index< i_, m_ >;
    using imk  = Fastor::Index< i_, m_, k_ >;
    using imL  = Fastor::Index< i_, m_, L_ >;
    using imLk = Fastor::Index< i_, m_, L_, k_ >;
    using k    = Fastor::Index< k_ >;
    using kl   = Fastor::Index< k_, l_ >;
    using km   = Fastor::Index< k_, m_ >;
    using kL   = Fastor::Index< k_, L_ >;
    using kM   = Fastor::Index< k_, M_ >;
    using l    = Fastor::Index< l_ >;
    using L    = Fastor::Index< L_ >;
    using Lm   = Fastor::Index< L_, m_ >;
    using LN   = Fastor::Index< L_, N_ >;
    using LmN  = Fastor::Index< L_, m_, N_ >;
    using m    = Fastor::Index< m_ >;
    using mj   = Fastor::Index< m_, j_ >;
    using mjL  = Fastor::Index< m_, j_, L_ >;

    using IJKL = Fastor::Index< I_, J_, K_, L_ >;
    using LI   = Fastor::Index< L_, I_ >;
    using JK   = Fastor::Index< J_, K_ >;
    using iK   = Fastor::Index< i_, K_ >;
    using iKjL = Fastor::Index< i_, K_, j_, L_ >;
    using im   = Fastor::Index< i_, m_ >;
    using jL   = Fastor::Index< j_, L_ >;
    using jLm  = Fastor::Index< j_, L_, m_ >;
    using KL   = Fastor::Index< K_, L_ >;
    using KLMP = Fastor::Index< K_, L_, M_, P_ >;
    using ML   = Fastor::Index< M_, L_ >;
    using LK   = Fastor::Index< L_, K_ >;
    using Pm   = Fastor::Index< P_, m_ >;
    using Nm   = Fastor::Index< N_, m_ >;
    using MPm  = Fastor::Index< M_, P_, m_ >;
    using KLPM = Fastor::Index< K_, L_, P_, M_ >;
    using KLMN = Fastor::Index< K_, L_, M_, N_ >;

    using to_ijm  = Fastor::OIndex< i_, j_, m_ >;
    using to_ijmM = Fastor::OIndex< i_, j_, m_, M_ >;

    using to_ijk  = Fastor::OIndex< i_, j_, k_ >;
    using to_ijLk = Fastor::OIndex< i_, j_, L_, k_ >;
    using to_ijL  = Fastor::OIndex< i_, j_, L_ >;
    using to_ijLm = Fastor::OIndex< i_, j_, L_, m_ >;

    using to_ijm  = Fastor::OIndex< i_, j_, m_ >;
    using to_ijmM = Fastor::OIndex< i_, j_, m_, M_ >;

    using to_Ii   = Fastor::OIndex< I_, i_ >;
    using to_IJKL = Fastor::OIndex< I_, J_, K_, L_ >;
    using to_iIjJ = Fastor::OIndex< i_, I_, j_, J_ >;
    using to_iIKL = Fastor::OIndex< i_, I_, K_, L_ >;
    using to_ijKL = Fastor::OIndex< i_, j_, K_, L_ >;
    using to_ijkL = Fastor::OIndex< i_, j_, k_, L_ >;
    using to_ijkK = Fastor::OIndex< i_, j_, k_, K_ >;
    using to_kK   = Fastor::OIndex< k_, K_ >;
    using to_kL   = Fastor::OIndex< k_, L_ >;

    using jkB = Fastor::Index< j_, k_, B_ >;
    using jB  = Fastor::Index< j_, B_ >;

    using to_jAkB = Fastor::OIndex< j_, A_, k_, B_ >;
    using to_jAB  = Fastor::OIndex< j_, A_, B_ >;
    using to_jkiB = Fastor::OIndex< j_, k_, i_, B_ >;
  } // namespace FastorIndices

  template < typename T,
             size_t nRows,
             size_t nCols,
             typename = std::enable_if< !std::is_const< std::remove_reference< T > >::value > >
  auto inline mapEigenToFastor( const Fastor::Tensor< T, nRows, nCols >& fastor )
  {
    return Eigen::Map< const Eigen::Matrix< T, nRows, nCols, Eigen::RowMajor > >( fastor.data() );
  }

  template < typename T, size_t nRows, size_t nCols, typename = void >
  auto inline mapEigenToFastor( Fastor::Tensor< T, nRows, nCols >& fastor )
  {
    return Eigen::Map< Eigen::Matrix< T, nRows, nCols, Eigen::RowMajor > >( fastor.data() );
  }

  template < template < typename, size_t... > class TensorType, typename T, size_t... Rest >
  void inline copyFastorToColumnMajor( T* target, const TensorType< T, Rest... >& source )
  {
    ( Fastor::TensorMap< T, Rest... >( target ) ) = Fastor::torowmajor( source );
  }

  enum DimensionType { U, W };
  template < DimensionType... dims, typename T, size_t... dims3D >
  auto inline reduceTo2D( const Fastor::Tensor< T, dims3D... >& theTensor3D )
  {

    static_assert( sizeof...( dims ) == sizeof...( dims3D ),
                   "CONVERSION PACK (1) AND TENSOR-ORDER PACK (2) MUST HAVE SAME LENGTH" );

    return theTensor3D( Fastor::fseq < dims == U ? 0 : 2, dims == U ? 2 : 3 > ()... );
  }

  template < typename Derived, size_t order >
  auto inline reduceTo2D( const Fastor::AbstractTensor< Derived, order >& theTensor3D )
  {
    using result_type = typename Derived::result_type;
    return reduceTo2D( result_type( theTensor3D ) );
  }

  constexpr int const3( size_t x ) { return 3; }

  template < typename T, size_t... dims2D >
  auto inline expandTo3D( const Fastor::Tensor< T, dims2D... >& theTensor2D )
  {
    static_assert( ( ( ( dims2D > 0 ) && ( dims2D < 3 ) ) && ... ),
                   "INPUT TENSOR IS NOT A VALID 2D MIXED DISPLACEMENT/ROTATION TENSOR" );

    Fastor::Tensor< double, const3( dims2D )... > theTensor3D( 0.0 );

    theTensor3D( Fastor::fseq < dims2D == 2 ? 0 : 2, dims2D == 2 ? 2 : 3 > ()... ) = theTensor2D;

    return theTensor3D;
  }

  template < typename Derived, size_t order >
  auto inline expandTo3D( const Fastor::AbstractTensor< Derived, order >& theTensor2D )
  {
    using result_type = typename Derived::result_type;
    return expandTo3D( result_type( theTensor2D ) );
  }

} // namespace Marmot
