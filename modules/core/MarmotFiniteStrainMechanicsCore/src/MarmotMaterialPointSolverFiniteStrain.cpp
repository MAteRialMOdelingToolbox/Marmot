#include "Marmot/MarmotMaterialPointSolverFiniteStrain.h"
#include "Marmot/MarmotMaterialFiniteStrainFactory.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>

using namespace Marmot;
using namespace Marmot::Solvers;
using namespace Fastor;
using namespace FastorStandardTensors;

MarmotMaterialPointSolverFiniteStrain::MarmotMaterialPointSolverFiniteStrain( std::string&         materialName,
                                                                              double*              materialProperties,
                                                                              int                  nMaterialProperties,
                                                                              const SolverOptions& options )
  : options( options )
{
  using namespace MarmotLibrary;

  // create material instance
  material = MarmotMaterialFiniteStrainFactory::createMaterial( materialName,
                                                                materialProperties,
                                                                nMaterialProperties,
                                                                1 );

  // get number of state variables
  nStateVars = material->getNumberOfRequiredStateVars();
  // initialize state variables
  stateVars = Eigen::VectorXd::Zero( nStateVars );
  // initialize material
  material->initializeYourself( stateVars.data(), nStateVars );
  // store initial state
  _initialStateVars = stateVars;
  stateVarsTemp     = stateVars;
}

void MarmotMaterialPointSolverFiniteStrain::addStep( const Step& step )
{
  step.checkControl();
  steps.push_back( step );
}

void MarmotMaterialPointSolverFiniteStrain::solve()
{
  for ( const auto& step : steps ) {
    std::cout << "Solving step from " << step.timeStart << " to " << step.timeEnd << std::endl;
    solveStep( step );
    std::cout << "+" + std::string( 78, '-' ) + "+" << std::endl;
  }
}

void MarmotMaterialPointSolverFiniteStrain::setInitialState( const Tensor33d&       initialStress,
                                                             const Eigen::VectorXd& initialStateVars )
{
  _initialStress    = initialStress;
  _initialStateVars = initialStateVars;
  stress            = _initialStress;
  stateVars         = _initialStateVars;
}

void MarmotMaterialPointSolverFiniteStrain::resetToInitialState()
{
  stress    = _initialStress;
  gradU     = Tensor33d( 0.0 );
  stateVars = _initialStateVars;
  stateVarsTemp.setZero();
  history.clear();
}

void MarmotMaterialPointSolverFiniteStrain::solveStep( const Step& step )
{
  double time     = step.timeStart;
  double dT       = step.dTStart;
  double stepTime = step.timeEnd - step.timeStart;

  int counter = 0;

  while ( time < step.timeEnd && counter <= step.maxIncrements ) {

    // adjust time step if overshooting
    if ( time + dT > step.timeEnd )
      dT = step.timeEnd - time;

    // setup increment
    Increment increment;
    increment.timeOld                     = time;
    increment.dT                          = dT;
    increment.gradUIncrement              = dT / stepTime * reshape< 9 >( step.gradUIncrementTarget );
    increment.stressIncrement             = dT / stepTime * reshape< 9 >( step.stressIncrementTarget );
    increment.isGradUComponentControlled  = reshape< 9 >( step.isGradUComponentControlled );
    increment.isStressComponentControlled = reshape< 9 >( step.isStressComponentControlled );

    // solve increment
    try {
      std::cout << "+" + std::string( 78, '-' ) + "+" << std::endl;
      std::cout << std::scientific << "  Solving increment " << counter + 1 << ", time: " << time << " to " << time + dT
                << ", dT: " << dT << std::endl;
      solveIncrement( increment );
      history.back().print();
      time += dT;
      stateVars = stateVarsTemp;
      counter++;
    }
    catch ( std::runtime_error& e ) {
      // if failed, reduce time step and retry
      if ( dT <= step.dTMin )
        throw std::runtime_error( "Minimum time step reached, cannot proceed." );
      dT = std::max( dT / 2.0, step.dTMin );
    }
  }

  if ( std::abs( time - step.timeEnd ) > 1e-12 )
    throw std::runtime_error( "Maximum number of increments reached, cannot proceed." );
}

void MarmotMaterialPointSolverFiniteStrain::solveIncrement( const Increment& increment )
{
  // set initial strain increment, set not controlled components to zero
  // use element-wise multiplication
  Tensor9d dGradU = increment.gradUIncrement * increment.isGradUComponentControlled.cast< double >();

  // create the mixed-control target
  Tensor9d target = reshape< 9 >( increment.stressIncrement );
  for ( int i = 0; i < 9; i++ )
    if ( increment.isGradUComponentControlled[i] )
      target[i] = increment.gradUIncrement[i];

  int    counter = 0;
  double resNorm = 1e12;
  double corNorm = 0.0;

  // temporary variables for material state
  MarmotMaterialFiniteStrain::ConstitutiveResponse< 3 > state;
  MarmotMaterialFiniteStrain::Deformation< 3 >          deformation;
  MarmotMaterialFiniteStrain::AlgorithmicModuli< 3 >    algorithmicModuli;
  MarmotMaterialFiniteStrain::TimeIncrement             timeInfo = { increment.timeOld + increment.dT, increment.dT };

  // store identical rows
  std::vector< std::pair< size_t, size_t > > identicalRows;

  // Newton-Raphson iteration
  while ( counter < options.maxIterations ) {
    std::cout << "    Iteration " << counter;

    // assign state variables to material
    stateVarsTemp = stateVars;

    // set state to previously converged state
    state.tau                  = stress;
    state.elasticEnergyDensity = 0.0;
    state.stateVars            = stateVarsTemp.data();

    // set deformation gradient
    deformation.F = Spatial3D::I + gradU + reshape< 3, 3 >( dGradU );

    // compute stress and tangent
    material->computeStress( state, algorithmicModuli, deformation, timeInfo );

    // initialize residual with stress increment
    Tensor9d residual = computeResidual( state.tau - stress, target, increment );

    // set tangent
    Tensor99d fullTangent = reshape< 9, 9 >( algorithmicModuli.dTau_dF );

    // modify tangent for mixed control
    modifyTangent( fullTangent, increment );

    // find identical rows which will be removed for inversion
    // only do this in the first iteration
    if ( counter == 0 ) {
      Tensor9d rand;
      rand.random();

      Tensor9d identicalRowCheck = fullTangent % rand;
      for ( size_t i = 0; i < 9; ++i ) {
        Tensor9d diff = abs( identicalRowCheck - identicalRowCheck[i] );
        for ( size_t j = i + 1; j < 9; ++j ) {
          if ( diff( j ) < 1e-14 ) {
            identicalRows.push_back( std::make_pair( i, j ) );
          }
        }
      }
    }

    // set the rows of the tangent matrix to zero for the identical rows and diagonal to 1
    // this ensures that the inversion will work correctly
    for ( const auto& rowPair : identicalRows ) {
      fullTangent( rowPair.second, all )            = 0.0;
      fullTangent( rowPair.second, rowPair.second ) = 1.0;
      residual( rowPair.second )                    = 0.0;
    }

    // compute residual norm
    resNorm = norm( residual );

    std::cout << std::scientific << ", ||ddGradU||: " << corNorm << ", ||R||: " << resNorm << std::endl;

    // if nan encountered
    if ( std::isnan( resNorm ) || std::isnan( corNorm ) )
      throw std::runtime_error( "NaN encountered in Newton-Raphson iteration." );

    // convergence check
    if ( corNorm < options.correctionTolerance && resNorm < options.residualTolerance )
      break;

    // if not converged,
    if ( counter == options.maxIterations - 1 )
      throw std::runtime_error( "Maximum number of iterations reached, no convergence." );

    // solve for correction
    Tensor9d correction = Fastor::solve( fullTangent, residual );

    // set the correction for the identical rows to the correction of the first row
    for ( const auto& rowPair : identicalRows ) {
      correction( rowPair.second ) = correction( rowPair.first );
    }

    // compute correction norm
    corNorm = norm( correction );

    // update displacement gradient increment
    dGradU -= correction;
    counter++;
  }

  std::cout << "    Converged after " << counter << " iterations." << std::endl;

  stress = state.tau;
  gradU += dGradU;
  dTau_dF = algorithmicModuli.dTau_dF;

  // copy back updated state variables
  stateVars = stateVarsTemp;
  history.push_back(
    HistoryEntry{ increment.timeOld + increment.dT, stress, Spatial3D::I + gradU, dTau_dF, stateVars } );
}

Tensor9d MarmotMaterialPointSolverFiniteStrain::computeResidual( const Tensor9d&  stressIncrement,
                                                                 const Tensor9d&  target,
                                                                 const Increment& increment )
{
  Tensor9d residual = stressIncrement;
  // replace displacement gradient controlled components
  for ( int i = 0; i < 9; i++ )
    if ( increment.isGradUComponentControlled[i] )
      residual[i] = increment.gradUIncrement[i];

  residual -= target;

  return residual;
}

void MarmotMaterialPointSolverFiniteStrain::modifyTangent( Tensor99d& tangent, const Increment& increment )
{
  // modify the tangent matrix based on control type
  for ( int i = 0; i < 9; i++ ) {
    if ( increment.isGradUComponentControlled[i] ) {
      tangent( i, all ) = 0.0;
      tangent( i, i )   = 1.0;
    }
  }
}

void MarmotMaterialPointSolverFiniteStrain::printHistory()
{
  std::cout << "Material Point History:" << std::endl;
  for ( const auto& entry : history ) {
    entry.print();
  }
}

void MarmotMaterialPointSolverFiniteStrain::exportHistoryToCSV( const std::string& filename )
{
  std::ofstream file( filename );
  if ( !file.is_open() )
    throw std::runtime_error( "Could not open file for writing: " + filename );

  // write header with fixed-width formatting
  const int w           = 15;
  const int ordering[9] = { 11, 12, 13, 21, 22, 23, 31, 32, 33 };
  file << std::scientific << "#" << std::setw( w - 1 ) << "Time,";

  // get ordering from array and write stress and strain components
  for ( int i = 0; i < 9; i++ )
    file << std::setw( w ) << "tau" + std::to_string( ordering[i] ) + ",";

  for ( int i = 0; i < 9; i++ )
    file << std::setw( w ) << "F" + std::to_string( ordering[i] ) + ",";

  for ( int i = 0; i < nStateVars; i++ )
    file << std::setw( w - ( i < nStateVars - 1 ? 1 : 2 ) ) << "SV" << i + 1 << ( i < nStateVars - 1 ? "," : "\n" );
  if ( nStateVars == 0 )
    file << "\n";

  // write data with fixed-width formatting
  for ( const auto& entry : history ) {
    file << std::scientific << std::setw( w - 1 ) << entry.time << ",";

    for ( int i = 0; i < 9; i++ )
      file << std::setw( w - 1 ) << flatten( entry.stress )[i] << ( i < 8 ? "," : "," );

    for ( int i = 0; i < 9; i++ )
      file << std::setw( w - 1 ) << flatten( entry.F )[i] << ( i < 8 ? "," : "," );

    for ( int i = 0; i < nStateVars; i++ )
      file << std::setw( w - 1 ) << entry.stateVars[i] << ( i < nStateVars - 1 ? "," : "\n" );
    if ( nStateVars == 0 )
      file << "\n";
  }

  file.close();
}
