#include "Marmot/MarmotMaterialPointSolverHypoElastic.h"
#include "Marmot/Marmot.h"
#include <fstream>

MarmotMaterialPointSolverHypoElastic::MarmotMaterialPointSolverHypoElastic( std::string&         materialName,
                                                                            double*              materialProperties,
                                                                            int                  nMaterialProperties,
                                                                            const SolverOptions& options )
  : options( options )
{
  using namespace MarmotLibrary;
  // get material code from name
  auto materialCode = MarmotMaterialFactory::getMaterialCodeFromName( materialName );

  // create material instance
  material = dynamic_cast< MarmotMaterialHypoElastic* >(
    MarmotMaterialFactory::createMaterial( materialCode, materialProperties, nMaterialProperties, 1 ) );
  // get number of state variables
  nStateVars = material->getNumberOfRequiredStateVars();
  // initialize state variables
  stateVars         = Eigen::VectorXd::Zero( nStateVars );
  _initialStateVars = Eigen::VectorXd::Zero( nStateVars );
  stateVarsTemp     = Eigen::VectorXd::Zero( nStateVars );

  material->assignStateVars( stateVarsTemp.data(), nStateVars );
}

void MarmotMaterialPointSolverHypoElastic::addStep( const Step& step )
{
  step.checkControl();
  steps.push_back( step );
}

void MarmotMaterialPointSolverHypoElastic::solve()
{
  for ( const auto& step : steps ) {
    std::cout << "Solving step from " << step.timeStart << " to " << step.timeEnd << std::endl;
    solveStep( step );
  }
}

void MarmotMaterialPointSolverHypoElastic::setInitialState( const Marmot::Vector6d& initialStress,
                                                            const Eigen::VectorXd&  initialStateVars )
{
  _initialStress    = initialStress;
  _initialStateVars = initialStateVars;
  stress            = _initialStress;
  stateVars         = _initialStateVars;
}

void MarmotMaterialPointSolverHypoElastic::resetToInitialState()
{
  stress = _initialStress;
  strain.setZero();
  stateVars = _initialStateVars;
  stateVarsTemp.setZero();
  history.clear();
}

void MarmotMaterialPointSolverHypoElastic::solveStep( const Step& step )
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
    increment.strainIncrement             = dT / stepTime * step.strainIncrementTarget;
    increment.stressIncrement             = dT / stepTime * step.stressIncrementTarget;
    increment.isStrainComponentControlled = step.isStrainComponentControlled;
    increment.isStressComponentControlled = step.isStressComponentControlled;

    // solve increment
    try {
      std::cout << "  Solving increment " << counter + 1 << ", time: " << time << " to " << time + dT << ", dT: " << dT
                << std::endl;
      solveIncrement( increment );
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

void MarmotMaterialPointSolverHypoElastic::solveIncrement( const Increment& increment )
{
  // set initial strain increment, set not controlled components to zero
  // use a Eigen Vector multiplication

  Marmot::Vector6d dStrain = increment.strainIncrement.array() *
                             increment.isStrainComponentControlled.array().cast< double >();

  // create the mixed-control target
  Marmot::Vector6d target = increment.stressIncrement;
  // replace strain controlled components
  // TODO: use Eigen vector operations
  for ( int i = 0; i < 6; i++ )
    if ( increment.isStrainComponentControlled[i] )
      target[i] = increment.strainIncrement[i];

  Marmot::Vector6d stressTemp = stress;
  Marmot::Matrix6d tangent, dStressDStrain;
  tangent.setZero();
  dStressDStrain.setZero();

  int    counter = 0;
  double pNewDT  = 1.0;
  double resNorm = 1e12;
  double corNorm = 0.0;

  std::cout << " dStrain initial: " << dStrain.transpose() << std::endl;

  // Newton-Raphson iteration
  while ( counter < options.maxIterations ) {
    std::cout << "    Iteration " << counter;

    // assign state variables to material
    stateVarsTemp = stateVars;

    // set stress to previous converged value
    stressTemp = stress;

    // set Abaqus style time array TODO: get rid of this
    double time[2] = { -1, increment.timeOld };

    std::cout << ", stress: " << stressTemp.transpose();

    // compute stress and tangent
    material->computeStress( stressTemp.data(), dStressDStrain.data(), dStrain.data(), &time[0], increment.dT, pNewDT );

    std::cout << ", stress: " << stressTemp.transpose();

    if ( pNewDT < 1.0 )
      throw std::runtime_error( "Material model requested time step reduction." );

    // initialize residual with stress increment
    Marmot::Vector6d residual = computeResidual( stressTemp - stress, target, increment );

    // set tangent
    tangent = dStressDStrain;

    // modify tangent for mixed control
    modifyTangent( tangent, increment );

    // compute residual norm
    resNorm = residual.norm();

    std::cout << std::scientific << ", ||ddE||: " << corNorm << ", ||R||: " << resNorm << std::endl;

    // convergence check
    if ( corNorm < options.correctionTolerance && resNorm < options.residualTolerance )
      break;

    // solve for correction
    Marmot::Vector6d correction = tangent.fullPivLu().solve( residual );
    corNorm                     = correction.norm();

    // update strain increment
    dStrain -= correction;
    counter++;
  }

  std::cout << "    Converged after " << counter << " iterations." << std::endl;

  stress = stressTemp;
  strain += dStrain;

  // copy back updated state variables
  stateVars = stateVarsTemp;
  history.push_back( HistoryEntry{ increment.timeOld + increment.dT, stress, strain, dStressDStrain, stateVars } );
}

Marmot::Vector6d MarmotMaterialPointSolverHypoElastic::computeResidual( const Marmot::Vector6d& stressIncrement,
                                                                        const Marmot::Vector6d& target,
                                                                        const Increment&        increment )
{
  Marmot::Vector6d residual = stressIncrement;
  // replace strain controlled components
  for ( int i = 0; i < 6; i++ )
    if ( increment.isStrainComponentControlled[i] )
      residual[i] = increment.strainIncrement[i];

  residual -= target;

  return residual;
}

void MarmotMaterialPointSolverHypoElastic::modifyTangent( Eigen::Matrix< double, 6, 6 >& tangent,
                                                          const Increment&               increment )
{
  // modify the tangent matrix based on control type
  for ( int i = 0; i < 6; i++ ) {
    if ( increment.isStrainComponentControlled[i] ) {
      tangent.row( i ).setZero();
      tangent( i, i ) = 1.0;
    }
  }
}

void MarmotMaterialPointSolverHypoElastic::printHistory()
{
  std::cout << "Material Point History:" << std::endl;
  for ( const auto& entry : history ) {
    entry.print();
  }
}

void MarmotMaterialPointSolverHypoElastic::exportHistoryToCSV( const std::string& filename )
{
  std::ofstream file( filename );
  if ( !file.is_open() )
    throw std::runtime_error( "Could not open file for writing: " + filename );

  // write header
  file << std::scientific << "# Time,";
  file << "Stress_11,Stress_22,Stress_33,Stress_12,Stress_13,Stress_23,";
  file << "Strain_11,Strain_22,Strain_33,Strain_12,Strain_13,Strain_23,";
  for ( int i = 0; i < nStateVars; i++ )
    file << "StateVar_" << i + 1 << ( i < nStateVars - 1 ? "," : "\n" );

  // write data
  for ( const auto& entry : history ) {
    file << entry.time << ",";
    for ( int i = 0; i < 6; i++ )
      file << entry.stress[i] << ( i < 5 ? "," : "," );
    for ( int i = 0; i < 6; i++ )
      file << entry.strain[i] << ( i < 5 ? "," : "," );
    for ( int i = 0; i < nStateVars; i++ )
      file << entry.stateVars[i] << ( i < nStateVars - 1 ? "," : "\n" );
  }

  file.close();
}
