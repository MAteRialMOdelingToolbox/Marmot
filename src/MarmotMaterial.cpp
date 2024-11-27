#include "Marmot/MarmotMaterial.h"
#include "Marmot/MarmotJournal.h"
#include <stdexcept>

MarmotMaterial::MarmotMaterial( const double* materialProperties_, int nMaterialProperties_, int materialNumber_ )
  : materialProperties( materialProperties_ ),
    nMaterialProperties( nMaterialProperties_ ),
    stateVars( nullptr ),
    nStateVars( 0 ),
    materialNumber( materialNumber_ )
{
}

void MarmotMaterial::assignStateVars( double* stateVars, int nStateVars )
{
  this->stateVars  = stateVars;
  this->nStateVars = nStateVars;
}

double* MarmotMaterial::getAssignedStateVars()
{
  return stateVars;
}

int MarmotMaterial::getNumberOfAssignedStateVars()
{
  return nStateVars;
}

void MarmotMaterial::initializeYourself()
{
  for ( auto i = 0; i < this->getNumberOfAssignedStateVars(); i++ )
    this->stateVars[i] = 0;
}

double MarmotMaterial::getDensity()
{

  throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << "not yet implemented" );
}

MarmotMaterial::~MarmotMaterial() {}
