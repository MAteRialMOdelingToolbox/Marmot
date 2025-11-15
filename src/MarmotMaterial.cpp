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

void MarmotMaterial::initializeYourself( double* stateVars, int nStateVars )
{
  for ( auto i = 0; i < nStateVars; i++ )
    stateVars[i] = 0;
}

double MarmotMaterial::getDensity()
{
  return 1.0;
  throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << "not yet implementedXXX" );
}

MarmotMaterial::~MarmotMaterial() {}
