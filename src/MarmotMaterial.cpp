#include "Marmot/MarmotMaterial.h"

MarmotMaterial::MarmotMaterial( const double*materialProperties_, int nMaterialProperties_, int materialNumber_ )
    : 
      materialProperties( materialProperties_ ),
      nMaterialProperties ( nMaterialProperties_ ),
      stateVars( nullptr ),
      nStateVars( 0 ),
      materialNumber ( materialNumber_ )
{}


void MarmotMaterial::assignStateVars( double* stateVars, int nStateVars )
{
    this->stateVars = stateVars;
    this->nStateVars = nStateVars;
}

const double* MarmotMaterial::getAssignedStateVars()
{
    return stateVars;
}
int MarmotMaterial::getNumberOfAssignedStateVars()
{
    return nStateVars;
}

MarmotMaterial::~MarmotMaterial(){}
