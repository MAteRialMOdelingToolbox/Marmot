#include "bftMaterial.h"

BftMaterial::BftMaterial( const double* materialProperties, int nMaterialProperties, int elementID, int gaussPt )
    : stateVars( nullptr ),
      nStateVars( 0 ),
      materialProperties( materialProperties ),
      nMaterialProperties( nMaterialProperties ),
      elementID( elementID ),
      gaussPt( gaussPt ){}

BftMaterial::~BftMaterial(){}
