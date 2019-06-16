#pragma once
#include "bftTypedefs.h"
#include "bftElement.h"
#include "userLibrary.h"
#include <string>

namespace UelDisplacementFactory {

    BftElement* generateUelT2D2( int elementID );
    BftElement* generateUelCPS4( int elementID );
    BftElement* generateUelCPS8( int elementID );
    BftElement* generateUelCPE4( int elementID );
    BftElement* generateUelCPE8( int elementID );
    BftElement* generateUelCPS4R( int elementID );
    BftElement* generateUelCPS8R( int elementID );
    BftElement* generateUelCPE4R( int elementID );
    BftElement* generateUelCPE8R( int elementID );
    BftElement* generateUelC3D8( int elementID );
    BftElement* generateUelC3D8R( int elementID );
    BftElement* generateUelC3D20( int elementID );
    BftElement* generateUelC3D20R( int elementID );

} // namespace UelDisplacementFactory
