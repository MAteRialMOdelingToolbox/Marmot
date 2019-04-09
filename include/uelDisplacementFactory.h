#pragma once
#include "bftTypedefs.h"
#include "bftElement.h"
#include "userLibrary.h"
#include <string>

namespace UelDisplacementFactory {

    BftElement* generateUelT2D2( int noEl );
    BftElement* generateUelCPS4( int noEl );
    BftElement* generateUelCPS8( int noEl );
    BftElement* generateUelCPE4( int noEl );
    BftElement* generateUelCPE8( int noEl );
    BftElement* generateUelCPS4R( int noEl );
    BftElement* generateUelCPS8R( int noEl );
    BftElement* generateUelCPE4R( int noEl );
    BftElement* generateUelCPE8R( int noEl );
    BftElement* generateUelC3D8( int noEl );
    BftElement* generateUelC3D8R( int noEl );
    BftElement* generateUelC3D20( int noEl );
    BftElement* generateUelC3D20R( int noEl );

} // namespace UelDisplacementFactory
