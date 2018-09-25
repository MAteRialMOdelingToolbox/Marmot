#pragma once
#include "userLibrary.h"
#include "bftUel.h"
#include "bftTypedefs.h"
#include <string>

namespace UelDisplacementFactory{

    BftUel* generateUelT2D2( 
                            int noEl);
    
    BftUel* generateUelCPS4( 
                            int noEl);

    BftUel* generateUelCPS8( 
                            int noEl);

    BftUel* generateUelCPE4( 
                            int noEl);

    BftUel* generateUelCPE8( 
                            int noEl);

    BftUel* generateUelCPS4R( 
                            int noEl);

    BftUel* generateUelCPS8R( 
                            int noEl);

    BftUel* generateUelCPE4R( 
                            int noEl);

    BftUel* generateUelCPE8R( 
                            int noEl);

    BftUel* generateUelC3D8( 
                            int noEl);

    BftUel* generateUelC3D8R( 
                            int noEl);

    BftUel* generateUelC3D20( 
                            int noEl);

    BftUel* generateUelC3D20R( 
                            int noEl);

}
