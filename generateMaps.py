#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 23 19:55:26 2017

@author: matthias
"""

t = """
            
            //Displacement
            //Plane Stress
            UelCPS4 =       402,
            UelCPS8 =       802,
            UelCPS8R =      805,

            // Plane Strain
            UelCPE4 =       407,
            UelCPE8 =       807,
            UelCPE8R =       808,

            // Plane Strain - EAS
            UelCPE4EAS2 = 40702,
            UelCPE4EAS4 = 40704,
            UelCPE4EAS5 = 40705,
            
            // Solid
            UelC3D8 =       803,
            UelC3D8R =      806,
            UelC3D20 =      2003,
            UelC3D20R =     2006,

            // Solid EAS
            UelC3D8EAS3 =   80303,
            UelC3D8EAS9 =   80309,
                     
            // Nonlocal 
            // Plane Stress
            
            UelCPS4NonLocal = 412,
            UelCPS8NonLocal = 812,
            UelCPS8RNonLocal = 815,

            // Plane Stress - EAS
            UelCPS4NonLocalEAS2 =41202,
            UelCPS4NonLocalEAS4 =41204,

            // Plane Strain
            UelCPE4NonLocal =417,
            UelCPE4RNonLocal = 418,
            UelCPE8NonLocal =817,
            UelCPE8RNonLocal = 818,

            // Plane Strain - EAS
            UelCPE4NonLocalEAS2 =41702,
            UelCPE4NonLocalEAS4 =41704,
            UelC3D8NonLocalEAS3 = 81303,
            UelC3D8NonLocalEAS9 = 81309,

            // Solid
            UelC3D8NonLocal = 813,
            UelC3D8RNonLocal = 816,
            UelC3D20NonLocal = 2013,
            UelC3D20RNonLocal = 2016,

"""

lines = t.split('\n')

elNames = []

for l in lines:
    l = l.strip()
    if len(l) == 0:
        continue
    
    if '//' in l:
        continue
    
    Key, elCode = l.split('=')
    Key = Key.strip()
    elNames.append(Key)
    
for el in elNames:
    print('{{ "{:}", {:} }}, '.format(el, el))
    

