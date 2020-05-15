#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 15 11:47:51 2020

@author: magdalena
"""

import os

cwd = os.getcwd()

def cloneGitUibk(projectName,user='afbDevelopers'):
    os.system("git clone git@git.uibk.ac.at:{:}/{:}.git".format(user,projectName))

def installEigen():
    os.chdir(os.path.join(os.path.expanduser("~"),"Downloads"))
    os.system("git clone https://gitlab.com/libeigen/eigen.git")
    os.system("sudo cp -r eigen/unsupported /usr/include && sudo cp -r eigen/Eigen/ /usr/include")

def basicInstall(argument):    
    if argument == '1' or 'all':
         installEigen()        
        
    if argument in ['2','3','all']:
        os.chdir(os.path.join(cwd,'modules'))
        print("Installing to folder {:}".format(os.path.join(cwd,'modules')))
        if argument == '2' or 'all':    
            cloneGitUibk("bftFiniteElementCore","c8441141"),
        if argument == '3' or 'all':
            cloneGitUibk("bftCosseratCore","c8441141"),
    
def materialInstall(argument):
    os.chdir(os.path.join(cwd,'modules/materials'))
    print("Installing to folder {:}".format(os.path.join(cwd,'modules/materials')))
    if argument == '1' or 'all':   # basic materials
        cloneGitUibk("linearElastic") 
        cloneGitUibk("DruckerPrager","c8441141")
        cloneGitUibk("MohrCoulomb","c8441146")
    if argument == '2' or 'all':   # concrete materials
        cloneGitUibk("CosseratLinearElastic","c8441141")
        cloneGitUibk("ModLeon","c8441141")  
        cloneGitUibk("CDP","c8441141")  
        cloneGitUibk("CDPM2","c8441146")
    if argument == '3' or 'all':   # shotcrete materials
        cloneGitUibk("ModLeon","c8441141") 
        cloneGitUibk("ModLeonNonLocal","c8441141") 
        cloneGitUibk("ShotLeonV2","c8441141")  
        cloneGitUibk("SolidificationCreep","c8441141")  
    if argument == '4' or 'all':   # rock materials
        cloneGitUibk("HoekBrown","c8441146") 
        cloneGitUibk("RockDamagePlasticity","c8441146")  
        cloneGitUibk("RockDamagePlasticityNonLocal","c8441146") 
    if argument == '5' or 'all':   # rock materials
        cloneGitUibk("PorousElastic","c8441146") 
        cloneGitUibk("ModifiedCamClay","c8441146")  
        cloneGitUibk("Barodesy","c8441146") 
        cloneGitUibk("SandHypo","c8441146") 
    else:
        print("ERROR: No material for your input provided")
        
def elementInstall(argument):
    os.chdir(os.path.join(cwd,'modules/elements'))
    print("Installing to folder {:}".format(os.path.join(cwd,'modules/elements')))
    if argument == '1' or 'all':   # displacement elements
        cloneGitUibk("uelDisplacement", "c8441141") 
    if argument == '2' or 'all':   # large displacement
        cloneGitUibk("uelDisplacementUL", "c8441141") 
        cloneGitUibk("uelDisplacementTL", "c8441141") 
    if argument == '3' or 'all':   # gradient-enhanced 
        cloneGitUibk("uelNonLocal", "c8441141")
        cloneGitUibk("uelNonLocalMixed", "c8441141") 
        cloneGitUibk("uelNonLocalEAS", "c8441141")
    if argument == '4' or 'all':   # gradient-enhanced largedisplacement 
        cloneGitUibk("uelNonLocalUL", "c8441141") 
        cloneGitUibk("uelNonLocalULFBar", "c8441141") 
    if argument == '5' or 'all':   # Cosserat
        cloneGitUibk("uelCosserat", "c8441141") 
        cloneGitUibk("uelNonLocalCosserat", "c8441141") 
    else:
        print("ERROR: No element for your input provided")
    
if __name__=="__main__":
    
    print("""
+------------------------------------------------------------------------------
| Installation of the bftUserLibrary including the following components       |
+------------------------------------------------------------------------------

 bftUserlibrary/

     |-----modules/
             |---- bftMechanics/
             |---- materials/
             |---- elements/                (opt)
             |---- bftFiniteElementCore/    (opt)
             |---- bftCosseratCore/         (opt)

Additional external dependencies        (opt)

      |-----usr/include/
             |---- Eigen/        
             |---- unsupported/Eigen/ 
""")
    
    while True: 
        print("""
+-----------------------------------------------------------------------------+
| Basic Installation                                                          |
+-----------------------------------------------------------------------------+            
        """)
        txt = input("""Choose the following installation options:
            (1) Eigen library
            (2) basic bftFiniteElement library
            (3) basic bftCosserat library
            
Input e.g. 1, all, 1 2 3 (hit enter to continue): """)
        userList = txt.split()
        
        if not txt:
            break
        
        try:
            for i in userList: basicInstall(i)
        except:
            print(" no option chosen")
    
    while True: 
        print("""
+-----------------------------------------------------------------------------+
| Installation of materials                                                   |
+-----------------------------------------------------------------------------+       """)
        txt = input("""Choose the material libraries you would like to install:
            (1) basic materials
            (2) concrete materials
            (3) shotcrete materials
            (4) soil materials
            (5) rock materials
            
Input e.g. 1, all, 1 2 3 (hit enter to continue): """)
        userList = txt.split()
        
        if not txt:
            break
        
        try:
            for i in userList: materialInstall(i)
        except:
            print(" no option chosen")
    
    print("""
+-----------------------------------------------------------------------------+
| Installation of elements                                                     |
+-----------------------------------------------------------------------------+            """)
   
    while True: 
        txt = input("""Choose the element libraries you would like to install (e.g. 1, 1 2 3, or Enter to quit):
            (1) basic elements
            (2) large displacement elements (total/updated Lagrangian)
            (3) gradient-enhanced elements
            (4) gradient-enhanced large displacement elements
            (5) Cosserat elements
            
Input e.g. 1, all, 1 2 3 (hit enter to continue): """)
        userList = txt.split()
        
        if not txt:
            break
        
        try:
            for i in userList: elementInstall(i)
        except:
            print(" no option chosen")

