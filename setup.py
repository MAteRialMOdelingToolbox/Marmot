#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 15 11:47:51 2020

@author: magdalena
"""

import os
import shutil
from distutils.util import strtobool

cwd = os.getcwd()
defaultUser = 'afbDevelopers'


dictMaterials = { # add material: 
                  # name: ("basic/concrete/rock/shotcrete", "afbDevelopers/c8441141/c8441146")
                "LinearElastic":                     ("basic", "afbDevelopers"),
                "DruckerPrager":                     ("basic", "c8441141"),
                "MohrCoulomb":                       ("basic", "c8441146"),
                "CosseratLinearElastic":             ("cosserat", "c8441141"),
                "ModLeon":                           ("concrete", "c8441141"),
                "CDP":                               ("concrete", "c8441141"),
                "CDPM2":                             ("concrete", "c8441141"),
                "ModLeonNonLocal":                   ("concrete", "c8441141"),      
                "ShotLeonV2":                        ("concrete", "c8441141"),        
                "SolidificationCreep":               ("concrete", "c8441141"), 
                "HoekBrown":                         ("rock", "c8441146"),        
                "RockDamagePlasticity":              ("rock", "c8441146"),        
                "RockDamagePlasticityNonLocal":      ("rock", "c8441146"), 
                "ORDP":                              ("rock", "c8441192"),
                "JointedHoekBrown":                  ("rock", "c8441192"),
                "PorousElastic":                     ("soil", "c8441146"),        
                "ModifiedCamClay":                   ("soil", "c8441146"),        
                "Barodesy":                          ("soil", "c8441146"), 
                "SandHypo":                          ("soil", "c8441146"), 
                "SandHypoMicropolar":                ("soil", "c8441146"), 
                "BarodesGradientVoid":               ("soil", "c8441146"), 
                "BarodesGradientDeformationModulus": ("soil", "c8441146"), 
                }


dictElements = { # add element: 
                  # name: ("displacement/large displacement/gradient-enhanced/ 
                  #          gradient-enhanced large displacement/cosserat", 
                # "afbDevelopers/c8441141/c8441146")
                "DisplacementFiniteElement":       ("displacement",       "c8441141"),
                "DisplacementULFiniteElement":     ("large displacement", "c8441141"),
                "DisplacementTLFiniteElement":     ("large displacement", "c8441141"),
                "GradientEnhancedDisplacementFiniteElement":           ("gradient-enhanced",  "c8441141"),
                "GradientEnhancedDisplacementMixedFiniteElement":      ("gradient-enhanced",  "c8441146"),
                "GradientEnhancedDisplacementEASFiniteElement":        ("gradient-enhanced",  "c8441141"),
                "GradientEnhancedDisplacementULFiniteElement":         ("gradient-enhanced large displacement", "c8441141"),      
                "GradientEnhancedDisplacementULFiniteElement":         ("gradient-enhanced large displacement", "c8441141"),        
                "GradientEnhancedDisplacementULFBarFiniteElement":     ("gradient-enhanced large displacement", "c8441141"), 
                "CosseratFiniteElement":           ("cosserat", "c8441141"),        
                "GradientEnhancedCosseratFiniteElement":   ("cosserat", "c8441141"),      
                }

lookUpMaterialCategories = {  "basic": "1",
                              "concrete": "2",
                              "shotcrete": "3",
                              "soil": "4",
                              "rock": "5",
                              }

lookUpElementCategories = {  "displacement": "1",
                              "large displacement": "2",
                              "gradient-enhanced": "3",
                              "gradient-enhanced large displacement": "4",
                              "cosserat": "5",
                              }




def cloneGitUibk(projectName,user='afbDevelopers'):
    os.system("git clone git@git.uibk.ac.at:{:}/{:}.git".format(user,projectName))

def installEigen():
    os.chdir(os.path.join(os.path.expanduser("~"),"Downloads"))
    os.system("git clone https://gitlab.com/libeigen/eigen.git")
    os.system("sudo cp -r eigen/unsupported /usr/include && sudo cp -r eigen/Eigen/ /usr/include")

def basicInstall(argument):    
    if argument == '1' or 'all':
         installEigen()        
        
    if argument in ['2','3','4','all']:
        os.chdir(os.path.join(cwd,'modules'))
        print("Installing to folder {:}".format(os.path.join(cwd,'modules')))
        if argument == '2' or 'all':
            cloneGitUibk("MarmotMechanics")
        if argument == '3' or 'all':    
            cloneGitUibk("MarmotFiniteElementCore","c8441141")
        if argument == '4' or 'all':
            cloneGitUibk("MarmotCosseratCore","c8441141"),
    
def materialInstall(argument):
    os.chdir(os.path.join(cwd,'modules/materials'))
    print("Installing to folder {:}".format(os.path.join(cwd,'modules/materials')))

    for material, props in dictMaterials.items():
        category = lookUpMaterialCategories.get(props[0])
        if argument == 'all':
            cloneGitUibk(material, props[1])
        elif category == argument:
            cloneGitUibk(material, props[1])
        elif material.lower() == argument.lower():
            cloneGitUibk(material, props[1])
        
def elementInstall(argument):
    os.chdir(os.path.join(cwd,'modules/elements'))
    print("Installing to folder {:}".format(os.path.join(cwd,'modules/elements')))
    for element, props in dictElements.items():
        category = lookUpMaterialCategories.get(props[0])
        if argument == 'all':
            cloneGitUibk(element, props[1])
        elif category == argument:
            cloneGitUibk(element, props[1])
        elif element.lower() == argument.lower():
            cloneGitUibk(element, props[1])
    
if __name__=="__main__":
    
    print("""
+------------------------------------------------------------------------------
| Installation of the bftUserLibrary including the following components       |
+------------------------------------------------------------------------------

 bftUserlibrary/

     |-----modules/
             |---- bftMechanics/            (opt)
             |---- materials/               (opt)
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
            (2) bftMechanics
            (3) bftFiniteElementCore
            (4) bftCosseratCore
            
Input e.g. 1, all, 1 2 3 (hit enter to continue): """)
        userList = txt.split()
        
        if not txt:
            break
        
        # try:
        for i in userList: basicInstall(i)
        # except:
            # print(" no option chosen")
    
    while True: 
        print("""
+-----------------------------------------------------------------------------+
| Installation of materials                                                   |
+-----------------------------------------------------------------------------+       """)
        print("""Choose the material libraries you would like to install:
        
        either as complete packages
            (1) basic materials
            (2) concrete materials
            (3) shotcrete materials
            (4) soil materials
            (5) rock materials
        
        or simply type the name of your desired material from the following list: """)
        
        for i in dictMaterials.keys():
            print("              * {:}".format(i))
        
        txt = input("""Input e.g. 1, ModLeon, LinearElastic, all, 1 2 3 (hit enter to continue): """)
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
        print("""Choose the element libraries you would like to install (e.g. 1, 1 2 3, or Enter to quit):
            (1) displacement elements
            (2) large displacement elements (total/updated Lagrangian)
            (3) gradient-enhanced elements
            (4) gradient-enhanced large displacement elements
            (5) Cosserat elements
            
        
        or simply type the name of your desired element from the following list: """)
        for i in dictElements.keys():
            print("              * {:}".format(i))
            
        txt = input("""Input e.g. 1, uelDisplacement, all, 1 2 3 (hit enter to continue): """)
        userList = txt.split()
        
        if not txt:
            break
        
        try:
            for i in userList: elementInstall(i)
        except:
            print(" no option chosen")
            
    print("""
+-----------------------------------------------------------------------------+
| Would you like to build the bftUserLibrary?                                 |
+-----------------------------------------------------------------------------+            """)
    while True: 
        txt = input("""yes or no (hit enter to continue): """)
      
        if not txt:
            break
        if strtobool(txt):
            shutil.rmtree("CMakeFiles")
            os.remove("CMakeCache.txt")
            os.system("cmake -DCMAKE_CXX_COMPILER=g++ .")
            os.system("make")
            
    print("""
+-----------------------------------------------------------------------------+
| Thank you for setting up the bftUserLibrary !                               |
+-----------------------------------------------------------------------------+            """)
