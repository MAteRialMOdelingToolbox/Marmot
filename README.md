This is the centra library of team Matthias/Magdalena providing access to
- BftMaterials:  Constutitive Models compatible with
    - EdelweissFE
    - Abaqus (via Abaqus-UserLibraryInterface)
    - Plaxis
    - OpenSees (via OpenSeesBftMaterialWrapper)
    - mpFEM 
- BftUels: Finite Elements compatible with
    - EdelweissFE
    - Abaqus (via Abaqus-UserLibraryInterface)
    - mpFEM 
    
For compilation, please make sure that
- bftMechanics is present
- desired Materials and Elements are present
- a proper version of the Eigen library is at hand
- and CMakeLists.txt is configured