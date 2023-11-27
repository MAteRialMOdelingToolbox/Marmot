g++ \
    -o a.out \
    -I/home/ad/mambaforge3/include \
    -I/home/ad/constitutiveModelling/Marmot/modules/core/MarmotMechanicsCore/include/ \
    -I/home/ad/constitutiveModelling/Marmot/modules/core/MarmotMathCore/include/ \
    -I/home/ad/constitutiveModelling/Marmot/modules/core/MarmotFiniteStrainMechanicsCore/include \
    -I/home/ad/constitutiveModelling/Marmot/modules/core/MarmotMicromorphicCore/include \
    -I../include \
    -std=c++17 \
    test.cpp \
    ../src/CompressibleNeoHooke.cpp \
    /home/ad/constitutiveModelling/Marmot/modules/core/MarmotFiniteStrainMechanicsCore/src/MarmotMaterialFiniteStrain.cpp \
    /home/ad/constitutiveModelling/Marmot/src/*.cpp \
    /home/ad/constitutiveModelling/Marmot/modules/core/MarmotMathCore/src/*.cpp \


