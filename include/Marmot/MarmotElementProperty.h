#pragma once

namespace MarmotLibrary {
    enum MaterialCode : int;
}

class MarmotMaterialSection  {
  public:
    MarmotLibrary::MaterialCode materialCode;
    const double*             materialProperties;
    int                       nMaterialProperties;

    MarmotMaterialSection( MarmotLibrary::MaterialCode materialCode,
                        const double*             materialProperties,
                        int                       nMaterialProperties )
        : materialCode( materialCode ),
          materialProperties( materialProperties ),
          nMaterialProperties( nMaterialProperties ){};
};

class ElementProperties  {
  public:
    const double* elementProperties;
    int           nElementProperties;

    ElementProperties( const double* elementProperties, int nElementProperties )
        : elementProperties( elementProperties ),
          nElementProperties( nElementProperties ){};
};
