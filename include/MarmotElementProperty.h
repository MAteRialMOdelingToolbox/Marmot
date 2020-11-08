#pragma once

namespace userLibrary {
    enum MaterialCode : int;
}

class MarmotMaterialSection  {
  public:
    userLibrary::MaterialCode materialCode;
    const double*             materialProperties;
    int                       nMaterialProperties;

    MarmotMaterialSection( userLibrary::MaterialCode materialCode,
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
