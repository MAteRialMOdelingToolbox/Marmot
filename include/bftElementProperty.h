#pragma once

namespace userLibrary {
    enum MaterialCode : int;
}

class BftMaterialSection  {
  public:
    userLibrary::MaterialCode materialCode;
    const double*             materialProperties;
    int                       nMaterialProperties;

    BftMaterialSection( userLibrary::MaterialCode materialCode,
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
