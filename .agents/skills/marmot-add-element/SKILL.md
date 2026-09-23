---
name: marmot-add-element
description: >-
  Procedure for adding, registering, testing, and documenting a new finite element formulation in Marmot.
  Use when implementing continuum elements, structural elements, mixed/EAS elements, finite strain formulations, or micropolar elements.
---

# Implementing a Finite Element in Marmot

Elements in Marmot are **templated on spatial dimension and node count** (`< int nDim, int nNodes >`), inheriting from `MarmotElement` and `MarmotGeometryElement< nDim, nNodes >`.

---

## 1. Topologies & Node Numbering
- **1D Bar2**: `1, 2` | **2D Quad4**: `1..4` CCW | **2D Quad8**: Corners `1..4`, mid-sides `5..8`
- **3D Hexa8**: Bottom `1..4` CCW, top `5..8` CCW | **3D Hexa20**: Corners `1..8`, mid-sides `9..20`

---

## 2. Templated Header Skeleton (`include/Marmot/<MyElement>.h`)

```cpp
#pragma once
#include "Marmot/MarmotElement.h"
#include "Marmot/MarmotGeometryElement.h"
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotStateVarVectorManager.h"
#include "Marmot/MarmotTypedefs.h"
#include <memory>
#include <vector>
#include <string>

namespace Marmot::Elements {

  template < int nDim, int nNodes >
  class MyElement : public MarmotElement, public MarmotGeometryElement< nDim, nNodes > {
  public:
    enum SectionType { UniaxialStress, PlaneStress, PlaneStrain, Solid };
    using ParentGeometry = MarmotGeometryElement< nDim, nNodes >;
    using typename ParentGeometry::BSized;
    using typename ParentGeometry::JacobianSized;
    using typename ParentGeometry::XiSized;

    static constexpr int sizeLoadVector = nNodes * nDim;
    using RhsSized      = Eigen::Matrix< double, sizeLoadVector, 1 >;
    using KeSizedMatrix = Eigen::Matrix< double, sizeLoadVector, sizeLoadVector >;

    struct QuadraturePoint {
      const XiSized xi;
      const double  weight;
      double        detJ{ 0.0 };
      BSized        B{ BSized::Zero() };
      std::unique_ptr< MarmotMaterialHypoElastic > material;
      QuadraturePoint( XiSized xi_, double w_ ) : xi( xi_ ), weight( w_ ) {}
    };

    std::vector< QuadraturePoint > qps;
    const int elLabel;
    const SectionType sectionType;

    MyElement( int elementID, SectionType section = Solid );
    ~MyElement() override = default;

    int getNumberOfRequiredStateVars() override;
    std::vector< std::vector< std::string > > getNodeFields() override;
    std::vector< int > getDofIndicesPermutationPattern() override;
    int getNNodes() override { return nNodes; }
    int getNSpatialDimensions() override { return nDim; }
    int getNDofPerElement() override { return sizeLoadVector; }
    std::string getElementShape() override { return ParentGeometry::getElementShape(); }

    void assignStateVars( double* stateVars, int nStateVars ) override;
    void assignNodeCoordinates( const double* coordinates ) override {
      ParentGeometry::assignNodeCoordinates( coordinates );
    }
    void initializeYourself() override;
    void setInitialConditions( StateTypes state, const double* values ) override;

    void computeYourself( const double* QTotal, const double* dQ, double* Pint, double* K,
                          const double* time, double dT, double& pNewdT ) override;
    void computeDistributedLoad( DistributedLoadTypes loadType, double* Pext, double* K, int elementFace,
                                 const double* load, const double* QTotal, const double* time, double dT ) override;
    void computeBodyForce( double* Pext, double* K, const double* load, const double* QTotal,
                           const double* time, double dT ) override;
  };

} // namespace Marmot::Elements
```

---

## 3. Best Practices & Factory Registration

- **In-Place `stateVars`**: Slice contiguous buffers into quadrature points using `QPStateVarManager` / `MarmotStateVarVectorManager`.
- **Zero Allocations**: Use `BSized`, `JacobianSized`, `RhsSized`, `KeSizedMatrix` from `MarmotGeometryElement` in quadrature loops.
- **Registration (`src/<MyElement>Registration.cpp`)**:
  ```cpp
  #include "Marmot/MyElement.h"
  #include "Marmot/MarmotElementFactory.h"

  namespace Marmot::Elements::Registration {
    using namespace MarmotLibrary;
    template < class T >
    MarmotElementFactory::elementFactoryFunction makeFactoryFunction() {
      return []( int id ) -> MarmotElement* { return new T( id ); };
    }
    const static bool CPS4_reg = MarmotElementFactory::registerElement( "CPS4", makeFactoryFunction< MyElement< 2, 4 > >() );
    const static bool C3D8_reg = MarmotElementFactory::registerElement( "C3D8", makeFactoryFunction< MyElement< 3, 8 > >() );
  }
  ```

---

## 4. Workflow Links
- **Test**: [`marmot-create-test`](../marmot-create-test/SKILL.md) (patch test, rigid body modes, analytical checks).
- **Docs**: [`marmot-documentation`](../marmot-documentation/SKILL.md) (`doc/pages/features/<element>.rst`).
- **QA**: [`marmot-code-review`](../marmot-code-review/SKILL.md).
