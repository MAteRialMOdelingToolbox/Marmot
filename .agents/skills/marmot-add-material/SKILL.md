---
name: marmot-add-material
description: >-
  Step-by-step procedure to implement, register, test, and document a new constitutive material model in Marmot.
  Use when creating or refactoring material formulations (small-strain hypoelastic, finite-strain, damage, plasticity, or gradient-enhanced models).
---

# Implementing a Constitutive Material Model in Marmot

## 1. Rules & Architecture
- **In-Place State Variables**: Store **strictly path-dependent** variables via `stateLayout.add("name", size)`. Access by in-place reference `stateLayout.getAs<double&>(state.stateVars, "name")`.
- **Early Exit on Zero Increments**: `if (dStrain.isZero(1e-14)) { dStress_dStrain = Cel; return; }`.
- **Exceptions & Logging**: Throw `MarmotExceptions.h` (`StressUpdateFailed`, `SolverConvergenceFailed`). Log via `MarmotJournal::warningToMSG(...)`.
- **Base Classes**:
  - Small-strain: `MarmotMaterialHypoElastic` (`MarmotMaterialHypoElasticFactory`)
  - Finite-strain: `MarmotMaterialFiniteStrain` (`MarmotMaterialFiniteStrainFactory`)
  - Nonlocal: `MarmotMaterialGeneralGradientEnhancedHypoElastic` (`MarmotMaterialGeneralGradientEnhancedHypoElasticFactory`)

---

## 2. Header Skeleton (`include/Marmot/<MyMaterial>.h`)

```cpp
#pragma once
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotTypedefs.h"

namespace Marmot::Materials {

  /**
   * @class MyMaterial
   * @brief Linear elastic / plastic constitutive formulation.
   */
  class MyMaterial : public MarmotMaterialHypoElastic {
  public:
    using MarmotMaterialHypoElastic::MarmotMaterialHypoElastic;
    MyMaterial( const double* materialProperties, int nMaterialProperties, int materialNumber );
    ~MyMaterial() override = default;

    double getDensity( const double* stateVars ) const override;
    void computeStress( state3D&                state,
                        Marmot::Matrix6d&       dStress_dStrain,
                        const Marmot::Vector6d& dStrain,
                        const timeInfo&         timeInfo ) const override;

  private:
    double m_E{ 0.0 };
    double m_nu{ 0.0 };
  };

} // namespace Marmot::Materials
```

---

## 3. Implementation Skeleton (`src/<MyMaterial>.cpp`)

```cpp
#include "Marmot/MyMaterial.h"
#include "Marmot/MarmotConstants.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotExceptions.h"

namespace Marmot::Materials {

  MyMaterial::MyMaterial( const double* props, int nProps, int label )
    : MarmotMaterialHypoElastic( props, nProps, label )
  {
    if ( this->nMaterialProperties < 2 ) {
      throw std::runtime_error( MakeString() << __PRETTY_FUNCTION__ << ": Insufficient material properties." );
    }
    m_E  = this->materialProperties[0];
    m_nu = this->materialProperties[1];

    stateLayout.add( "kappa", 1 );
    stateLayout.finalize();
  }

  double MyMaterial::getDensity( const double* /*stateVars*/ ) const
  {
    if ( this->nMaterialProperties < 3 ) throw std::runtime_error( MakeString() << __PRETTY_FUNCTION__ << ": No density." );
    return this->materialProperties[2];
  }

  void MyMaterial::computeStress( state3D& state, Matrix6d& dStress_dStrain, const Vector6d& dStrain, const timeInfo& /*tInfo*/ ) const
  {
    const auto Cel = ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( m_E, m_nu );
    if ( dStrain.isZero( 1e-14 ) ) {
      dStress_dStrain = Cel;
      return;
    }
    double& kappa = stateLayout.getAs< double& >( state.stateVars, "kappa" );
    state.stress += Cel * dStrain;
    dStress_dStrain = Cel;
  }

} // namespace Marmot::Materials
```

---

## 4. Factory Registration (`src/<MyMaterial>Registration.cpp`)

```cpp
#include "Marmot/MyMaterial.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"

namespace Marmot::Materials::Registration {
  using namespace MarmotLibrary;
  const static bool registered = MarmotMaterialHypoElasticFactory::registerMaterial< MyMaterial >( "MYMATERIAL" );
}
```

---

## 5. Workflow Links
- **Test**: [`marmot-create-test`](../marmot-create-test/SKILL.md) (analytical + `MarmotMathCore` tangent test).
- **Docs**: [`marmot-documentation`](../marmot-documentation/SKILL.md) (Sphinx `doc/pages/features/<name>.rst`).
- **QA**: [`marmot-code-review`](../marmot-code-review/SKILL.md).
