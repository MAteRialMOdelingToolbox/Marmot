/* ---------------------------------------------------------------------
 *                                       _
 *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_
 * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
 * | | | | | | (_| | |  | | | | | | (_) | |_
 * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
 *
 * Unit of Strength of Materials and Structural Analysis
 * University of Innsbruck,
 * 2020 - today
 *
 * festigkeitslehre@uibk.ac.at
 *
 * Alexandros Stathas alexandros.stathas@boku.ac.at
 * Matthias Neuner matthias.neuner@uibk.ac.at
 *
 * This file is part of the MAteRialMOdellingToolbox (marmot).
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of marmot.
 * ---------------------------------------------------------------------
 */

#pragma once
#include "Fastor/Fastor.h"
#include "Marmot/MarmotMaterial.h"

/**
 *
 * Derived abstract base class for elastic materials expressed purely in rate form.
 *
 * In general, the nominal stress rate tensor \f$ \sigRate \f$ can be written as a function of the nominal stress tensor
 * \f$ \sig \f$, the stretching rate tensor \f$ \epsRate \f$ and the time \f$ t \f$.
 *
 * \f[  \displaystyle \sigRate = f( \sig, \epsRate, t, ...) \f]
 *
 * In course of numerical time integration, this relation will be formulated incrementally as
 *
 * \f[  \displaystyle \Delta \sig = f ( \sig_n, \Delta\eps, \Delta t, t_n, ...) \f]
 *
 * with
 *
 * \f[  \displaystyle \Delta\eps =  \epsRate\, \Delta t \f]
 *
 * and the algorithmic tangent
 *
 * \f[ \displaystyle \frac{d \sig }{d \eps } =  \frac{d \Delta \sig }{d \Delta \eps } \f]
 *
 * This formulation is compatible with an Abaqus interface.
 */
class MarmotMaterialHypoElasticInterface : public MarmotMaterial {

public:
  using MarmotMaterial::MarmotMaterial;
  using Tensor1D = Fastor::Tensor< double, 3 >;
  using Tensor2D = Fastor::Tensor< double, 3, 3 >;

  /// Characteristic element length
  double characteristicElementLength;
  /**
   * Set the characteristic element length at the considered quadrature point.
   * It is needed for the regularization of materials with softening behavior based on the mesh-adjusted softening
   * modulus.
   *
   * @param[in] length characteristic length; will be assigned to @ref characteristicElementLength
   */
  void setCharacteristicElementLength( double length );

  /**
   * For a given linearized strain increment \f$\Delta\boldsymbol{\varepsilon}\f$ at the old and the current time,
   * compute the Cauchy stress and the algorithmic tangent
   * \f$\frac{\partial\boldsymbol{\sigma}^{(n+1)}}{\partial\boldsymbol{\varepsilon}^{(n+1)}}\f$.
   *
   * @param[in,out]	force           conjugate force due to displacement jump
   * @param[in,out]	surface_stress  conjugate surface_stress due to conjugate average surface strain
   * @param[in,out]	H_inv_ij	tangent representing the derivatives of the linearized displacement jump (H_inv_{ij})
   * @param[in,out]	Z_ijkl	  tangent representing the derivatives of the linearized average surface strain (Z_{ijkl})
   * @param[in,out]	H_inv_nF_ijk tangent of the consistency terms
   * @param[in,out]	Yn_H_inv_Fn_ijkl representing the derivatives of the linearized average surface strain
   * (Yn_H_inv_Fn_{ijkl})
   * @param[in]	dU linearized displacement increment on "top (0)" and "bottom (1)" sides of the interface
   * @param[in] dSurface_strain linearized surface strain increment on top and bottom sides of the interface
   * @param[in] normal The normal to the interface positive to the direction of the "top (0)" side
   * @param[in]	timeOld	Old (pseudo-)time
   * @param[in]	dt	(Pseudo-)time increment from the old (pseudo-)time to the current (pseudo-)time
   * @param[in,out]	pNewDT	Suggestion for a new time increment
   */
  virtual void computeStress( double*       force,
                              double*       surface_stress,
                              double*       H_inv_ij,
                              double*       Z_ijkl,
                              double*       H_inv_nF_ijk,
                              double*       Yn_H_inv_Fn_ijkl,
                              const double* dU,
                              const double* dSurface_strain,
                              const double* normal,
                              const double* timeOld,
                              const double  dT,
                              double&       pNewDT ) = 0;

  // virtual void computeStress(
  //                            Tensor1D&  force,
  //                            Tensor2D&  surface_stress,
  //                            Fastor::Tensor<double, 21,21>& dStress_dStrain,
  //                            const Fastor::Tensor<double, 6,1>& dU,
  //                            const Fastor::Tensor<double, 18,1>& dSurface_strain,
  //                            const Tensor1D& normal,
  //                            const double* timeOld,
  //                            const double  dT,
  //                            double&       pNewDT );
};

namespace MarmotLibrary {

  /**
   * @class MarmotMaterialHypoElasticInterfaceFactory
   * @brief Factory class for creating material instances.
   *
   * This class provides a mechanism to register materials by name,
   * and to create material instances based on their properties.
   * It allows for dynamic material creation without hardcoding specific material types.
   */
  class MarmotMaterialHypoElasticInterfaceFactory {
  public:
    using materialFactoryFunction = MarmotMaterialHypoElasticInterface* (*)( const double* materialProperties,
                                                                             int           nMaterialProperties,
                                                                             int           materialNumber );
    MarmotMaterialHypoElasticInterfaceFactory() = delete;

    /**
     * @brief Create a material instance based on its code and properties.
     * @param[in] materialName Name of the material.
     * @param[in] materialProperties Array of material properties.
     * @param[in] nMaterialProperties Number of properties in the array.
     * @param[in] materialNumber Unique identifier for the material instance.
     * @return Pointer to the created MarmotMaterialHypoElasticInterface instance, or nullptr if creation failed.
     */
    static MarmotMaterialHypoElasticInterface* createMaterial( const std::string& materialName,
                                                               const double*      materialProperties,
                                                               int                nMaterialProperties,
                                                               int                materialNumber );

    /**
     * @brief Register a material with its code and factory function.
     * @param[in] materialName Name of the material.
     * @param[in] factoryFunction Function to create material instances.
     * @return True if registration was successful, false if the code already exists.
     */
    static bool registerMaterial( const std::string& materialName, materialFactoryFunction factoryFunction );

  private:
    static std::unordered_map< std::string, materialFactoryFunction > materialNameToModelAssociation;
  };
} // namespace MarmotLibrary

/** @brief Creates a default factory function for MarmotMaterialHypoElasticInterface instances.
 * @tparam T The type of the material to be created.
 * @return A factory function that creates an instance of T using the provided material properties.
 */
template < typename T >
MarmotLibrary::MarmotMaterialHypoElasticInterfaceFactory::materialFactoryFunction makeDefaultMarmotMaterialHypoElasticInterfaceFactoryFunction()
{
  return []( const double* materialProperties,
             int           nMaterialProperties,
             int           materialNumber ) -> MarmotMaterialHypoElasticInterface* {
    return new T( materialProperties, nMaterialProperties, materialNumber );
  };
}
