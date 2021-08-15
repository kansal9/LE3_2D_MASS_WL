/*
 * Copyright (C) 2012-2020 Euclid Science Ground Segment
 *
 * This library is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation; either version 3.0 of the License, or (at your option)
 * any later version.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

/**
 * @file LE3_2D_MASS_WL_SPHERICAL/SphericalParam.h
 * @date 07/14/19
 * @author user
 */

#ifndef _LE3_2D_MASS_WL_SPHERICAL_SPHERICALPARAM_H
#define _LE3_2D_MASS_WL_SPHERICAL_SPHERICALPARAM_H

#include <string>
#include <boost/filesystem.hpp>
#include "ElementsKernel/Logging.h"

#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"

#include "ST_DataModelBindings/dpd/le3/wl/twodmass/inp/euc-test-le3-wl-twodmass-ParamsConvergenceSphere.h"
#include "ST_DataModelBindings/dpd/le3/wl/twodmass/inp/euc-test-le3-wl-twodmass-ParamsConvergencePatchesToSphere.h"

namespace LE3_2D_MASS_WL_SPHERICAL {

/**
 * @class SphericalParam
 * @brief
 *
 */
class SphericalParam {

public:

  /**
   * @brief   Default empty constructor
   * @return  Empty object of class SphericalParam
  */
  SphericalParam();

  /**
   * @brief    Constructor from parameters
   * @param    <nside> map resolution Nside
   * @param    <NItReducedShear> Number of iterations to perform ReducedShear algorithm [default: 0 (Runs KS)]
   * @param    <NInpaint> Number of iterations to perform Inpainting algorithm [default: 0 (Runs KS)]
   * @param    <BmodesZeros> Set 1: to force B-modes to be zero inside the mask
   * @param    <EqualVarPerScale> Set 1: to force the same variance in the inpainted part and in the original data
   * @param    <NInpScales> number of scales for the sigmaBounded constraint [Inpainting option: default is automatic]
   * @param    <Nbins> Number of redshift bins
   * @param    <zMin> Minimum value or redshift value
   * @param    <zMax> Maximum value or redshift value
   * @param    <balancedBins> to force balanced bins inside maps
   * @param    <ExtensionName> name of extension to write in outut fits file header (either KAPPA_SPHERE or 
               GALCOUNT_SPHERE or SNR_SPHERE)
   * @param    <RSsigmaGauss> Value of Sigma for gaussian filter in case of reduced shear
   * @param    <sigmaGauss> Value of Sigma for gaussian filter which will be applied to Map
   * @param    <threshold> Threshold for filtering other than gaussian
   * @param    <NResamples> to create noise maps
   * @return   Fully allocated parameters
  */
  SphericalParam (int nside, int NItReducedShear, int NInpaint, long BmodesZeros, long EqualVarPerScale,
           int NInpScales, int Nbins, double Zmin, double Zmax, long balancedBins, std::string ExtensionName,
           float RSsigmaGauss, float sigmaGauss, float threshold, int NResamples);

  /**
   * @brief Destructor
   */
  virtual ~SphericalParam() = default;

  /**
   * @brief   function to read the parameter file in XML wrt dpd for convergence sphere
   * @param   <paramFile> parameter filename in <string> format
   * @return  parameters from the file
  */
  SphericalParam getConvergenceSphereParam(const std::string& paramConvFile);

  /**
   * @brief   function to read the parameter file in XML wrt dpd for Convergence Patches to sphere
   * @param   <paramFile> parameter filename in <string> format
   * @return  parameters from the file
  */
  //SphericalParam getConvergencePatchesToSphereParam(const std::string& paramConvToSphereFile);

  /**
   * @brief   function to return zMin value
   * @return  Minimum Redshift (Z) value from Parameter file
  */
  double getZMin();

  /**
   * @brief   function to return zMax value
   * @return  Maximum Redshift (Z) value from Parameter file
  */
  double getZMax();

  /**
   * @brief   function to return the BmodesZeros
   * @return  Map to/not to force B-modes to be zero inside the mask
  */
  long getBmodesZeros();

  /**
   * @brief   function to return the variance
   * @return  not/to be the same in the inpainted part and in the original data
  */
  long getEqualVarPerScale();

  /**
   * @brief   function to return the Nside
   * @return  resolution of the map
  */
  int getNside();

  /**
   * @brief   function to return the number of iterations for reduced shear
  */
  int getNItReducedShear();

  /**
   * @brief   function to return the number of iterations for Inpainting
  */
  int getNInpaint();

  /**
   * @brief   function to return number of Scales from Parameter file
  */
  int getNInpScales();

  /**
   * @brief   function to return Sigma value for gaussian filtering from Parameter file
  */
  float getSigmaGauss();

  /**
   * @brief   function to return Sigma value for gaussian filtering from Parameter file
  */
  float getRSSigmaGauss();

  /**
   * @brief   function to return threshold for filtering
  */
  float getThresholdFDR();

  /**
   * @brief   function to return number of redshift bins
  */
  int getNbins();

  /**
   * @brief   function to return BalancedBin
   * @return  to/not to force balanced bins inside maps
  */
  long getBalancedBins();

  /**
   * @brief   function to return number of samples for SNR maps
  */
  int getNResamples();

  /**
   * @brief   function to return Extension name to write on output fits file
   * @return  Extension name (KAPPA_SPHERE or GALCOUNT_SPHERE or SNR_SPHERE ....)
  */
  std::string getExtName();

  /**
   * @brief   function to set Extension name to write on output fits file
   * @return  Extension name (KAPPA_SPHERE or GALCOUNT_SPHERE or SNR_SPHERE ....)
  */
  void setExtName(std::string& name);
private:

int m_nside, m_NInpaint, m_NItReducedShear, m_NInpScales, m_Nbins, m_NResamples;
long m_BmodesZeros, m_EqualVarPerScale, m_balancedBins;
float m_sigmaGauss, m_thresholdFDR, m_RSsigmaGauss;
double m_Zmin, m_Zmax;
std::string ExtName;

}; /* End of SphericalParam class */

} /* namespace LE3_2D_MASS_WL_SPHERICAL */

#endif
