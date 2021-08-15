/**
 * @file LE3_2D_MASS_WL_CARTESIAN/CartesianParam.h
 * @date 05/13/19
 * @author user
 *
 * @copyright (C) 2012-2020 Euclid Science Ground Segment
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
 *
 */

#ifndef _LE3_2D_MASS_WL_CARTESIAN_CARTESIANPARAM_H
#define _LE3_2D_MASS_WL_CARTESIAN_CARTESIANPARAM_H

#include <string>
#include <boost/filesystem.hpp>
#include "ElementsKernel/Logging.h"

#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"

// Datamodel for INPUT products
#include "ST_DataModelBindings/dpd/le3/wl/twodmass/inp/euc-test-le3-wl-twodmass-ParamsConvergencePatch.h"
#include "ST_DataModelBindings/dpd/le3/wl/twodmass/inp/euc-test-le3-wl-twodmass-ParamsConvergenceClusters.h"
#include "ST_DataModelBindings/dpd/le3/wl/twodmass/inp/euc-test-le3-wl-twodmass-ParamsConvergencePatchesToSphere.h"
#include "ST_DataModelBindings/dictionary/pro/le3/wl/twodmass/euc-test-le3-wl-twodmass.h"

namespace LE3_2D_MASS_WL_CARTESIAN {

/**
 * @class CartesianParam
 * @brief
 *
 */
class CartesianParam {

public:
  /**
   * @brief   Default empty constructor
   * @return  Empty object of class cartesianParam
  */
  CartesianParam();

  /**
   * @brief    Constructor from parameters
   * @param    <NItReducedShear> Number of iterations to perform ReducedShear algorithm [default: 0 (Runs KS)]
   * @param    <nbPatches> number of patches
   * @param    <PixelSize> pixelSize in arcminute
   * @param    <PatchWidth> MapSize/PatchWidth in degrees
   * @param    <mapCenterX> Center of Map at X-axis in degrees
   * @param    <mapCenterY> Center of Map at Y-axis in degrees
   * @param    <nbZBins> number of redshift bins
   * @param    <zMin> zMin
   * @param    <zMax> zMax
   * @param    <zMargin> zMargin
   * @param    <BalancedBins> parameter to trigger the equal number of gamaxies in a patch
   * @param    <NInpaint> Number of iterations to perform Inpainting algorithm [default: 0 (Runs KS)]
   * @param    <EqualVarPerScale> 1: set variance inside the mask is inforced to be equal to the variance outside
   *            the mask at different scales
   * @param    <ForceBModes> Set 1: to force B-modes to be zero inside the mask
   * @param    <nbScales> number of scales for the sigmaBounded constraint [Inpainting option: default is automatic]
   * @param    <add_borders> Set 0: not to add borders to the map / 1: to add borders to the map
   * @param    <RSsigmaGauss> Reduced shear denoising: Value of Sigma for gaussian filter
   * @param    <sigmaGauss> Final denoising: Value of Sigma for gaussian filter
   * @param    <ExtensionName> name of extension to write in outut fits file header (either KAPPA_PATCH or SNR_PATCH)
   * @param    <ParaFileType> Type of parameter file type used to write outut xml file
   *            (Conv_Patch or Conv_Cluster or Conv_PatchesToSphere)
   * @param    <massThreshold> Threshold in Mass
   * @param    <thresholdFDR> Final denoising: Threshold for MRLens filtering
   * @param    <nbSamples> to create noise maps
   * @param    <removeOffset> Set 0: not to remove offset from map / 1: to remove offset from map
   * @param    <squareMap> set 1: to force the gnomonic projection to be square

   * @details  <mapCenterY, mapCenterX, PatchWidth> are the substitute of <raMin, raMax, decMin, decMax, nbBinsX and
   *           nbBinsY> to extract the Map.
   * @return  Fully allocated parameters
  */
  CartesianParam(int NItReducedShear, int NPatches, float PixelSize, float PatchWidth, std::vector<double> mapCenterX,
           std::vector<double> mapCenterY, int nbZBins, std::vector<double> zMin, double zMax, double zMargin,
           long BalancedBins, int NInpaint, long EqualVarPerScale, long ForceBMode, int nbScales, long add_borders,
           float RSsigmaGauss, float sigmaGauss, std::string ExtensionName, std::string ParaFileType,
           double massThreshold=0.0, float thresholdFDR = 0.0, int nbSamples=0,
           long removeOffset = 0, bool squareMap = true);
  /**
   * @brief Destructor
   */
  virtual ~CartesianParam() = default;

  /**
   * @brief   function to read the convergence patches parameter file in XML wrt dpd
   * @param   <paramConvergencePatch> parameter filename in <string> format
   * @return  parameters from the file
  */
  CartesianParam ReadConvPatchXMLFile (const std::string& paramConvergencePatch);

  /**
   * @brief   function to read the convergence Clusters parameter file in XML wrt dpd
   * @param   <paramConvClusters> parameter filename in <string> format
   * @return  parameters from the file
  */
  CartesianParam readConvClustersXMLFile (const std::string& paramConvClusters);

  /**
   * @brief   function to read the parameter file to create convergence Patches to sphere in XML wrt dpd
   * @param   <paramConvClusters> parameter filename in <string> format
   * @param   <catRamin> Minimum RA value from input survey in double
   * @param   <catRamax> Maximum RA value from input survey in double
   * @param   <catDecmin> Minimum DEC value from input survey in double
   * @param   <catDecmax> Maximum DEC value from input survey in double
   * @return  parameters from the file
  */
  CartesianParam readConvPatchesToSphereXMLFile (const std::string& paramConvPatchesToSphere, double& catRamin,
                                                 double& catRamax, double& catDecmin, double& catDecmax);

public:

  /**
   * @brief   function to return zMin value
   * @return  Minimum Redshift (Z) value from Parameter file
  */
  std::vector<double> getZMin();

  /**
   * @brief   function to return zMax value
   * @return  Maximum Redshift (Z) value from Parameter file
  */
  double getZMax();

  /**
   * @brief   function to return PixelSize value
   * @return  PixelSize from Parameter file
  */
  float getPixelsize();

  /**
   * @brief   function to return Sigma value for Map from Parameter file
   * @return  Sigma value
  */
  float getSigmaGauss();

  /**
   * @brief   function to set Sigma value for Map
  */
  void SetSigmaGauss(float val);

  /**
   * @brief   function to return Sigma value for gaussian filtering in case of reduce shear
  */
  float getRSSigmaGauss();

  /**
   * @brief   function to return Threshold from Parameter file
   * @return  Threshold
  */
  float getThreshold();

  /**
   * @brief   function to return threshold in mass
  */
  double getMassThreshold();

  /**
   * @brief   function to return number of catalogs based on Redshift bins from Parameter file
   * @return  Redshift bins
  */
  int getnbZBins();

  /**
   * @brief   function to return number of catalogs based on number of galaxies from Parameter file
   * @return  number of catalogs based on number of galaxies
  */
  int getnbPatches();

  /**
   * @brief   function to return number of Scales from Parameter file
   * @return  number of scales
  */
  int getnbScales();

  /**
   * @brief   function to return Redshift margin
  */
  double getZMargin();

  /**
   * @brief   function to return Iteration numbers from Parameter file
   * @return  number of iteration to perform inpainting
  */
  int getNInpaint();

  /**
   * @brief   function to return Iteration numbers from Parameter file
   * @return  number of iteration to perform inpainting
  */
  int getNItReducedShear();

  /**
   * @brief   function to return remove offset
   * @return  remove offset setting to/not to remove offset
  */
  long get_removeOffset();

  /**
   * @brief   function to return add borders setting from Parameter file
   * @return  add borders setting to/not to add borders to map
  */
  long get_addBorders();

  /**
   * @brief   function to return sigma bound setting from Parameter file
   * @return  sigma bound setting to force variance inside the mask is to be equal to the variance outside
   *           the mask at different scales
  */
  long getEqualVarPerScale();

  /**
   * @brief   function to return BmodesZeros value from Parameter file
   * @return  SquareMap to/not to force B-modes to be zero inside the mask
  */
  long getForceBMode();

  /**
   * @brief   function to return SquareMap value (by default it's true)
   * @return  SquareMap to/not to force the gnomonic projection to be square
  */
  bool getSquareMap();

  /**
   * @brief   function to return MapSize/PatchWidth in degrees from Parameter file
   * @return  MapSize in degrees
  */
  float getPatchWidth();

  /**
   * @brief   function to return MapCenter at X-axis (Ra) from Parameter file
   * @return  MapCenter at X-axis
  */
  std::vector<double> getMapCenterX();

  /**
   * @brief   function to return MapCenter at Y-axis (Dec) from Parameter file
   * @return  MapCenter at Y-axis
  */
  std::vector<double> getMapCenterY();

  /**
   * @brief   function to return number of SNR maps needed from Parameter file
   * @return  Number of samples
  */
  int getNSamples();

  /**
   * @brief   function to return balanced number of galaxies
   * @return  1: balanced galaxies needed 0: not needed
  */
  long get_BalancedBins();

  /**
   * @brief   function to return Extension name to write on output fits file
   * @return  Extension name (SNR_PATCH or KAPPA_PATCH)
  */
  std::string getExtName();

  /**
   * @brief   function to set Extension name to write on output fits file
   * @return  Extension name (SNR_PATCH or KAPPA_PATCH ...)
  */
  void setExtName(std::string& name);

  /**
   * @brief   function to return ParaFileType to write an output xml file
   * @return  ParaFileType (Conv_Patch or Conv_Cluster or Conv_PatchesToSphere)
  */
  std::string getParaFileType();

  /**
  * @brief Returns the value of the RaMax
  */
 double getRaMax(double& ramin);
  /**
  * @brief Returns the value of the RaMin
  */
 double getRaMin(double& MapCenterX);
  /**
  * @brief Returns the value of the DecMin calculated using mapCenter
  */
 double getDecMin(double& MapCenterY);
  /**
  * @brief Returns the value of the DecMax
  */
 double getDecMax(double& decmin);
  /**
  * @brief Returns the value of the X axis
  */
 int getXaxis();
  /**
  * @brief Returns the value of the Y axis
  */
 int getYaxis();
  /**
  * @brief Returns the value of Nside for patches to sphere
  */
 int getNside();

private:
double m_zMargin, raMin, raMax, decMin, decMax, m_zMax, m_massThreshold, m_overlap;
float m_sigmaGauss, m_thresholdFDR, m_PatchWidth, m_PixelSize;
float m_RSsigmaGauss;
std::vector<double> mapCenterX;
std::vector<double> mapCenterY;
std::vector<double> m_zMin;
//std::vector<double> m_zMax;
int m_nbZBins, m_NInpaint, m_nbScales, m_NItReducedShear, m_nbPatches, m_nbSamples, xbin, ybin, m_Nside;
long m_removeOffset, m_add_borders, m_ForceBMode, m_EqualVarPerScale, m_balancedBin;
bool squareMap;
std::string ExtName, m_ParaFileType;

};  // End of CartesianParam class

  /**
  * @brief     reads the parameter file
  * @param     <ParamFile>, <boost::filesystem::path> Parameter filename with path
  * @param     <params>, <LE3_2D_MASS_WL_CARTESIAN::CartesianParam> object to return parameters
  * @param     <catRamin> Minimum RA value from input survey in double
  * @param     <catRamax> Maximum RA value from input survey in double
  * @param     <catDecmin> Minimum DEC value from input survey in double
  * @param     <catDecmax> Maximum DEC value from input survey in double
  */
 void readParameterFile (const boost::filesystem::path& ParamFile, LE3_2D_MASS_WL_CARTESIAN::CartesianParam &params,
                         double catRamin=0., double catRamax=0., double catDecmin=0., double catDecmax=0.);
/*
template<typename T>
T* initializeArray(const unsigned int nbBinsX, const unsigned int nbBinsY, const unsigned int nbBinsZ)
 {
   // Create the map array
   T *array = new T[nbBinsX*nbBinsY*nbBinsZ];
   // Initialize it to zeros
   for (unsigned int i = 0; i<nbBinsX*nbBinsY*nbBinsZ; i++)
    {
      array[i] = 0.;
    }
   return array;
 }*/
}  // namespace LE3_2D_MASS_WL_CARTESIAN

#endif
