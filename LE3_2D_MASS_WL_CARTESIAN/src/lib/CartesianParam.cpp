/**
 * @file src/lib/CartesianParam.cpp
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

#include "LE3_2D_MASS_WL_CARTESIAN/CartesianParam.h"
#include <cmath>

// Input namespace and classes
using namespace dpd::le3::wl::twodmass::inp::paramsconvergencepatch;
using namespace dpd::le3::wl::twodmass::inp::paramsconvergenceclusters;
using namespace dpd::le3::wl::twodmass::inp::paramsconvergencempatchestosphere;
using namespace Euclid::WeakLensing::TwoDMass;

static Elements::Logging logger = Elements::Logging::getLogger("CartesianParam");
namespace LE3_2D_MASS_WL_CARTESIAN {

  /**
   * @brief   constructor (to create an object with default values)
  */
 CartesianParam::CartesianParam():m_NItReducedShear(10), m_nbPatches(1), m_PixelSize( 0.586/60.), m_PatchWidth(10.),
           mapCenterX(0.), mapCenterY(0.), m_nbZBins(1), m_zMin(0.), m_zMax(10.), m_balancedBin(0),
           m_NInpaint(100), m_EqualVarPerScale(0), m_ForceBMode(1), m_nbScales(0), m_add_borders(0),
           m_sigmaGauss(0.), m_thresholdFDR(0.), m_nbSamples(0), m_removeOffset(0), squareMap(true),
           raMin(0.0), raMax(0.0), decMin(0.0), decMax(0.0), xbin(1024), ybin(1024), ExtName("KAPPA_PATCH"),
           m_zMargin(0.), m_RSsigmaGauss(0.), m_massThreshold(0.), m_ParaFileType("Conv_Patch")
 { }

 CartesianParam::CartesianParam(int NItReducedShear, int NPatches, float PixelSize, float PatchWidth,
           std::vector<double> mapCenterX, std::vector<double> mapCenterY, int nbZBins, std::vector<double> zMin,
           double zMax, double zMargin, long BalancedBins, int NInpaint, long EqualVarPerScale, long ForceBMode,
           int nbScales, long add_borders, float RSsigmaGauss, float sigmaGauss, std::string ExtensionName,
           std::string ParaFileType, double massThreshold, float thresholdFDR, int nbSamples,
           long removeOffset, bool squareMap): m_NItReducedShear(NItReducedShear), m_nbPatches(NPatches),
           m_PixelSize(0.586/60.), m_PatchWidth(10.), mapCenterX(mapCenterX), mapCenterY(mapCenterY),
           m_nbZBins(nbZBins), m_zMin(zMin), m_zMax(zMax), m_zMargin(zMargin), m_balancedBin(BalancedBins),
           m_NInpaint(NInpaint), m_EqualVarPerScale(EqualVarPerScale), m_ForceBMode(ForceBMode),
           m_nbScales(nbScales), m_add_borders(add_borders), m_RSsigmaGauss(RSsigmaGauss), m_sigmaGauss(sigmaGauss),
           m_thresholdFDR(thresholdFDR), m_nbSamples(nbSamples),
           m_removeOffset(removeOffset), squareMap(squareMap), m_massThreshold(massThreshold), raMin(0.0), raMax(0.0),
           decMin(0.0), decMax(0.0), xbin(1024), ybin(1024), ExtName(ExtensionName), m_ParaFileType(ParaFileType) { }

  /**
   * @brief   function to read Convergence Patches parameter XML file with respect to Data Model
  */
 CartesianParam CartesianParam::ReadConvPatchXMLFile (const std::string& paramConvergencePatch){
  // print the filename from which parameters are going to be read
  logger.info() << "Getting information from input Parameter XML file " << paramConvergencePatch << " ...";
  try {

   auto ConvParam_xml = dpd::le3::wl::twodmass::inp::paramsconvergencepatch::DpdTwoDMassParamsConvergencePatch
                        (paramConvergencePatch, xml_schema::flags::dont_validate);

   m_NInpaint = ConvParam_xml->Data().GapsParams().NInpaint();

   m_NItReducedShear = ConvParam_xml->Data().ReducedShear().NItReducedShear();
   m_RSsigmaGauss = ConvParam_xml->Data().ReducedShear().GaussSTD();

    // sometimes parameter can appear more than one time, its corresponding method returns a sequence
   // Iterate over individual records.
   pro::le3::wl::twodmass::twoDMassParamsConvergencePatch::PatchParams_sequence& pp
                                             (ConvParam_xml->Data().PatchParams());
   for (pro::le3::wl::twodmass::twoDMassParamsConvergencePatch::PatchParams_iterator i
                                                       (pp.begin()); i!=pp.end(); ++i) {
     pro::le3::wl::twodmass::patchParams& pt (*i);
     pro::le3::wl::twodmass::patchParams::PatchList_sequence& pl (pt.PatchList());
     for (pro::le3::wl::twodmass::patchParams::PatchList_iterator ii (pl.begin()); ii!=pl.end(); ++ii) {
       pro::le3::wl::twodmass::patchDefinition& ptr (*ii);
// In this case, it is list of patch width and pixel size
       m_PatchWidth = ptr.PatchWidth();
       logger.info()<< "PatchWidth: "<< ptr.PatchWidth();
       m_PixelSize = (ptr.PixelSize())/60.;
       mapCenterX.push_back(ptr.ProjCtr().Longitude());
       mapCenterY.push_back(ptr.ProjCtr().Latitude());
     }
     if (pt.NPatches().present()) {
        m_nbPatches = pt.NPatches().get();
     }
   }

   m_EqualVarPerScale = ConvParam_xml->Data().GapsParams().EqualVarPerScale();

   if (true == ConvParam_xml->Data().GapsParams().ForceBMode()) {
   m_ForceBMode = 1;
   } else {
   m_ForceBMode = 0;
   }

   m_nbScales = ConvParam_xml->Data().GapsParams().NInpScale();

   if (true == ConvParam_xml->Data().GapsParams().AddBorder()) {
   m_add_borders = 1;
   } else {
   m_add_borders = 0;
   }

   m_sigmaGauss = ConvParam_xml->Data().DenoiseParams().GaussSTD();

  // The present() method used to check if the element exists in the XML
  // The get() method must be used to retrieve the element.
   if (ConvParam_xml->Data().DenoiseParams().ThresholdFDR().present()) {
      m_thresholdFDR = ConvParam_xml->Data().DenoiseParams().ThresholdFDR().get();
   }

   m_nbZBins = ConvParam_xml->Data().RedshiftBins().Nbins();

   pro::le3::wl::redshiftBinList::RedshiftBin_sequence& rbin
                                             (ConvParam_xml->Data().RedshiftBins().RedshiftBin());

   for (pro::le3::wl::redshiftBinList::RedshiftBin_iterator i
                                                       (rbin.begin()); i!=rbin.end(); ++i) {
     pro::le3::wl::redshiftBin& pt (*i);
     m_zMax = pt.ZMax();
//     m_zMin = pt.ZMin();
//     m_zMax.push_back(pt.ZMax());
     m_zMin.push_back(pt.ZMin());
   }

   if (true == ConvParam_xml->Data().RedshiftBins().BalancedBins()) {
   m_balancedBin = 1;
   } else {
   m_balancedBin = 0;
   }
   m_nbSamples = ConvParam_xml->Data().NResamples();;
  } catch (const xml_schema::exception& e) {
   std::cerr << e << std::endl;
  }
 ExtName = "KAPPA_PATCH";
 m_ParaFileType = "Conv_Patch";
 m_massThreshold = 0.;
 m_zMargin = 0.;
 return CartesianParam (m_NItReducedShear, m_nbPatches, m_PixelSize, m_PatchWidth, mapCenterX, mapCenterY, m_nbZBins,
                  m_zMin, m_zMax, m_zMargin, m_balancedBin, m_NInpaint, m_EqualVarPerScale, m_ForceBMode, m_nbScales,
                  m_add_borders, m_RSsigmaGauss, m_sigmaGauss, ExtName, m_ParaFileType, m_massThreshold,
                  m_thresholdFDR, m_nbSamples, m_removeOffset, true);
 }

  /**
   * @brief   function to read the convergence Clusters parameter file with respect to Data Model
  */
 CartesianParam CartesianParam::readConvClustersXMLFile (const std::string& paramConvClusters){
  // print the filename from which parameters are going to be read
  logger.info() << "Getting information from input Parameter XML file " << paramConvClusters << " ...";
  // Parse the XML file and create the binding object
  logger.debug() << "Parsing file " << paramConvClusters << " ...";
  try {

   auto ConvCluster_xml = dpd::le3::wl::twodmass::inp::paramsconvergenceclusters::DpdTwoDMassParamsConvergenceClusters
                                            (paramConvClusters, xml_schema::flags::dont_validate);
   logger.debug("Binding object created successfully");

   m_NInpaint = ConvCluster_xml->Data().GapsParams().NInpaint();

   m_NItReducedShear = ConvCluster_xml->Data().ReducedShear().NItReducedShear();
   m_RSsigmaGauss = ConvCluster_xml->Data().ReducedShear().GaussSTD();

   m_EqualVarPerScale = ConvCluster_xml->Data().GapsParams().EqualVarPerScale();

   if (true == ConvCluster_xml->Data().GapsParams().ForceBMode()) {
   m_ForceBMode = 1;
   } else {
   m_ForceBMode = 0;
   }

   m_nbScales = ConvCluster_xml->Data().GapsParams().NInpScale();

   if (true == ConvCluster_xml->Data().GapsParams().AddBorder()) {
   m_add_borders = 1;
   } else {
   m_add_borders = 0;
   }
   m_sigmaGauss = ConvCluster_xml->Data().DenoiseParams().GaussSTD();

  // The present() method used to check if the element exists in the XML
  // The get() method must be used to retrieve the element.
   if (ConvCluster_xml->Data().DenoiseParams().ThresholdFDR().present()) {
      m_thresholdFDR = ConvCluster_xml->Data().DenoiseParams().ThresholdFDR().get();
   }
   m_nbSamples = ConvCluster_xml->Data().NResamples();
   m_massThreshold = ConvCluster_xml->Data().MassThreshold();
   m_zMargin = ConvCluster_xml->Data().ZMargin();
   m_zMax = ConvCluster_xml->Data().ZMax();
   //m_zMax.push_back(ConvCluster_xml->Data().ZMax());
   m_PixelSize = (ConvCluster_xml->Data().PixelSize())/60.;
   m_PatchWidth = ConvCluster_xml->Data().PatchWidth();
// In this case, it is just single values for patch width and pixel size
//   m_PatchWidth.push_back(ConvCluster_xml->Data().PatchWidth());
//   m_PixelSize.push_back((ConvCluster_xml->Data().PixelSize())/60.);
   m_zMin.push_back(0.0); // zMargin + zCluster(from catalog)
   mapCenterX [0.0]; // RA Cluster (from catalog)
   mapCenterY [0.0]; // DEC Cluster (from catalog)
   m_nbPatches = 0;
   ExtName = "KAPPA_PATCH";
   m_ParaFileType = "Conv_Cluster";
  } catch (const xml_schema::exception& e) {
   std::cerr << e << std::endl;
  }

  return CartesianParam (m_NItReducedShear, m_nbPatches, m_PixelSize,
                      m_PatchWidth, mapCenterX, mapCenterY, m_nbZBins,
                      m_zMin, m_zMax, m_zMargin, m_balancedBin, m_NInpaint,
                      m_EqualVarPerScale, m_ForceBMode, m_nbScales,
                      m_add_borders, m_RSsigmaGauss, m_sigmaGauss, ExtName, m_ParaFileType, m_massThreshold,
                      m_thresholdFDR, m_nbSamples, m_removeOffset, true);
 }

 CartesianParam CartesianParam::readConvPatchesToSphereXMLFile (const std::string& paramConvPatchesToSphere, 
               double& catRamin, double& catRamax, double& catDecmin, double& catDecmax){
  // print the filename from which parameters are going to be read
  logger.info() << "Getting information from input Parameter XML file " << paramConvPatchesToSphere << " ...";
  // Parse the XML file and create the binding object
  logger.debug() << "Parsing file " << paramConvPatchesToSphere << " ...";
  try {

   auto ConvToSphere_xml = 
        dpd::le3::wl::twodmass::inp::paramsconvergencempatchestosphere::DpdTwoDMassParamsConvergencePatchesToSphere
                                            (paramConvPatchesToSphere, xml_schema::flags::dont_validate);
   logger.debug("Binding object created successfully");

   m_NInpaint = ConvToSphere_xml->Data().GapsParams().NInpaint();

   m_NItReducedShear = ConvToSphere_xml->Data().ReducedShear().NItReducedShear();
   m_RSsigmaGauss = ConvToSphere_xml->Data().ReducedShear().GaussSTD();

   m_EqualVarPerScale = ConvToSphere_xml->Data().GapsParams().EqualVarPerScale();

   if (true == ConvToSphere_xml->Data().GapsParams().ForceBMode()) {
   m_ForceBMode = 1;
   } else {
   m_ForceBMode = 0;
   }

   m_nbScales = ConvToSphere_xml->Data().GapsParams().NInpScale();

   if (true == ConvToSphere_xml->Data().GapsParams().AddBorder()) {
   m_add_borders = 1;
   } else {
   m_add_borders = 0;
   }
   m_sigmaGauss = ConvToSphere_xml->Data().DenoiseParams().GaussSTD();

  // The present() method used to check if the element exists in the XML
  // The get() method must be used to retrieve the element.
   if (ConvToSphere_xml->Data().DenoiseParams().ThresholdFDR().present()) {
      m_thresholdFDR = ConvToSphere_xml->Data().DenoiseParams().ThresholdFDR().get();
   }
   m_nbSamples = ConvToSphere_xml->Data().NResamples();

   m_nbZBins = ConvToSphere_xml->Data().RedshiftBins().Nbins();

   pro::le3::wl::redshiftBinList::RedshiftBin_sequence& rbin
                                             (ConvToSphere_xml->Data().RedshiftBins().RedshiftBin());

   for (pro::le3::wl::redshiftBinList::RedshiftBin_iterator i
                                                       (rbin.begin()); i!=rbin.end(); ++i) {
     pro::le3::wl::redshiftBin& pt (*i);
     m_zMax = pt.ZMax();
     //m_zMin = pt.ZMin();
     //m_zMax.push_back(pt.ZMax());
     m_zMin.push_back(pt.ZMin());
   }

   m_PatchWidth = ConvToSphere_xml->Data().MultiPatchParams().PatchWidth();
   m_PixelSize = (ConvToSphere_xml->Data().MultiPatchParams().PixelSize())/60.;
 // In this case, it is just single values for patch width and pixel size
 //  m_PatchWidth.push_back(ConvToSphere_xml->Data().MultiPatchParams().PatchWidth());
 //  m_PixelSize.push_back((ConvToSphere_xml->Data().MultiPatchParams().PixelSize())/60.);
   m_Nside = ConvToSphere_xml->Data().MultiPatchParams().Nside();
   m_overlap = ConvToSphere_xml->Data().MultiPatchParams().Overlap();

   double d_overlap = 1 + m_overlap;
   double ddec = catDecmax - catDecmin;
   double dra = catRamax - catRamin;

   int nlat = int((ddec*d_overlap)/m_PatchWidth) + 1;
   double dtheta_lat = ddec/nlat;
   m_nbPatches = 0;

   for (size_t i=0; i<nlat; i++) {
     double lat = i*dtheta_lat + dtheta_lat/2. + catDecmin;
     // number of patches in longitude wrt overlap and size of the survey
     int nlon = int ((dra * std::cos((i*m_PatchWidth + m_PatchWidth/2)* M_PI/180.) * d_overlap) / m_PatchWidth) + 1;
     double dtheta_lon = (dra * std::cos((i*m_PatchWidth + m_PatchWidth/2)* M_PI/180.)) / nlon;
     for (size_t j = 0; j<nlon; j++) {
        double lon = j*dtheta_lon + dtheta_lon/2. + catRamin;
        mapCenterX.push_back(lon);
        mapCenterY.push_back(lat);
        logger.info()<< "CenterX: "<< lon <<" "<< "CenterY: "<< lat;
        m_nbPatches = m_nbPatches + 1;
     }
   }
   logger.info()<< "Number of Patches: "<< m_nbPatches;
   ExtName = "KAPPA_PATCHES";
   m_ParaFileType = "Conv_PatchesToSphere";
  } catch (const xml_schema::exception& e) {
   std::cerr << e << std::endl;
  }

  return CartesianParam (m_NItReducedShear, m_nbPatches, m_PixelSize,
                      m_PatchWidth, mapCenterX, mapCenterY, m_nbZBins,
                      m_zMin, m_zMax, m_zMargin, m_balancedBin, m_NInpaint,
                      m_EqualVarPerScale, m_ForceBMode, m_nbScales,
                      m_add_borders, m_RSsigmaGauss, m_sigmaGauss, ExtName, m_ParaFileType, m_massThreshold,
                      m_thresholdFDR, m_nbSamples, m_removeOffset, true);
 }

  /**
   * @brief   functions return values of the parameters based on object
  */
 std::vector<double> CartesianParam::getZMin(){ return m_zMin; }
 double CartesianParam::getZMax(){ return m_zMax; }
 float CartesianParam::getPixelsize(){ return m_PixelSize; }
 float CartesianParam::getSigmaGauss(){ return m_sigmaGauss; }
 void CartesianParam::SetSigmaGauss(float val) {
   m_sigmaGauss = val;
 }
 float CartesianParam::getThreshold(){return m_thresholdFDR; }
 float CartesianParam::getRSSigmaGauss(){ return m_RSsigmaGauss; }

 int CartesianParam::getnbZBins(){ return m_nbZBins; }
 int CartesianParam::getnbPatches(){ return m_nbPatches; }
 int CartesianParam::getNInpaint(){ return m_NInpaint; }
 int CartesianParam::getNItReducedShear(){ return m_NItReducedShear; }
 int CartesianParam::getnbScales(){ return m_nbScales; }
 float CartesianParam::getPatchWidth(){ return m_PatchWidth; }
 std::vector<double> CartesianParam::getMapCenterX(){ return mapCenterX; }
 std::vector<double> CartesianParam::getMapCenterY(){ return mapCenterY; }
 int CartesianParam::getNSamples() { return m_nbSamples;}
 double CartesianParam::getMassThreshold(){ return m_massThreshold; }
 double CartesianParam::getZMargin(){ return m_zMargin; }

 long CartesianParam::get_removeOffset() {
  return m_removeOffset;
 }
 long CartesianParam::get_addBorders(){
  return m_add_borders;
 }
 long CartesianParam::getEqualVarPerScale(){
  return m_EqualVarPerScale;
 }
 long CartesianParam::getForceBMode(){
  return m_ForceBMode;
 }
 long CartesianParam::get_BalancedBins(){
  return m_balancedBin;
 }
 bool CartesianParam::getSquareMap(){
  return squareMap;
 }
 std::string CartesianParam::getExtName(){
  return ExtName;
 }
 void CartesianParam::setExtName(std::string& name) {
  ExtName = name;
 }
 std::string CartesianParam::getParaFileType(){
  return m_ParaFileType;
 }
  /**
  * @brief Returns the value of the RaMax
  */
 double CartesianParam::getRaMax(double& ramin) {
  raMax = ramin + getPatchWidth();
  return raMax;
 }
  /**
  * @brief Returns the value of the RaMin
  */
 double CartesianParam::getRaMin(double& MapCenterX) {
  if (fabs(MapCenterX) > 0 ) {
   raMin = MapCenterX-(getPatchWidth()/2.);
  }
  return raMin;
 }
  /**
  * @brief Returns the value of the DecMin calculated using mapCenter
  */
 double CartesianParam::getDecMin(double& MapCenterY) {
  if (fabs(MapCenterY) > 0 ) {
   decMin = MapCenterY-(getPatchWidth()/2.);
  }
  return decMin;
 }
  /**
  * @brief Returns the value of the DecMax
  */
 double CartesianParam::getDecMax(double& decmin) {
  decMax = decmin + getPatchWidth();
  return decMax;
 }

 int CartesianParam::getXaxis(){
    if (fabs(getPatchWidth()) > 0 && fabs(getPixelsize()) > 0){
     xbin =  ceil(getPatchWidth()/getPixelsize());
    }
  return xbin;
 }
  /**
  * @brief Returns the value of the Y axis
  */
 int CartesianParam::getYaxis(){
    if (fabs(getPatchWidth()) > 0 && fabs(getPixelsize()) > 0){
     ybin = ceil(getPatchWidth()/getPixelsize());
    }
   return ybin;
 }

 int CartesianParam::getNside() {
  return m_Nside;
 }

void readParameterFile (const boost::filesystem::path& ParamFile, LE3_2D_MASS_WL_CARTESIAN::CartesianParam &params,
                        double catRamin, double catRamax, double catDecmin, double catDecmax) {
   logger.info() << "Using Paramfile: " << ParamFile << " as DM input product";
   //check file is of XML format
   if (true == checkFileType(ParamFile.native(), Euclid::WeakLensing::TwoDMass::signXML)) {
    if (true == fileHasField(ParamFile.native(), "DpdTwoDMassParamsConvergenceClusters")) {
      logger.info()<<"Parameter file is for Cluster Convergence Patches..";
        params.readConvClustersXMLFile (ParamFile.native());
        logger.info()<< "Done reading Convergence Clusters Parameter file. . . ";
    }
    // read XML, get parameters
    // check if it's to get a patch
    if (true == fileHasField(ParamFile.native(), "DpdTwoDMassParamsConvergencePatch")) {
        logger.info()<<"Parameter file is for Convergence Patch..";
        params.ReadConvPatchXMLFile (ParamFile.native());
        logger.info()<< "Done reading Convergence Patch Parameter file. . . ";
    }
    // check if parameter file is for patches to sphere
    if (true == fileHasField(ParamFile.native(), "DpdTwoDMassParamsConvergencePatchesToSphere")) {
        logger.info()<<"Parameter file is for Convergence Patches To Sphere..";
        params.readConvPatchesToSphereXMLFile (ParamFile.native(), catRamin, catRamax, catDecmin, catDecmax);
        logger.info()<< "Done reading Convergence Patches To Sphere Parameter file. . . ";
    }
   }
}
}  // namespace LE3_2D_MASS_WL_CARTESIAN
