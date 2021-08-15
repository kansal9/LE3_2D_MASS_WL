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
 * @file src/lib/SphericalParam.cpp
 * @date 07/14/19
 * @author user
 */

#include "LE3_2D_MASS_WL_SPHERICAL/SphericalParam.h"

// Input namespace and classes
using namespace dpd::le3::wl::twodmass::inp::paramsconvergencesphere;
using namespace dpd::le3::wl::twodmass::inp::paramsconvergencempatchestosphere;

using namespace Euclid::WeakLensing::TwoDMass;

static Elements::Logging logger = Elements::Logging::getLogger("SphericalParam");

namespace LE3_2D_MASS_WL_SPHERICAL {

 SphericalParam::SphericalParam():m_nside(2048), m_NItReducedShear(0), m_NInpaint(10), m_BmodesZeros(0),
                                  m_EqualVarPerScale(0), m_NInpScales(1), m_Nbins(1), m_Zmin(0.0), m_Zmax(10.0),
                                  m_balancedBins(1), m_sigmaGauss(0.0), m_thresholdFDR(0.0), m_NResamples(0),
                                  m_RSsigmaGauss(0.0), ExtName("KAPPA_SPHERE")
 {}

 SphericalParam::SphericalParam (int nside, int NItReducedShear, int NInpaint, long BmodesZeros, long EqualVarPerScale,
                   int NInpScales, int Nbins, double Zmin, double Zmax, long balancedBins, std::string ExtensionName,
                   float RSsigmaGauss, float sigmaGauss, float threshold, int NResamples): m_nside(nside),
                   m_NItReducedShear(NItReducedShear), m_NInpaint(NInpaint), m_BmodesZeros(BmodesZeros),
                   m_EqualVarPerScale(EqualVarPerScale), m_NInpScales(NInpScales), m_Nbins(Nbins), m_Zmin(Zmin),
                   m_Zmax(Zmax), m_balancedBins(balancedBins), m_RSsigmaGauss(RSsigmaGauss), m_sigmaGauss(sigmaGauss),
                   m_thresholdFDR(threshold), m_NResamples(NResamples), ExtName(ExtensionName)
 {}

 SphericalParam SphericalParam::getConvergenceSphereParam(const std::string& paramConvFile) {
  logger.info() << "Getting information from input Parameter XML file " << paramConvFile << " ...";
  // Parse the XML file and create the binding object
  logger.debug() << "Parsing file " << paramConvFile << " ...";
//std::auto_ptr <dpd::le3::wl::twodmass::inp::paramsconvergencesphere::DpdTwoDMassParamsConvergenceSphere>
//                                                                        paramConv_xml;
  try {
   auto paramConv_xml = dpd::le3::wl::twodmass::inp::paramsconvergencesphere::DpdTwoDMassParamsConvergenceSphere
                                                          (paramConvFile, xml_schema::flags::dont_validate);
   logger.debug("Binding object created successfully");

   m_NItReducedShear = paramConv_xml->Data().ReducedShear().NItReducedShear();
   m_RSsigmaGauss = paramConv_xml->Data().ReducedShear().GaussSTD();

   m_nside = paramConv_xml->Data().Nside();

   m_Nbins = paramConv_xml->Data().RedshiftBins().Nbins();

   pro::le3::wl::redshiftBinList::RedshiftBin_sequence& rbin
                                             (paramConv_xml->Data().RedshiftBins().RedshiftBin());

   for (pro::le3::wl::redshiftBinList::RedshiftBin_iterator i
                                                       (rbin.begin()); i!=rbin.end(); ++i) {
     pro::le3::wl::redshiftBin& pt (*i);
     m_Zmax = pt.ZMax();
     m_Zmin = pt.ZMin();
   }

   if (true == paramConv_xml->Data().RedshiftBins().BalancedBins()) {
     m_balancedBins = 1;
   } else {
     m_balancedBins = 0;
   }

   m_NInpaint = paramConv_xml->Data().GapsParams().NInpaint();

   if (true == paramConv_xml->Data().GapsParams().ForceBMode()) {
     m_BmodesZeros = 1;
   } else {
     m_BmodesZeros = 0;
   }
   if (true == paramConv_xml->Data().GapsParams().EqualVarPerScale()) {
     m_EqualVarPerScale = 1;
   } else {
     m_EqualVarPerScale = 0;
   }
   m_NInpScales = paramConv_xml->Data().GapsParams().NInpScale();

   m_NResamples = paramConv_xml->Data().NResamples();

   m_sigmaGauss = paramConv_xml->Data().DenoiseParams().GaussSTD();

   if (paramConv_xml->Data().DenoiseParams().ThresholdFDR().present()) {
      m_thresholdFDR = paramConv_xml->Data().DenoiseParams().ThresholdFDR().get();
   }

   ExtName = "KAPPA_SPHERE";

  } catch (const xml_schema::exception& e) {
   std::cerr << e << std::endl;
  }
  return SphericalParam (m_nside, m_NItReducedShear, m_NInpaint, m_BmodesZeros, m_EqualVarPerScale, m_NInpScales,
     m_Nbins,  m_Zmin, m_Zmax, m_balancedBins, ExtName, m_RSsigmaGauss, m_sigmaGauss, m_thresholdFDR, m_NResamples);
 }

  double SphericalParam::getZMin() {
   return m_Zmin;
  }
  double SphericalParam::getZMax() {
   return m_Zmax;
  }
  long SphericalParam::getBmodesZeros() {
   return m_BmodesZeros;
  }

  long SphericalParam::getEqualVarPerScale() {
   return m_EqualVarPerScale;
  }

  int SphericalParam::getNside() {
   return m_nside;
  }

  int SphericalParam::getNItReducedShear() {
   return m_NItReducedShear;
  }
  float SphericalParam::getRSSigmaGauss() {
   return m_RSsigmaGauss;
  }

  int SphericalParam::getNInpScales() {
   return m_NInpScales;
  }

  float SphericalParam::getSigmaGauss() {
   return m_sigmaGauss;
  }

  float SphericalParam::getThresholdFDR() {
   return m_thresholdFDR;
  }
  int SphericalParam::getNbins() {
   return m_Nbins;
  }
  long SphericalParam::getBalancedBins() {
   return m_balancedBins;
  }
  int SphericalParam::getNInpaint() {
   return m_NInpaint;
  }
  int SphericalParam::getNResamples() {
   return m_NResamples;
  }
 std::string SphericalParam::getExtName(){
  return ExtName;
 }
 void SphericalParam::setExtName(std::string& name) {
  ExtName = name;
 }
} // LE3_2D_MASS_WL_SPHERICAL namespace
