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
 * @file src/lib/PeakParam.cpp
 * @date 09/17/19
 * @author user
 */

#include "LE3_2D_MASS_WL_PEAK_COUNT/PeakParam.h"

using namespace Euclid::WeakLensing::TwoDMass;
// Input namespace and classes
using namespace dpd::le3::wl::twodmass::inp::paramspeakcountconvergence;
using namespace dpd::le3::wl::twodmass::inp::paramspeakcountmassap;

static Elements::Logging logger = Elements::Logging::getLogger("PeakParam");

namespace LE3_2D_MASS_WL_PEAK_COUNT {

PeakParam::PeakParam(): m_nbScales(0), m_MinPeakThresh(0.0){ }

PeakParam::PeakParam(double minPeakThreshold, std::vector <double> &ApPeakRadius, int nbScales):
                  m_MinPeakThresh(minPeakThreshold), m_nbScales(nbScales), m_ApPeakRadius(ApPeakRadius)
{ }

PeakParam PeakParam::readPeakConvXML(const std::string& paramPeakConvergence){
 // print the filename from which parameters are going to be read
  logger.info() << "Getting information from input Parameter XML file " << paramPeakConvergence << " ...";
  // Parse the XML file and create the binding object
  logger.debug() << "Parsing file " << paramPeakConvergence << " ...";
  try {

   auto PeakConvParam_xml = dpd::le3::wl::twodmass::inp::paramspeakcountconvergence::
       DpdTwoDMassParamsPeakCatalogConvergence (paramPeakConvergence, xml_schema::flags::dont_validate);
   logger.debug("Binding object created successfully");
   m_nbScales = PeakConvParam_xml->Data().NPeakScale();
   m_MinPeakThresh = PeakConvParam_xml->Data().MinPeakThresh();
  } catch (const xml_schema::exception& e) {
   std::cerr << e << std::endl;
  }
return PeakParam( m_MinPeakThresh,  m_ApPeakRadius, m_nbScales );
}

PeakParam PeakParam::readPeakMassApertXML(const std::string& paramPeakMassAper){
// print the filename from which parameters are going to be read
  logger.info() << "Getting information from input Parameter XML file " << paramPeakMassAper << " ...";
  // Parse the XML file and create the binding object
  logger.debug() << "Parsing file " << paramPeakMassAper << " ...";
  try {
   auto PeakMassApParam_xml = dpd::le3::wl::twodmass::inp::paramspeakcountmassap::
            DpdTwoDMassParamsPeakCatalogMassAperture (paramPeakMassAper, xml_schema::flags::dont_validate);
   logger.debug("Binding object created successfully");

   auto& ApPeakRadius_list = PeakMassApParam_xml->Data().ApPeakRadius();
   logger.info() << "The list contains " << ApPeakRadius_list.size()<< " items";
   for (auto &item : ApPeakRadius_list) {
    logger.info() << "item value: " << item;
    m_ApPeakRadius.push_back(item);
   }

   m_MinPeakThresh = PeakMassApParam_xml->Data().MinPeakThresh();
  } catch (const xml_schema::exception& e) {
   std::cerr << e << std::endl;
  }
   m_nbScales =0;
return PeakParam( m_MinPeakThresh, m_ApPeakRadius, m_nbScales);
}

int PeakParam::getnbScales(){ return m_nbScales; }

double PeakParam::getMinPeakThresh(){ return m_MinPeakThresh; }

std::vector<double> PeakParam::getApPeakRadius(){ return m_ApPeakRadius; }

void readPeakParamFile (const boost::filesystem::path& PeakParamFile, LE3_2D_MASS_WL_PEAK_COUNT::PeakParam &params) {
  // Case 1: check file is of XML format
   if (true == checkFileType(PeakParamFile.native(), Euclid::WeakLensing::TwoDMass::signXML)) {
    // read XML, get parameter file type
    if (true == fileHasField(PeakParamFile.native(), "DpdTwoDMassParamsPeakCatalogMassAperture")) {
      logger.info() << "Parameter file found for Aperture Mass Peak count";
      params.readPeakMassApertXML(PeakParamFile.native());
    }
    if (true == fileHasField(PeakParamFile.native(), "DpdTwoDMassParamsPeakCatalogConvergence")) {
      logger.info() << "Parameter file found for Wavelet Peak count";
      logger.info()<<"Parameter file is for using Wavelet filters..";
      params.readPeakConvXML(PeakParamFile.native());
    }
   } else {
      throw Elements::Exception() << "Parameter file is not in XML format. . . .";
   }
}

} // LE3_2D_MASS_WL_PEAK_COUNT namespace
