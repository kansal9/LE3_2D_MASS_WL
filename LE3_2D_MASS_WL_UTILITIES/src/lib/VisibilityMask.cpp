/**
 * @file src/lib/VisibilityMask.cpp
 * @date 06/16/21
 * @author Vanshika Kansal
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

#include "LE3_2D_MASS_WL_UTILITIES/VisibilityMask.h"

static Elements::Logging logger = Elements::Logging::getLogger("VisibilityMask");

namespace Euclid {
 namespace WeakLensing {
  namespace TwoDMass {

VisibilityMask::VisibilityMask(int newNside): m_newNside(newNside) {}

void VisibilityMask::readVisibilityMask(fs::path& workdir, fs::path& MaskFile,
                      std::vector < std::vector < double> >& data) {
   if (true == MaskFile.string().empty()) {
     throw Elements::Exception() << "Input Mask filename is not found . . . ";
   }
   fs::path datadir {workdir / "data"};
   std::string inputFileName;
   if (true == checkFileType(datadir/MaskFile, Euclid::WeakLensing::TwoDMass::signFITS)) {
    if (false == fs::exists(datadir/MaskFile)) {
       throw Elements::Exception() << "Input data product " << datadir/MaskFile << " not found";
    }
     logger.info("Input Visibility Mask is in Fits format..");
     inputFileName = (datadir/MaskFile).string();
     getMaskData(inputFileName, data);
     logger.info("Done reading Input Visibility Mask");
   } else { //Case 2: When input file is in XML format
     if (false == fs::exists(workdir/MaskFile)) {
       throw Elements::Exception() << "Input data product " << workdir/MaskFile << " not found";
     }
     logger.info("Input Visibility Mask is in XML format..");
     getMaskData(workdir, MaskFile, data);
     logger.info("Done reading Input Visibility Mask");
   }
}

void VisibilityMask::getMaskData(fs::path& workdir, fs::path& InFile,
                        std::vector<std::vector<double> >& data) {
   fs::path datadir {workdir / "data"};
   // Case 1: When input file is in XML format (fetch filename from xml file)
   if (false == checkFileType(workdir /InFile, Euclid::WeakLensing::TwoDMass::signFITS)) {
    if (true == checkFileType(workdir /InFile, Euclid::WeakLensing::TwoDMass::signXML)) {
     logger.info(" Input file is in XML format..");
     logger.info() << "Using file " << workdir / InFile << " as DM input product";
     if (true == fileHasField((workdir /InFile).native(), "DpdTwoDMassVisibilityMask")) {
       logger.info() << "Input File type: FITS: Visibility Mask";
       // Read filename of the FITS file from the input XML file
       DmInput in_xml = DmInput::readVisibilityMaskXMLFile(workdir / InFile);
       // Get fits filename
       m_inputFile = (datadir /  in_xml.getFitsCatalogFilename()).string();
       //logger.info() << "Input File : " << m_inputFile;
       getMaskData(m_inputFile, data);
     }
    }
   }
}

std::string VisibilityMask::getMaskFitsFilename() {
  return m_inputFile;
}

void VisibilityMask::getMaskData(const std::string& filename, std::vector<std::vector<double> >& data) {
  if (true == checkFileType(fs::path(filename), Euclid::WeakLensing::TwoDMass::signFITS)){
     std::pair<Healpix_Map<double>, Healpix_Map<double> > mapPair;
     //mapPair = readHealpixMap(filename);
     mapPair = changeResolution(filename);
     int nside = mapPair.first.Nside();
     int npix = mapPair.first.Npix();
     logger.info() << "nside: " << nside;
     arr<double> myarr = mapPair.first.Map();
     std::vector<double> Indices;
     myarr.copyTo(Indices);
     arr<double> mValues = mapPair.second.Map();
     std::vector<double> mask_values;
     mValues.copyTo(mask_values);
     std::vector<double> ra;
     std::vector<double> dec;
     for (int id_pix=0; id_pix<npix; id_pix++){
          std::pair<double, double> radec = getIndex2RaDec(mapPair.first, id_pix);
          dec.push_back(radec.second);
          ra.push_back(radec.first);
     }
     data.push_back(ra);
     data.push_back(dec);
     data.push_back(Indices); //in case we need to use indicies
     data.push_back(mask_values);
  }
}

std::pair<Healpix_Map<double>, Healpix_Map<double> > VisibilityMask::changeResolution
                                     (const std::string& filename) {
   std::pair<Healpix_Map<double>, Healpix_Map<double> > mapPair;
   mapPair = readHealpixMap(filename);
   int nside = mapPair.first.Nside();

   Healpix_Map<double> outMapE;
   Healpix_Map<double> outMapB;
   outMapE.SetNside(m_newNside, RING);
   outMapB.SetNside(m_newNside, RING);
   outMapE.fill(0.);
   outMapB.fill(0.);

   if (nside < m_newNside) {
     logger.info() << "The upgrade of resolution begin . .";
     outMapE.Import_upgrade(mapPair.first);
     outMapB.Import_upgrade(mapPair.second);
     logger.info() << "The upgrade of resolution End . .";
   } else {
     if (nside > m_newNside) {
     logger.info() << "The degrade of resolution begin . .";
       outMapE.Import_degrade(mapPair.first);
       outMapB.Import_degrade(mapPair.second);
     logger.info() << "The degrade of resolution end . .";
     } else {
       outMapE.Import_nograde(mapPair.first);
       outMapB.Import_nograde(mapPair.second);
       logger.info() << "The current resolution is same as desired resolution . .";
     }
   }

  return std::pair<Healpix_Map<double>, Healpix_Map<double> > (outMapE, outMapB);
}

  } /* namespace TwoDMass */
 } /* namespace WeakLensing */
} /* namespace Euclid */
