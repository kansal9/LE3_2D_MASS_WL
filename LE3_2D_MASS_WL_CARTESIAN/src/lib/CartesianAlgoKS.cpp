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
 * @file src/lib/CartesianAlgoKS.cpp
 * @date 10/21/19
 * @author user
 */

#include "LE3_2D_MASS_WL_CARTESIAN/CartesianAlgoKS.h"
#include "ElementsKernel/Temporary.h"
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

using namespace Euclid::WeakLensing::TwoDMass;
using namespace LE3_2D_MASS_WL_CARTESIAN;
// Output namespace and classes
using namespace dpd::le3::wl::twodmass::out::convergencepatch;
 // handle on created path names
 //boost::filesystem::path temp_path;
static Elements::Logging logger = Elements::Logging::getLogger("CartesianAlgoKS");
namespace Euclid {
 namespace WeakLensing {
  namespace TwoDMass {
   namespace CartesianKS {

CartesianAlgoKS::CartesianAlgoKS(LE3_2D_MASS_WL_CARTESIAN::CartesianParam &params):
                          m_cartesianParam(params) {
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
 // Extract ShearMap function
////////////////////////////////////////////////////////////////////////////////////////////////////////
bool CartesianAlgoKS::extractShearMap(const std::string& shearMap, std::vector<std::vector<double> >& Data,
                                      LE3_2D_MASS_WL_CARTESIAN::CoordinateBound& CB) {
 LE3_2D_MASS_WL_CARTESIAN::ShearMap *m_ShearMap = nullptr;
 MapMaker map(Data, m_cartesianParam);
 m_ShearMap = map.getShearMap(CB);
 // Pixelate X and Y axis
 if ((m_ShearMap->getXdim()) == 2048 && (m_ShearMap->getYdim()) == 2048) {
  m_ShearMap->pixelate(1, 1);
 }
 if ((m_ShearMap->getXdim()) == 4096 && (m_ShearMap->getYdim()) == 4096) {
  m_ShearMap->pixelate(2, 2);
 }
 // Writing Shear Map
 if (m_ShearMap != nullptr) {
  std::string name="SHEAR_PATCH";
  m_cartesianParam.setExtName(name);
  m_ShearMap->writeMap(shearMap, m_cartesianParam);
 } else {
  return false;
 }

 delete m_ShearMap;
 m_ShearMap = nullptr;

 return true;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////
 // Extract ConvergenceMap function NOT REQUIRED (Delete)
////////////////////////////////////////////////////////////////////////////////////////////////////////
bool CartesianAlgoKS::extractConvergnceMap (const std::string& convergenceMap,
                      std::vector<std::vector<double> >& Data, LE3_2D_MASS_WL_CARTESIAN::CoordinateBound& CB) {
 LE3_2D_MASS_WL_CARTESIAN::ConvergenceMap *m_ConvergenceMap = nullptr;
 MapMaker map(Data, m_cartesianParam);

 m_ConvergenceMap = map.getConvMap(CB);
 // Pixelate X and Y axis
 if ((m_ConvergenceMap->getXdim()) == 2048 && (m_ConvergenceMap->getYdim()) == 2048) {
  m_ConvergenceMap->pixelate(1, 1);
 }
 if ((m_ConvergenceMap->getXdim()) == 4096 && (m_ConvergenceMap->getYdim()) == 4096) {
  m_ConvergenceMap->pixelate(2, 2);
 }
 // Writing convergence Map
 if (m_ConvergenceMap != nullptr) {
  std::string name="KAPPA_PATCH";
  m_cartesianParam.setExtName(name);
  m_ConvergenceMap->writeMap(convergenceMap, m_cartesianParam);
 } else {
  return false;
 }

 delete m_ConvergenceMap;
 m_ConvergenceMap = nullptr;

 return true;

}
////////////////////////////////////////////////////////////////////////////////////////////////////////
 // Perform KS Mass Mapping
////////////////////////////////////////////////////////////////////////////////////////////////////////
bool CartesianAlgoKS::performKSMassMapping(const std::string& shearMap, const std::string& outConvMap) {
 LE3_2D_MASS_WL_CARTESIAN::ConvergenceMap *m_ConvergenceMap = nullptr;
 MassMapping mass(m_cartesianParam);

 m_ConvergenceMap = mass.getSheartoConv(shearMap);

 // Writing Map
 if (m_ConvergenceMap != nullptr) {
  std::string name="KAPPA_PATCH";
  m_cartesianParam.setExtName(name);
  if (m_cartesianParam.getNInpaint() != 0 || outConvMap.find("NReSample") != std::string::npos) {
    m_ConvergenceMap->writeMap(outConvMap, m_cartesianParam);
  } else {
    m_ConvergenceMap->writeMap(outConvMap, m_cartesianParam);
  }
 } else {
  return false;
 }
 delete m_ConvergenceMap;
 m_ConvergenceMap = nullptr;

 return true;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
 // Perform Inverse KS Mass Mapping
////////////////////////////////////////////////////////////////////////////////////////////////////////
  bool CartesianAlgoKS::performInverseKSMassMapping(const std::string& convMap, const std::string& outShearMap) {
   LE3_2D_MASS_WL_CARTESIAN::ShearMap *m_ShearMap = nullptr;
   MassMapping mass(m_cartesianParam);

   m_ShearMap = mass.getConvtoShear(convMap);
   // Writing Shear Map
   if (m_ShearMap != nullptr) {
    m_ShearMap->writeMap(outShearMap, m_cartesianParam);
   } else {
    return false;
   }

   delete m_ShearMap;
   m_ShearMap = nullptr;
   return true;
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////
 // Perform Reduced Shear Computation
////////////////////////////////////////////////////////////////////////////////////////////////////////

bool CartesianAlgoKS::performReducedShear(boost::filesystem::path& InShearMap, boost::filesystem::path&
      outConvergenceMap, boost::filesystem::path& workdir, boost::filesystem::path outputConvMaps){
  std::vector<fs::path> filenames;
  fs::path datadir = workdir / "data" ;
  if (true == checkFileType((datadir/InShearMap).native(), Euclid::WeakLensing::TwoDMass::signFITS)) {
    if ((InShearMap.string()).empty() == false) {
      filenames.push_back(InShearMap);
  } } else {
    filenames = read_filenames(workdir, InShearMap);
  }
    int reducedIter = m_cartesianParam.getNItReducedShear();
    std::vector<fs::path> rsfilenames;
    for (size_t i = 0; i<filenames.size(); i++) {
      size_t pos = (filenames[i].string()).find("ShearMap");

      outConvergenceMap = fs::path("EUC_LE3_WL_ConvergenceMapKS_" + (filenames[i].string()).substr (pos));

      LE3_2D_MASS_WL_CARTESIAN::ConvergenceMap *m_ConvergenceMap = nullptr;
      LE3_2D_MASS_WL_CARTESIAN::ShearMap *m_reducedShear=nullptr;
      //m_ConvergenceMap = new LE3_2D_MASS_WL_CARTESIAN::ConvergenceMap(outConvergenceMap);
      m_reducedShear = new LE3_2D_MASS_WL_CARTESIAN::ShearMap((datadir/filenames[i]).native());

      for (int it = 0; it < reducedIter; it++) {
        fs::path reducedShearMap {};
        //reducedShearMap = fs::path("EUC_LE3_WL_ShearMap_reduced_"+ getDateTimeString() + ".fits");
        reducedShearMap = fs::path("EUC_LE3_WL_ShearMap_" + std::to_string(it) + "_Reduced" +
                           (filenames[i].string()).substr (pos));
       // if (m_cartesianParam.get_addBorders()){
        //  m_reducedShear->add_borders();
      //  }
        if (it > 0) {
          getReducedShearMap(*m_reducedShear, *m_ConvergenceMap);
        //  m_reducedShear->computeReducedShear(*m_ConvergenceMap);
        }

        /*if (m_cartesianParam.get_addBorders()) {
          m_reducedShear->remove_borders();
        }*/

        // Writing reduced Shear Map
        if (m_reducedShear != nullptr) {
          //reducedShearMap.clear();
          m_reducedShear->writeMap((datadir/reducedShearMap).native(), m_cartesianParam);
        } else {
          logger.info()<<"Reduced Shear Map array is empty . . .";
          return false;
        }
        fs::path RSMap = datadir /reducedShearMap;
        if (it == reducedIter-1) {
          //if ((reducedShearMap.string()).find("NReSample") != std::string::npos) {
          //  rsfilenames.push_back(reducedShearMap);
          //} else {
            std::ofstream outfile;
            outfile.open ((workdir / outputConvMaps).string(), std::ios_base::app);
            outfile << "[";
            outConvergenceMap.clear();
            getNoisyConvergenceMap(workdir, reducedShearMap, outConvergenceMap);
            outfile << outConvergenceMap.filename();
            outConvergenceMap.clear();
            getDenoisedConvergenceMap(workdir, reducedShearMap, outConvergenceMap);
            outfile << ",";
            outfile << outConvergenceMap.filename();
            outfile << "]";
            outfile.close();
          //}
        } else {
          perform_MassMapping_Function( RSMap, outConvergenceMap, workdir);

          delete m_ConvergenceMap;
          m_ConvergenceMap = nullptr;

          m_ConvergenceMap = new LE3_2D_MASS_WL_CARTESIAN::ConvergenceMap((workdir/outConvergenceMap).native());
          /*if (fabs(m_cartesianParam.getRSSigmaGauss())>0.001){//Apply gaussian filter on the map
            m_ConvergenceMap->applyGaussianFilter(m_cartesianParam.getRSSigmaGauss());
          }
          double std_dev = m_ConvergenceMap->getSigma();
          logger.info() <<"std dev: " << std_dev;
          m_ConvergenceMap->thresholding(5.*std_dev);*/
          getTildeConvergence(*m_ConvergenceMap);
        }
      } // end of reduced shear iteration numbers
     } //end of number of maps iterations

 return true;
 }

void CartesianAlgoKS::getReducedShearMap(LE3_2D_MASS_WL_CARTESIAN::ShearMap& reducedShear,
                LE3_2D_MASS_WL_CARTESIAN::ConvergenceMap& ConvMap) {
        if (m_cartesianParam.get_addBorders()){
          reducedShear.add_borders();
        }

        reducedShear.computeReducedShear(ConvMap);

        if (m_cartesianParam.get_addBorders()) {
          reducedShear.remove_borders();
        }
}

void CartesianAlgoKS::getTildeConvergence(LE3_2D_MASS_WL_CARTESIAN::ConvergenceMap& ConvMap) {
  if (fabs(m_cartesianParam.getRSSigmaGauss())>0.001){//Apply gaussian filter on the map
     ConvMap.applyGaussianFilter(m_cartesianParam.getRSSigmaGauss());
  }
  double std_dev = ConvMap.getSigma();
  logger.info() <<"std dev: " << std_dev;
  ConvMap.thresholding(5.*std_dev);

}
////////////////////////////////////////////////////////////////////////////////////////////////////////
 // Perform Inpainting
////////////////////////////////////////////////////////////////////////////////////////////////////////
  bool CartesianAlgoKS::performInPainting(const std::string& convMap, const std::string& shearMap,
                                          const std::string& outConvMap) {
   // const std::string workdir = "/run/media/user/Backup_Drive/WL_Results/";
    LE3_2D_MASS_WL_CARTESIAN::ConvergenceMap *m_ConvergenceMap=nullptr;
    LE3_2D_MASS_WL_CARTESIAN::ShearMap *m_ShearMap=nullptr;
    m_ConvergenceMap = new LE3_2D_MASS_WL_CARTESIAN::ConvergenceMap(convMap);
    m_ShearMap = new LE3_2D_MASS_WL_CARTESIAN::ShearMap(shearMap);
    remove(convMap.c_str());
    if (m_cartesianParam.get_addBorders()){
     m_ShearMap->add_borders();
   //  const std::string shearBorder = workdir + "shearMap_withBorder.fits";
   //  m_ShearMap->writeMap(shearBorder, m_cartesianParam);
    }
    InpaintingAlgo myIPalgo(*m_ShearMap, *m_ConvergenceMap, m_cartesianParam);
    LE3_2D_MASS_WL_CARTESIAN::ConvergenceMap *myconvMap = myIPalgo.performInPaintingAlgo();

    if (myconvMap==nullptr)
    {
     delete m_ShearMap;
     delete m_ConvergenceMap;
     m_ShearMap = nullptr;
     m_ConvergenceMap = nullptr;
     return false;
    }
   //  const std::string convKSPBorder = workdir + "convergenceMapKSPlus_withBorder.fits";
   //  myconvMap->writeMap(convKSPBorder, m_cartesianParam);
    // In case borders were added, remove them
    if (m_cartesianParam.get_addBorders()) {
     myconvMap->remove_borders();
    }

    // Writing Convergence Map
    if (myconvMap != nullptr) {
     if (outConvMap.find("NReSample") != std::string::npos) {
       myconvMap->writeMap(outConvMap, m_cartesianParam);
     } else {
       myconvMap->writeMap(outConvMap, m_cartesianParam);
     }
    } else {
     return false;
    }

    delete m_ShearMap;
    delete m_ConvergenceMap;
    delete myconvMap;
    m_ShearMap = nullptr;
    m_ConvergenceMap = nullptr;
    myconvMap = nullptr;
    return true;
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////
 // Perform Mass Mapping Function
////////////////////////////////////////////////////////////////////////////////////////////////////////
 bool CartesianAlgoKS::perform_MassMapping_Function(boost::filesystem::path& InShearMap, boost::filesystem::path&
      outConvergenceMap, boost::filesystem::path& workdir){

     logger.info("KS conversion from Shear to Convergence Map");
     //performKSMassMapping((workdir/"data"/InShearMap).native(), (workdir/"data"/outConvergenceMap).native());
     performKSMassMapping(InShearMap.native(), (workdir/outConvergenceMap).native());

     // Perform inpainting

     if (m_cartesianParam.getNInpaint() != 0) {
      logger.info("entering inPainting()");
      fs::path IntermediateConvMap = outConvergenceMap;
      std::string str = outConvergenceMap.string();
      boost::replace_all(str, "KS", "KSPlus");
      remove(outConvergenceMap);
      outConvergenceMap = fs::path(str);
      // perform InPainting
//      if (performInPainting((workdir/"data"/IntermediateConvMap).native(), (workdir/"data"/InShearMap).native(),
//                                            (workdir/"data"/outConvergenceMap).native()) == false)
      if (performInPainting((workdir/IntermediateConvMap).native(), InShearMap.native(),
                                            (workdir/outConvergenceMap).native()) == false)
      {
        //return Elements::ExitCode::OK;
        return false;
      }
     }
 return true;
 }

////////////////////////////////////////////////////////////////////////////////////////////////////////
 // set name of the output files and perform Mass Mapping
////////////////////////////////////////////////////////////////////////////////////////////////////////
  bool CartesianAlgoKS::getNoisyConvergenceMap(boost::filesystem::path& workdir, boost::filesystem::path& InShearMap,
                                               boost::filesystem::path& outConvergenceMap) {
    float SigmaGauss = 0.;
    fs::path datadir = workdir / "data" ;
    fs::path shearMap = datadir / InShearMap;
    if(m_cartesianParam.getSigmaGauss() != 0) {
     SigmaGauss = m_cartesianParam.getSigmaGauss();
    }
    m_cartesianParam.SetSigmaGauss(0.);
    size_t pos;

    if ((InShearMap.string()).empty() == false) {
     if (outConvergenceMap.string().empty() == true) {
      if (false == ((InShearMap.string()).find("NReSample") != std::string::npos)) {
       pos = (InShearMap.string()).find("ShearMap");
       outConvergenceMap = fs::path("EUC_LE3_WL_NoisyConvergenceMapKS_" + (InShearMap.string()).substr (pos));
       perform_MassMapping_Function(shearMap, outConvergenceMap, datadir);
      }
     } else {
       perform_MassMapping_Function(shearMap, outConvergenceMap, workdir);
     }
    }
    m_cartesianParam.SetSigmaGauss(SigmaGauss);
   return true;
  }
////////////////////////////////////////////////////////////////////////////////////////////////////////
  bool CartesianAlgoKS::getDenoisedConvergenceMap(boost::filesystem::path& workdir, boost::filesystem::path& InShearMap,
                                                  boost::filesystem::path& outConvergenceMap) {
    size_t pos;
    fs::path datadir = workdir / "data" ;
    fs::path shearMap = datadir / InShearMap;

    if ((InShearMap.string()).empty() == false) {
     if (outConvergenceMap.string().empty() == true) {
      if (false == ((InShearMap.string()).find("NReSample") != std::string::npos)) {
       pos = (InShearMap.string()).find("ShearMap");
       if (m_cartesianParam.getSigmaGauss() != 0) {
         outConvergenceMap = fs::path("EUC_LE3_WL_DenoisedConvergenceMapKS_" + (InShearMap.string()).substr (pos));
         perform_MassMapping_Function(shearMap, outConvergenceMap, datadir);
       }
      }
     } else {
         perform_MassMapping_Function(shearMap, outConvergenceMap, workdir);
     }
    }
   return true;
  }
////////////////////////////////////////////////////////////////////////////////////////////////////////
  boost::filesystem::path CartesianAlgoKS::getResampledMaps(boost::filesystem::path& workdir,
                  std::vector<fs::path>& filenames) {
    fs::path datadir = workdir / "data" ;
    float SigmaGauss = 0.;
    if(m_cartesianParam.getSigmaGauss() != 0) {
     SigmaGauss = m_cartesianParam.getSigmaGauss();
    }
    m_cartesianParam.SetSigmaGauss(0.);
    boost::filesystem::path outputMaps = fs::path("NReSamplesConvergenceMapList_" + getDateTimeString() + ".json");
    std::ofstream outfile;
    outfile.open ((workdir /outputMaps).string(), std::ios_base::app);
    outfile << "[";
    for (size_t i = 0; i<filenames.size(); i++) {
      size_t pos;
      if ((filenames[i].string()).find("NReSample") != std::string::npos) {
        pos = (filenames[i].string()).find("NReSample"); // position of "NReSample" in input shear map name
        boost::filesystem::path outConvergenceMap {};
        outConvergenceMap = fs::path("EUC_LE3_WL_ConvergenceMapKS_" + (filenames[i].string()).substr (pos) );
        fs::path inFile = datadir / filenames[i];
        //perform_MassMapping_Function(filenames[i], outConvergenceMap, workdir);
        perform_MassMapping_Function(inFile, outConvergenceMap, datadir);
        outfile << outConvergenceMap.filename();
        if (i < filenames.size()-1) {
          outfile << ",";
        }
      }
    }
    outfile << "]";
    outfile.close();
    m_cartesianParam.SetSigmaGauss(SigmaGauss);
   return outputMaps;
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////
  bool CartesianAlgoKS::getSNRMap(boost::filesystem::path& workdir, boost::filesystem::path& NReMaps,
                                  boost::filesystem::path& outConvergenceMap) {
   std::vector<fs::path> filenames = read_filenames(workdir, NReMaps);
   logger.info () <<"number of sampled maps: " << filenames.size();
   std::vector<LE3_2D_MASS_WL_CARTESIAN::ConvergenceMap*> MapList;
   for (size_t i = 0; i<filenames.size(); i++) {
    logger.info () <<"filenames: " << filenames[i];
    LE3_2D_MASS_WL_CARTESIAN::ConvergenceMap *SampledMap = nullptr;
    SampledMap = new ConvergenceMap((workdir/"data"/filenames[i]).native());
    MapList.push_back(SampledMap);
   }
   int xbin = MapList[0]->getXdim();
   int ybin = MapList[0]->getYdim();
   int zbin = MapList[0]->getZdim();
   double* squr = new double[xbin*ybin*zbin];
   std::fill_n(squr, xbin*ybin*zbin, 0);
   for (size_t i = 0; i< MapList.size(); i++) {
    for (int z = 0; z < zbin; z++) {
     for (int y = 0; y < ybin; y++) {
      for (int x = 0; x < xbin; x++) {
       squr[xbin*ybin*z+xbin*y+x] += 
                  pow(MapList[i]->getBinValue(x, y, z), 2);
      }
     }
    }
   }
   double* rms = new double[xbin*ybin*zbin];
   std::fill_n(rms, xbin*ybin*zbin, 0);
   for (int i = 0; i< xbin*ybin*zbin; i++) {
     double mean = squr[i] / double (MapList.size());
     rms[i] = sqrt (mean);
   }
   delete [] squr;
   LE3_2D_MASS_WL_CARTESIAN::ConvergenceMap *myConvergenceMap = new ConvergenceMap(rms, xbin, ybin, zbin);
   if ((true == outConvergenceMap.string().empty())) {
    outConvergenceMap = fs::path("EUC_LE3_WL_SNRConvergenceMap_" + getDateTimeString() + ".fits");
   }
   if (myConvergenceMap != nullptr) {
    std::string name="KAPPA_PATCH";
    m_cartesianParam.setExtName(name);
    myConvergenceMap->writeMap((workdir/"data"/outConvergenceMap).native(), m_cartesianParam);
   }
   delete [] rms;
   delete myConvergenceMap;

   return true;
  }
////////////////////////////////////////////////////////////////////////////////////////////////////////
 // write xml files
////////////////////////////////////////////////////////////////////////////////////////////////////////
  bool CartesianAlgoKS::writeXMLfile (boost::filesystem::path& outConvergenceMapJson,
                                      boost::filesystem::path& outXMLConvergenceMap) {
    // Create a DMOutput object
    DmOutput dm;

    std::string file = (outXMLConvergenceMap).string() ;
    int NResamples = m_cartesianParam.getNSamples();

    auto filenameJson = outConvergenceMapJson.filename();
    auto workdir = outConvergenceMapJson.parent_path();

    std::vector<fs::path> filenames = read_filenames(workdir, filenameJson);
    const fs::path outfile {file};

    if ((m_cartesianParam.getParaFileType()).compare("Conv_Patch") == 0) {
      const std::string product_type = "DpdTwoDMassConvergencePatch";
      auto product = initProduct<dpd::le3::wl::twodmass::out::convergencepatch::dpdTwoDMassConvergencePatch,
    							   pro::le3::wl::twodmass::twoDMassCollectConvergencePatch,
			                       int>(product_type, outfile, NResamples);
      for (size_t i=0; i<filenames.size(); i++) {
        fs::path mapOut (filenames[i].native());
        if ((mapOut.string()).find("NoisyConvergence") != std::string::npos) {
          logger.info() << "writing file: " << mapOut;
          dm.createNoisedPatchXml(product, mapOut);
        }
        if ((mapOut.string()).find("DenoisedConvergence") != std::string::npos) {
          dm.createDenoisedPatchXml(product, mapOut);
        }
        if ((mapOut.string()).find("SNRConvergence") != std::string::npos) {
          dm.createSNRPatchOutputXml(product, mapOut);
        }
      }
      writeProduct<dpd::le3::wl::twodmass::out::convergencepatch::dpdTwoDMassConvergencePatch> (product, outfile);
    }

  /*  if ((m_cartesianParam.getParaFileType()).compare("Conv_Cluster") == 0) {
      const std::string product_type = "DpdTwoDMassConvergenceSingleCluster";
       dm.createSingleClusterXml(outfile, mapOut);
    }*/
    logger.info() << "DM output products created in: " << file;
   return true;
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////
 // END
////////////////////////////////////////////////////////////////////////////////////////////////////////
   } // namespace CartesianKS
  } // namespace TwoDMass
 } // namespace WeakLensing
} // namespace Euclid
