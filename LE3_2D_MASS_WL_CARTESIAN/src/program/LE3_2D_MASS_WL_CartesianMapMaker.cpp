/**
 * @file src/program/LE3_2D_MASS_WL_CartesianMapMaker.cpp
 * @date 10/22/19
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

#include <map>
#include <string>
#include <chrono>
#include <ctime>

// Datamodel for INPUT products
#include "ST_DataModelBindings/dpd/le3/wl/twodmass/inp/euc-test-le3-wl-twodmass-ParamsConvergencePatch.h"

#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/program_options.hpp>
#include "ElementsKernel/ProgramHeaders.h"

#include "LE3_2D_MASS_WL_CARTESIAN/CartesianAlgoKS.h"
#include "LE3_2D_MASS_WL_UTILITIES/DmInput.h"
#include "LE3_2D_MASS_WL_CARTESIAN/CoordinateBound.h"
#include "LE3_2D_MASS_WL_UTILITIES/ReadCatalog.h"
#include "LE3_2D_MASS_WL_UTILITIES/NoisyCatalogData.h"
#include "LE3_2D_MASS_WL_UTILITIES/CatalogData.h"
#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"

namespace po = boost::program_options;

using boost::program_options::options_description;
using boost::program_options::variable_value;
using std::chrono::duration;
using std::chrono::high_resolution_clock;
using std::string;
using std::vector;

using namespace Euclid::WeakLensing::TwoDMass;
using namespace Euclid::WeakLensing::TwoDMass::CartesianKS;
using LE3_2D_MASS_WL_CARTESIAN::CartesianParam;
using LE3_2D_MASS_WL_CARTESIAN::CoordinateBound;

namespace fs = boost::filesystem;

// Input namespace and classes
using namespace dpd::le3::wl::twodmass::inp::paramsconvergencepatch;

static Elements::Logging logger = Elements::Logging::getLogger("LE3_2D_MASS_WL_CartesianMapMaker");

std::vector < std::vector < double> > readClusterCatalog(fs::path& workdir, fs::path& filename) {
   fs::path datadir {workdir / "data"};
  // Variable to save input catalog name
   std::string inputCatalog;
   std::vector<std::vector<double> > data;
   CatalogData cdata("CLUSTER");

  // check whether input cluster catalog file is in fits format(standalone) / XML(pipeline) and then fetch file name
   if (true == checkFileType(datadir /filename, Euclid::WeakLensing::TwoDMass::signFITS)) {
     logger.info("Input Cluster Catalog is in Fits format..");
     inputCatalog = (datadir /filename).string();
     cdata.readCatalog(inputCatalog, data);
     logger.info("Done reading Input Cluster Catalog");
   } else { //Case 2: When input catalog file is in XML format
     logger.info("Input Cluster Catalog is in XML format..");
     cdata.getCatalogData(workdir, filename, data);
     logger.info("Done reading Input Cluster Catalog");
//     std::copy(tempData.begin(), tempData.end(), back_inserter(data));
   }
 return data;
}

class LE3_2D_MASS_WL_CartesianMapMaker : public Elements::Program {

public:

  options_description defineSpecificProgramOptions() override {
  
    options_description options {};
   // input working directory
   options.add_options()
   ("workdir", po::value<string>()->default_value(""), "Work Directory");

   // input log directory: default is none
   options.add_options()
   ("logdir", po::value<string>()->default_value(""), "logs Directory");

   // input catalog file
   options.add_options()
   ("input_ShearCatalog", po::value<string>()->default_value(""), "input shear Catalog in fits/xml format");

   // input cluster catalog file
   options.add_options()
   ("input_ClusterCatalog", po::value<string>()->default_value(""), "input cluster Catalog in fits/xml format");

   // input parameter file
   options.add_options()
   ("paramFile", po::value<string>()->default_value(""), "Input Parameter File in XML");

   // intermediate product (shear Map)
   options.add_options()
   ("outShearMap", po::value<string>()->default_value(""), "output shear Map in fits format");
   options.add_options()
   ("outShearMapList", po::value<string>()->default_value(""), "List of output shear Maps in text/json format");
   // product (Convergence Map) NOT REQUIRED FOR PIPELINE
   options.add_options()
   ("outConvergenceMap", po::value<string>()->default_value(""), "output convergence Map in fits format");

    return options;
  }

  Elements::ExitCode mainMethod(std::map<std::string, variable_value>& args) override {

    //Start time and intro
    auto start_time = std::chrono::system_clock::now();
    logger.info("#");
    logger.info("# Entering mainMethod()");
    logger.info("#");

////////////////////////////////////////////////////////////////////////////////////////////////////////
  // get out time running ID for filename part
////////////////////////////////////////////////////////////////////////////////////////////////////////
  //std::string timeNamePart = getDateTimeString();
  //logger.info("Run ID: " + getDateTimeString());

  logger.info() << "Using " << omp_get_max_threads() << " thread(s) max.";

////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Get the workdir, the data_dir and manage Input output
////////////////////////////////////////////////////////////////////////////////////////////////////////
   fs::path workdir {args["workdir"].as<string>()};
   fs::path InCatalog {args["input_ShearCatalog"].as<string>()};
   fs::path nameCluster {args["input_ClusterCatalog"].as<string>()};
   fs::path ParameterFile {args["paramFile"].as<string>()};
   fs::path logPath (args["logdir"].as<string>());
   fs::path ShearMap {args["outShearMap"].as<std::string>()};
   fs::path ConvergenceMap {args["outConvergenceMap"].as<std::string>()};
   fs::path output_shearMaps = args["outShearMapList"].as<string>();
   fs::path datadir {workdir / "data"};
/*
////////////////////////////////////////////////////////////////////////////////////////////////////////
  // get the filenames from json file
////////////////////////////////////////////////////////////////////////////////////////////////////////
   ReadCatalog readFilename;
   const auto input_filenames = readFilename.readFilenames((workdir/file).native());
   for(const auto& InCatalog: input_filenames) {
   logger.info() << "\t" << InCatalog.string();
*/

////////////////////////////////////////////////////////////////////////////////////////////////////////
 // Reading input catalog
////////////////////////////////////////////////////////////////////////////////////////////////////////
  std::vector<std::vector<double> > catData;
  std::vector<std::vector<double> > clusterData;
  double zMin, zMax; // get min max from catalog
  std::vector<double> centerX;
  std::vector<double> centerY;
  std::vector<double> zmin;
  std::vector<double> clusterId;
  double zmax;
  double z_halo_max = 1.;
////////////////////////////////////////////////////////////////////////////////////////////////////////
  // get Shear catalog Data
////////////////////////////////////////////////////////////////////////////////////////////////////////
   ReadCatalog read;
   read.readShearCatalog(workdir, InCatalog, catData);

////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Variables need in case of Patches to sphere Parameter file
////////////////////////////////////////////////////////////////////////////////////////////////////////
   double catRamin, catRamax, catDecmin, catDecmax;
   vecMinMax(catData[0], &catRamin, &catRamax);
   vecMinMax(catData[1], &catDecmin, &catDecmax);
////////////////////////////////////////////////////////////////////////////////////////////////////////
  // check parameter file exists
////////////////////////////////////////////////////////////////////////////////////////////////////////
  if (false == fs::is_regular(workdir /ParameterFile)) {
   logger.info()<< "Parameter file not found or not in XML format ....";
   return Elements::ExitCode::NOINPUT;
  }
////////////////////////////////////////////////////////////////////////////////////////////////////////
  // read parameter File
////////////////////////////////////////////////////////////////////////////////////////////////////////
  CartesianParam params;
  readParameterFile ((workdir/ParameterFile), params, catRamin, catRamax, catDecmin, catDecmax);
////////////////////////////////////////////////////////////////////////////////////////////////////////
  // get centers of the patches
////////////////////////////////////////////////////////////////////////////////////////////////////////
 // Case 1: if there is only one input catalog i.e. Shear
  if (true == nameCluster.string().empty()) {
   vecMinMax(catData[5], &zMin, &zMax);
   //auto minmax = std::minmax_element(catData[5].begin(), catData[5].end());
   //zMin = *minmax.first;
   //zMax = *minmax.second;
   centerX = params.getMapCenterX();
   centerY = params.getMapCenterY();
   centerX.shrink_to_fit();
   centerY.shrink_to_fit();
   for (size_t i = 0; i < centerX.size(); i++) {
     zmin.push_back(zMin);
     zmax=zMax;
   }
   zmin.shrink_to_fit();
  }
// case 2: if there are both input catalogs
  else /*if (false == nameCluster.string().empty())*/ {
   clusterData = readClusterCatalog(workdir, nameCluster);
   logger.info() << "cluster data size: " << clusterData[1].size();
   zmin.clear();
   centerX.clear();
   centerY.clear();
   zMax = *max_element(clusterData[3].begin(), clusterData[3].end());
   //logger.info()<< "zmax: "<< zMax;
   double zm = params.getZMax();
   //logger.info()<< "zm: "<< zm;
   if (zm < zMax && zm!= 0) {
      zmax=zm;
   } else {
      zmax=zMax;
   }
//   logger.info() << "cluster data size: " << clusterData[1].size();
//   logger.info() << "Mass threshold: " << params.getMassThreshold();
   for (size_t i = 0; i< clusterData[1].size(); ++i) {
    if ((clusterData[5][i] > params.getMassThreshold()) && (clusterData[3][i] < z_halo_max)) {
     centerX.push_back(clusterData[1][i]);
     centerY.push_back(clusterData[2][i]);
     zmin.push_back(clusterData[3][i] + params.getZMargin());
     clusterId.push_back(clusterData[0][i]);
    }
   }
   logger.info() << "Number of selected clusters: " << clusterId.size();
   centerX.shrink_to_fit();
   centerY.shrink_to_fit();
   zmin.shrink_to_fit();
   clusterId.shrink_to_fit();
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////
  if (true == output_shearMaps.string().empty()) {
      output_shearMaps = fs::path("ShearMapsList_" + getDateTimeString() + ".json");
  }
////////////////////////////////////////////////////////////////////////////////////////////////////////
 // Creating Patch of convergence Map from Input catalog
////////////////////////////////////////////////////////////////////////////////////////////////////////
   std::ofstream outfile;
   outfile.open ((workdir /output_shearMaps).string(), std::ios_base::app);
   outfile << "[";
 for (size_t it=0; it<centerX.size(); it++){
   //logger.info() << "Center: " << centerX[it] << "	" << centerY[it];
   //logger.info() << "Cluster Id: " << std::setprecision (15)<< clusterId[it];
   //logger.info() << "z_min: " << zmin[it] << "	z_max: " << zmax;
  // Object to perform Cartesian KS algorithm
   Euclid::WeakLensing::TwoDMass::CartesianKS::CartesianAlgoKS CartesainAlgo(params);
   double rmin, dmin;
   rmin = params.getRaMin(centerX[it]);
   dmin = params.getDecMin(centerY[it]);
   CoordinateBound m_CB(rmin, params.getRaMax(rmin), dmin, params.getDecMax(dmin), zmin[it], zmax);
   logger.info() << "ra_min: " << m_CB.getRaMin() << "	ra_max: " << m_CB.getRaMax();
   logger.info() << "dec_min: " << m_CB.getDecMin() << "	dec_max: " << m_CB.getDecMax();
   logger.info() << "z_min: " << m_CB.getZMin() << "	z_max: " << m_CB.getZMax();
   if ((ConvergenceMap.string()).empty() == false) {
    logger.info("creating convergence Map");
   // ConvergenceMap = data_dir / fs::path("EUC_LE3_WL_ConvergenceMap_"  + getDateTimeString() + ".fits");
    CartesainAlgo.extractConvergnceMap ((datadir /ConvergenceMap).native(), catData, m_CB);
   }

////////////////////////////////////////////////////////////////////////////////////////////////////////
 // Creating Patch of Shear Map from Input catalog
////////////////////////////////////////////////////////////////////////////////////////////////////////
   logger.info("creating shear Map");
   if ((true == ShearMap.string().empty())) {
    ShearMap = fs::path("ShearMap_0" + std::to_string(it) + "_" + getDateTimeString() + ".fits");
   }
   CartesainAlgo.extractShearMap ((datadir /ShearMap).native(), catData, m_CB);
   outfile << ShearMap.filename();
   if (it < centerX.size()-1) {
     outfile << ",";
   }
   ShearMap.clear();
  // Get associated SNR maps
  logger.info("creating associated snr shear Maps");
   if(params.getNSamples() > 0) {
    outfile << ",";
    for (int iter = 0; iter < params.getNSamples(); ++iter) {
      std::vector<std::vector<double> > RanData;
      NoisyCatalogData randomise;
      RanData = randomise.create_noisy_data(catData);
      ShearMap = fs::path("ShearMap_0" + std::to_string(it) + "_NReSample_0" + std::to_string(iter) + "_" +
                                         getDateTimeString() + ".fits");
      CartesainAlgo.extractShearMap ((datadir / ShearMap).native(), RanData, m_CB);
////////////////////////////////////////////////////////////////////////////////////////////////////////
 // Generate the json/txt product
////////////////////////////////////////////////////////////////////////////////////////////////////////
      outfile << ShearMap.filename();
      if (iter < params.getNSamples()-1) {
        outfile << ",";
      }
      ShearMap.clear();
    } //end of if
   }// end of for loop for Nsample
 } //end of for loop
   outfile << "]";
   outfile.close();
////////////////////////////////////////////////////////////////////////////////////////////////////////
 //End of patch extraction
////////////////////////////////////////////////////////////////////////////////////////////////////////
   auto end_time = std::chrono::system_clock::now();
   std::chrono::duration<double> elapsed_seconds = end_time-start_time;
   std::time_t endTime = std::chrono::system_clock::to_time_t(end_time);
   logger.info() << "finished Program at " << std::ctime(&endTime);
   logger.info() << "elapsed time: " << elapsed_seconds.count() << "s";

    logger.info("#");
    logger.info("# Exiting mainMethod()");
    logger.info("#");

    return Elements::ExitCode::OK;
  }

};

MAIN_FOR(LE3_2D_MASS_WL_CartesianMapMaker)
