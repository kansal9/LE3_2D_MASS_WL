/**
 * @file src/program/LE3_2D_MASS_WL_CartesianMCMaps.cpp
 * @date 10/25/19
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

#include "LE3_2D_MASS_WL_CARTESIAN/GetCartesianMCMaps.h"
#include "LE3_2D_MASS_WL_CARTESIAN/CartesianAlgoKS.h"
#include "LE3_2D_MASS_WL_CARTESIAN/CoordinateBound.h"
#include "LE3_2D_MASS_WL_UTILITIES/ReadCatalog.h"
#include "LE3_2D_MASS_WL_UTILITIES/CatalogData.h"
#include "LE3_2D_MASS_WL_UTILITIES/DmInput.h"
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
using LE3_2D_MASS_WL_CARTESIAN::GetCartesianMCMaps;

namespace fs = boost::filesystem;

// Input namespace and classes
using namespace dpd::le3::wl::twodmass::inp::paramsconvergencepatch;

static Elements::Logging logger = Elements::Logging::getLogger("LE3_2D_MASS_WL_CartesianMCMaps");

class LE3_2D_MASS_WL_CartesianMCMaps : public Elements::Program {

public:

  options_description defineSpecificProgramOptions() override {
  
    options_description options {};
    
   // input working directory
   options.add_options()
   ("workdir", po::value<string>()->default_value(""), "Working Directory Path");

   // input log directory
   options.add_options()
   ("logdir", po::value<string>()->default_value(""), "log Directory Path");

   // input shear catalog file
   options.add_options()
   ("inputShearCatalog", po::value<string>()->default_value(""), "input shear Catalog");

   // input parameter file
   options.add_options()
   ("PatchparameterFile", po::value<string>()->default_value(""), " Input Parameter File");

   // output product (list of MC Maps)
   options.add_options()
   ("MCConvergenceMaps", po::value<string>()->default_value(""), "MC convergence Maps Name in txt/jason file");

    return options;
  }

  Elements::ExitCode mainMethod(std::map<std::string, variable_value>& args) override {

    //Elements::Logging logger = Elements::Logging::getLogger("LE3_2D_MASS_WL_CartesianMCMaps");
    //Start time and intro
    auto start_MCMainTime = std::chrono::system_clock::now();

    logger.info("#");
    logger.info("# Entering mainMethod()");
    logger.info("#");

////////////////////////////////////////////////////////////////////////////////////////////////////////
  // get out time running ID for filename part
////////////////////////////////////////////////////////////////////////////////////////////////////////
  //logger.info("Run ID: " + getDateTimeString());

  logger.info() << "Using " << omp_get_max_threads() << " thread(s) max.";

////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Get the workdir, the data_dir and manage Input output
////////////////////////////////////////////////////////////////////////////////////////////////////////
   fs::path workdir {args["workdir"].as<string>()};
   fs::path datadir {workdir / "data"};
   fs::path inFilename {args["inputShearCatalog"].as<string>()};
   fs::path inPatchParamFile {args["PatchparameterFile"].as<std::string>()};
   fs::path logPath (args["logdir"].as<string>());
   std::vector<fs::path> MCConvergenceMaps;// = args["MCConvergenceMaps"].as<std::vector<std::string> >();
   std::vector<fs::path> MCShearMaps;
////////////////////////////////////////////////////////////////////////////////////////////////////////
  // check parameter file exists
////////////////////////////////////////////////////////////////////////////////////////////////////////
  if (false == fs::is_regular(workdir/inPatchParamFile)) {
   logger.info()<< "Parameter file not found or not in XML format ....";
   return Elements::ExitCode::NOINPUT;
  }
////////////////////////////////////////////////////////////////////////////////////////////////////////
  // check whether parameter file is of XML format or not
////////////////////////////////////////////////////////////////////////////////////////////////////////
  CartesianParam Patchparams;
  readParameterFile ((workdir/inPatchParamFile), Patchparams);

////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Variable to save input catalog name
////////////////////////////////////////////////////////////////////////////////////////////////////////
   std::string inCatalogFile;
////////////////////////////////////////////////////////////////////////////////////////////////////////
 // Reading input catalog
////////////////////////////////////////////////////////////////////////////////////////////////////////
  ReadCatalog read;
  logger.info("Reading Input Catalog");
  std::vector<std::vector<double> > catData;
  double zMin, zMax; // get min max from catalog
  std::vector<double> centerX;
  std::vector<double> centerY;
  std::vector<double> zmin;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // check whether input catalog file is in fits format(standalone) / XML(pipeline) and then fetch Catalogfile name
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 if (false == inFilename.string().empty()) {
   read.readShearCatalog(workdir, inFilename, catData);
   vecMinMax(catData[5], &zMin, &zMax);
   centerX = Patchparams.getMapCenterX();
   centerY = Patchparams.getMapCenterY();
   for (size_t i = 0; i < centerX.size(); i++) {
     zmin.push_back(zMin);
   }
 } else {
   logger.info() << "ERROR: There is no Input Shear Catalogue . . . .";
 }
////////////////////////////////////////////////////////////////////////////////////////////////////////
 // get Denoised Shear Map
////////////////////////////////////////////////////////////////////////////////////////////////////////
for (size_t it=0; it<centerX.size(); it++){
  double ramin, decmin;
  ramin = Patchparams.getRaMin(centerX[it]);
  decmin = Patchparams.getDecMin(centerY[it]);
  CoordinateBound m_CB(ramin, Patchparams.getRaMax(ramin), decmin,
                     Patchparams.getDecMax(decmin), zmin[it], zMax);

  GetCartesianMCMaps MCMaps(catData, Patchparams, m_CB);
  LE3_2D_MASS_WL_CARTESIAN::ShearMap *DeNoisedShearMap;
  DeNoisedShearMap = MCMaps.getDeNoisedShearMap();
  //DeNoisedShearMap->writeMap("denoisedShearMap.fits", params);

////////////////////////////////////////////////////////////////////////////////////////////////////////
 // get noised N Shear Maps
////////////////////////////////////////////////////////////////////////////////////////////////////////

  std::vector<LE3_2D_MASS_WL_CARTESIAN::ShearMap*> NoisedShearMapList;
  NoisedShearMapList = MCMaps.getNoisedShearMaps();

  for (size_t i = 0; i<NoisedShearMapList.size(); i++) {
   fs::path shearMapName =  (datadir / fs::path ("EUC_LE3_WL_ShearMap_0" +
                                 std::to_string(i)+ "_" + getDateTimeString() + ".fits"));
   MCShearMaps.push_back(shearMapName);
  }

  for (const auto &iter: NoisedShearMapList){
    //logger.info() << (*iter );
    logger.info() << (&iter - &NoisedShearMapList[0]);
//    logger.info() << "ngal: " <<(NoisedShearMapList[0]->getNumberOfGalaxies());
    auto index = (&iter - &NoisedShearMapList[0]);
    MCMaps.performAddition(*DeNoisedShearMap, *iter, MCShearMaps[index].native() );
  }

 delete DeNoisedShearMap;
 DeNoisedShearMap = nullptr;
////////////////////////////////////////////////////////////////////////////////////////////////////////
 // perform mass mapping on N shear Maps
////////////////////////////////////////////////////////////////////////////////////////////////////////
  logger.info("entering Mass mapping()");
  for (size_t i = 0; i<MCShearMaps.size(); i++) {
   fs::path ConvergenceMapName =  (datadir / fs::path ("EUC_LE3_WL_ConvergenceMap_0" +
                                 std::to_string(i)+ "_" + getDateTimeString() + ".fits"));

   LE3_2D_MASS_WL_CARTESIAN::ConvergenceMap *m_ConvergenceMap;
   LE3_2D_MASS_WL_CARTESIAN::MassMapping mass(Patchparams);
   m_ConvergenceMap = mass.getSheartoConv(MCShearMaps[i].native());
   m_ConvergenceMap->writeMap(ConvergenceMapName.native(), Patchparams);

   MCConvergenceMaps.push_back(ConvergenceMapName);
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////
 // Generate the jason/txt product
////////////////////////////////////////////////////////////////////////////////////////////////////////
  fs::path outConv_filenames = args["MCConvergenceMaps"].as<string>();
  std::ofstream outConvfile ((workdir /fs::path (std::to_string(it)) / outConv_filenames).string());
  outConvfile << "[";
  size_t i = 0;
  while (i != MCConvergenceMaps.size()) {
    fs::path outConvMapName(MCConvergenceMaps[i]);
    outConvfile << outConvMapName.filename();
    i++;
    if (i<MCConvergenceMaps.size()) {
      outConvfile << ", ";
    }
  }
  outConvfile << "]";
}
////////////////////////////////////////////////////////////////////////////////////////////////////////
 //End of MC Main
////////////////////////////////////////////////////////////////////////////////////////////////////////
    auto end_MCMainTime = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end_MCMainTime-start_MCMainTime;
    //std::time_t end_time = std::chrono::system_clock::to_time_t(end_MCMainTime);
    //logger.info() << "finished Program at " << std::ctime(&end_time);
    logger.info() << "elapsed time: " << elapsed_seconds.count() << "s";

    logger.info("#");
    logger.info("# Exiting mainMethod()");
    logger.info("#");

    return Elements::ExitCode::OK;
  }

};

MAIN_FOR(LE3_2D_MASS_WL_CartesianMCMaps)
