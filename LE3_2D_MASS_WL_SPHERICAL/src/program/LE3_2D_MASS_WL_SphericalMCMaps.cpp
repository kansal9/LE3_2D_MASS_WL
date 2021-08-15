/**
 * @file src/program/LE3_2D_MASS_WL_SphericalMCMaps.cpp
 * @date 01/10/20
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

#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/program_options.hpp>
#include "ElementsKernel/ProgramHeaders.h"

#include "LE3_2D_MASS_WL_UTILITIES/DmInput.h"
#include "LE3_2D_MASS_WL_UTILITIES/CatalogData.h"
#include "LE3_2D_MASS_WL_SPHERICAL/GetSphericalMCMaps.h"
#include "LE3_2D_MASS_WL_SPHERICAL/Sph_mass_mapping.h"
#include "LE3_2D_MASS_WL_SPHERICAL/SphericalUtils.h"

// inDatamodel for INPUT products
#include "ST_DataModelBindings/dpd/le3/wl/twodmass/inp/euc-test-le3-wl-twodmass-ParamsConvergenceSphere.h"

namespace po = boost::program_options;
namespace fs = boost::filesystem;
using boost::program_options::options_description;
using boost::program_options::variable_value;
using std::chrono::duration;
using std::chrono::high_resolution_clock;
using std::string;
using std::vector;

// Input namespace and classes
using namespace dpd::le3::wl::twodmass::inp::paramsconvergencesphere;

using namespace Euclid::WeakLensing::TwoDMass;
using namespace Euclid::WeakLensing::TwoDMass::Spherical;

using LE3_2D_MASS_WL_SPHERICAL::Sph_map_maker;
using LE3_2D_MASS_WL_SPHERICAL::Sph_mass_mapping;
using LE3_2D_MASS_WL_SPHERICAL::SphericalParam;
using Euclid::WeakLensing::TwoDMass::Spherical::SphericalIO;
using LE3_2D_MASS_WL_SPHERICAL::GetSphericalMCMaps;

static Elements::Logging logger = Elements::Logging::getLogger("LE3_2D_MASS_WL_SphericalMCMaps");

class LE3_2D_MASS_WL_SphericalMCMaps : public Elements::Program {

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
   ("inCatalogFile", po::value<string>()->default_value(""), "input Catalog in fits/xml format with path");

   // input parameter file
   options.add_options()
   ("sphericalParameterFile", po::value<string>()->default_value(""), "second Input Parameter File");

   // output product (list of MC Maps)
   options.add_options()
   ("MCConvergenceMaps", po::value<string>()->default_value(""), "MC convergence Maps Name in txt/jason file");

    return options;
  }

  Elements::ExitCode mainMethod(std::map<std::string, variable_value>& args) override {

    Elements::Logging logger = Elements::Logging::getLogger("LE3_2D_MASS_WL_SphericalMCMaps");

    //Start time and intro
    auto start_SphMCMainTime = std::chrono::system_clock::now();

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
   fs::path inCatalog {args["inCatalogFile"].as<string>()};
   fs::path inSphParamFile {args["sphericalParameterFile"].as<std::string>()};
   fs::path logPath (args["logdir"].as<string>());
   std::vector<fs::path> MCConvergenceMaps;// = args["MCConvergenceMaps"].as<std::vector<std::string> >();
   std::vector<fs::path> MCShearMaps;
////////////////////////////////////////////////////////////////////////////////////////////////////////
  // check parameter file exists
////////////////////////////////////////////////////////////////////////////////////////////////////////
  if (false == fs::is_regular(workdir/inSphParamFile)) {
   logger.info()<< "Parameter file not found or not in XML format ....";
   return Elements::ExitCode::NOINPUT;
  }
////////////////////////////////////////////////////////////////////////////////////////////////////////
  // check whether parameter file is of XML format or not
////////////////////////////////////////////////////////////////////////////////////////////////////////
   SphericalParam SphParam;
   readSphericalParameterFile((workdir/inSphParamFile), SphParam);
////////////////////////////////////////////////////////////////////////////////////////////////////////
 // Object to Write Spherical Map
////////////////////////////////////////////////////////////////////////////////////////////////////////
   SphericalIO SphericalIO(SphParam);
////////////////////////////////////////////////////////////////////////////////////////////////////////
 // Object to perform Spherical Mass Mapping
////////////////////////////////////////////////////////////////////////////////////////////////////////
    Sph_mass_mapping mapping;
////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Variable to save input catalog name
////////////////////////////////////////////////////////////////////////////////////////////////////////
   std::string inCatalogFile;
////////////////////////////////////////////////////////////////////////////////////////////////////////
 // Reading input catalog
////////////////////////////////////////////////////////////////////////////////////////////////////////
   CatalogData read;
   logger.info("# Reading input catalog");
   std::vector<std::vector<double> > inData;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // check whether input catalog file is in fits format(standalone) / XML(pipeline) and then fetch Catalogfile name
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Case 1: When input catalog file is in FITS format
   if (true == checkFileType(datadir /inCatalog, Euclid::WeakLensing::TwoDMass::signFITS)) {
     logger.info("Input Catalog is in Fits format..");
     inCatalogFile = (datadir /inCatalog).string();
     //inData = read.readCatalog(inCatalogFile);
     read.readCatalog(inCatalogFile, inData);
     logger.info("Done reading Input shear Catalog");
   }  else { // Case 2: When input catalog file is in XML format
     //inData = read.getCatalogData(workdir, inCatalog);
     read.getCatalogData(workdir, inCatalog, inData);
   }
////////////////////////////////////////////////////////////////////////////////////////////////////////
 // get Denoised Shear Map
////////////////////////////////////////////////////////////////////////////////////////////////////////
  GetSphericalMCMaps MCMaps(inData, SphParam);
  std::pair<Healpix_Map<double>, Healpix_Map<double> > DeNoisedShearPair;
  DeNoisedShearPair = MCMaps.getDeNoisedShearMap();

   /*fs::path gamma = "EUC_LE3_WL_DeNoisedGamma_" + getDateTimeString() + ".fits";
   SphericalIO.write_Map((workdir/gamma).native(), DeNoisedShearPair.first, "GAMMA1");
   SphericalIO.write_Map((workdir/gamma).native(), DeNoisedShearPair.second, "GAMMA2");*/

////////////////////////////////////////////////////////////////////////////////////////////////////////
 // get noised N Shear Maps
////////////////////////////////////////////////////////////////////////////////////////////////////////

  std::vector<std::pair<Healpix_Map<double>, Healpix_Map<double> > > NoisedShearMapList;
  NoisedShearMapList = MCMaps.getNoisedShearMaps();

  for (int i = 0; i<NoisedShearMapList.size(); i++) {
   fs::path shearMap =  (datadir / fs::path ("EUC_LE3_WL_MC_ShearMap_0" +
                                 std::to_string(i)+ "_" + getDateTimeString() + ".fits"));
   MCShearMaps.push_back(shearMap);
  }

  for (const auto &iter: NoisedShearMapList){
//  for (int iter = 0; iter<NoisedShearMapList.size(); iter++) {
    auto index = (&iter - &NoisedShearMapList[0]);
    logger.info() << index;
    MCMaps.performAddition(DeNoisedShearPair, NoisedShearMapList[index], MCShearMaps[index].native() );
  }
////////////////////////////////////////////////////////////////////////////////////////////////////////
 // perform mass mapping on N shear Maps
////////////////////////////////////////////////////////////////////////////////////////////////////////
  logger.info("Running Spherical KS Mass Mapping");
  for (int i = 0; i<MCShearMaps.size(); i++) {
    fs::path Kappa =  datadir / fs::path("EUC_LE3_WL_MC_Kappa_0" + std::to_string(i)+ "_" +
                                  getDateTimeString() + ".fits");
    Healpix_Map<double> mapShearE;
    read_Healpix_map_from_fits(MCShearMaps[i].native(), mapShearE);

    std::pair<Healpix_Map<double>, Healpix_Map<double> > kappaPair =
                      mapping.create_SheartoConvMap(mapShearE, DeNoisedShearPair.second);

   SphericalIO.write_Map(Kappa.native(), kappaPair.first, "KappaE");
   SphericalIO.write_Map(Kappa.native(), kappaPair.second, "KappaB");

   MCConvergenceMaps.push_back(Kappa);
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////
 // Generate the jason/txt product
////////////////////////////////////////////////////////////////////////////////////////////////////////
  fs::path outFilenames = args["MCConvergenceMaps"].as<string>();
  std::ofstream outSphConvfile ((workdir /outFilenames).string());
  outSphConvfile << "[";
  size_t i = 0;
  while (i != MCConvergenceMaps.size()) {
    fs::path outConvMapName(MCConvergenceMaps[i]);
    outSphConvfile << outConvMapName.filename();
    i++;
    if (i<MCConvergenceMaps.size()) {
      outSphConvfile << ", ";
    }
  }
  outSphConvfile << "]";

////////////////////////////////////////////////////////////////////////////////////////////////////////
 //End of MC Main
////////////////////////////////////////////////////////////////////////////////////////////////////////
    auto end_SphMCMainTime = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end_SphMCMainTime-start_SphMCMainTime;
    //std::time_t end_time = std::chrono::system_clock::to_time_t(end_SphMCMainTime);
    //logger.info() << "finished Program at " << std::ctime(&end_time);
    logger.info() << "elapsed time: " << elapsed_seconds.count() << "s";

    logger.info("#");
    logger.info("# Exiting mainMethod()");
    logger.info("#");

    return Elements::ExitCode::OK;
  }

};

MAIN_FOR(LE3_2D_MASS_WL_SphericalMCMaps)
