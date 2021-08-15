/**
 * @file src/program/LE3_2D_MASS_WL_SphericalMapMaker.cpp
 * @date 11/30/19
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

#include <boost/property_tree/ptree.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include "ElementsKernel/ProgramHeaders.h"

#include "LE3_2D_MASS_WL_SPHERICAL/Sph_map_maker.h"
#include "LE3_2D_MASS_WL_UTILITIES/ReadCatalog.h"
#include "LE3_2D_MASS_WL_SPHERICAL/SphericalParam.h"
#include "LE3_2D_MASS_WL_SPHERICAL/SphericalIO.h"
#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"
#include "LE3_2D_MASS_WL_SPHERICAL/SphericalUtils.h"

#include "LE3_2D_MASS_WL_UTILITIES/DmOutput.h"
#include "LE3_2D_MASS_WL_UTILITIES/DmInput.h"

// Datamodel for INPUT products
#include "ST_DataModelBindings/dpd/le3/wl/twodmass/inp/euc-test-le3-wl-twodmass-ParamsConvergenceSphere.h"

using std::chrono::duration;
using std::chrono::high_resolution_clock;
using std::string;
using std::vector;

using LE3_2D_MASS_WL_SPHERICAL::Sph_map_maker;
using LE3_2D_MASS_WL_SPHERICAL::SphericalParam;
using Euclid::WeakLensing::TwoDMass::Spherical::SphericalIO;

// Input namespace and classes
using namespace dpd::le3::wl::twodmass::inp::paramsconvergencesphere;

namespace po = boost::program_options;
namespace fs = boost::filesystem;
using boost::program_options::options_description;
using boost::program_options::variable_value;
using namespace Euclid::WeakLensing::TwoDMass;
using namespace Euclid::WeakLensing::TwoDMass::Spherical;

static Elements::Logging logger = Elements::Logging::getLogger("LE3_2D_MASS_WL_SphericalMapMaker");

class LE3_2D_MASS_WL_SphericalMapMaker : public Elements::Program {

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
   ("inputShearCatalog", po::value<string>()->default_value(""), "input Catalog in fits/xml format");

   // input parameter file
   options.add_options()
   ("sphericalParameterFile", po::value<string>()->default_value(""), "Input Parameter File in XML");

   // output product (shear Map)
   options.add_options()
   ("outShearMap", po::value<string>()->default_value(""), "output shear Map E & B mode in fits format");

   //options.add_options()
   //("outShearMapList", po::value<string>()->default_value(""), "output shear Map list in json format");
   options.add_options()
   ("GalCountMap", po::value<string>()->default_value(""),
    "output Galaxy count Map which conatins number of Galaxies per pixel for each redshift bin in json format");

    return options;
  }

  Elements::ExitCode mainMethod(std::map<std::string, variable_value>& args) override {

    //Elements::Logging logger = Elements::Logging::getLogger("LE3_2D_MASS_WL_SphericalMapMaker");

   // start time & intro
   auto SphMapMaker_start = std::chrono::system_clock::now();

    logger.info("#");
    logger.info("# Entering mainMethod()");
    logger.info("#");
////////////////////////////////////////////////////////////////////////////////////////////////////////
  // get out time running ID for filename part
////////////////////////////////////////////////////////////////////////////////////////////////////////
   logger.info() << "Using " << omp_get_max_threads() << " thread(s) max.";
////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Get the workdir and manage Input output
////////////////////////////////////////////////////////////////////////////////////////////////////////
   fs::path workdir {args["workdir"].as<string>()};
   fs::path datadir {workdir / "data"};
   fs::path fileParam {args["sphericalParameterFile"].as<std::string>()};
   fs::path InputCatalog {args["inputShearCatalog"].as<string>()};
   fs::path gamma {args["outShearMap"].as<string>()};
   fs::path GalCountFilename {};

////////////////////////////////////////////////////////////////////////////////////////////////////////
 // Reading input catalog
////////////////////////////////////////////////////////////////////////////////////////////////////////
   std::vector<std::vector<double> > Data;
   ReadCatalog read;
   read.readShearCatalog(workdir, InputCatalog, Data);
////////////////////////////////////////////////////////////////////////////////////////////////////////
  // check parameter file exists
////////////////////////////////////////////////////////////////////////////////////////////////////////
  if (false == fs::is_regular(workdir /fileParam)) {
   logger.info()<< "Parameter file not found or not in XML format ....";
   return Elements::ExitCode::NOINPUT;
  }
////////////////////////////////////////////////////////////////////////////////////////////////////////
  // check whether parameter file is of XML format or .ini format
////////////////////////////////////////////////////////////////////////////////////////////////////////
   SphericalParam SphParam;
   readSphericalParameterFile((workdir/fileParam), SphParam);
////////////////////////////////////////////////////////////////////////////////////////////////////////
 // Object to Write Spherical Map
////////////////////////////////////////////////////////////////////////////////////////////////////////
   SphericalIO SphericalIO(SphParam);
////////////////////////////////////////////////////////////////////////////////////////////////////////
 // Creating healpix Shear Map from input catalog
////////////////////////////////////////////////////////////////////////////////////////////////////////
   logger.info("# Running Map Maker: Creating Healpix Gamma Maps");
   Sph_map_maker mapMaker(SphParam);
   //std::pair<Healpix_Map<double>, Healpix_Map<double> > shearPair = mapMaker.create_ShearMap(Data);
   auto [ Shear1, Shear2, GalCount ] = mapMaker.create_ShearMap(Data);

////////////////////////////////////////////////////////////////////////////////////////////////////////
 // set FITS filenames
////////////////////////////////////////////////////////////////////////////////////////////////////////
  if ((gamma.string()).empty() == true) {
    gamma = fs::path("EUC_LE3_WL_Gamma_NSide" + std::to_string(SphParam.getNside()) + "_" +
           getDateTimeString() + ".fits");
   }
  fs::path GalCountMap = args["GalCountMap"].as<string>();
  if ((GalCountMap.string()).empty() == true) {
    GalCountMap = fs::path("EUC_LE3_WL_GalCount_NSide" + std::to_string(SphParam.getNside()) + "_" +
           getDateTimeString() + ".json");
   }

  GalCountFilename = fs::path("EUC_LE3_WL_GalCount_NSide" + std::to_string(SphParam.getNside()) + "_" +
           getDateTimeString() + ".fits");

  Data.clear();
////////////////////////////////////////////////////////////////////////////////////////////////////////
 // Writing Fits files
////////////////////////////////////////////////////////////////////////////////////////////////////////
   logger.info("# Writing healpix Gamma Maps");

   SphericalIO.write_Map((workdir/gamma).native(), Shear1, "GAMMA1");
   SphericalIO.write_Map((workdir/gamma).native(), Shear2, "GAMMA2");
   SphericalIO.write_Map((datadir/GalCountFilename).native(), GalCount, "GALCOUNT");

    std::ofstream outfile;
    outfile.open ((workdir / GalCountMap).string(), std::ios_base::app);
    outfile << "[";
    outfile << GalCountFilename.filename();
    outfile << "]";
    outfile.close();

   logger.info("Done!");
////////////////////////////////////////////////////////////////////////////////////////////////////////
 // End of the Program and Gives the final time of execution
////////////////////////////////////////////////////////////////////////////////////////////////////////
  auto SphMapMaker_end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = SphMapMaker_end-SphMapMaker_start;
  std::time_t SphMapMaker_end_time = std::chrono::system_clock::to_time_t(SphMapMaker_end);
  logger.info() << "finished Program at " << std::ctime(&SphMapMaker_end_time);
  logger.info() << "elapsed time: " << elapsed_seconds.count() << "s";

    logger.info("# Exiting mainMethod()");
    logger.info("#");

    return Elements::ExitCode::OK;
  }

};

MAIN_FOR(LE3_2D_MASS_WL_SphericalMapMaker)
