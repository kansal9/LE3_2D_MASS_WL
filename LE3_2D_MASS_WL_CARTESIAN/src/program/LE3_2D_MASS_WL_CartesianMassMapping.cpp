/**
 * @file src/program/LE3_2D_MASS_WL_CartesianMassMapping.cpp
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

// Datamodel for OUTPUT products
#include "ST_DataModelBindings/dpd/le3/wl/twodmass/out/euc-test-le3-wl-twodmass-ConvergencePatch.h"

// Datamodel for INPUT products
#include "ST_DataModelBindings/dpd/le3/wl/twodmass/inp/euc-test-le3-wl-twodmass-ParamsConvergencePatch.h"

#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/program_options.hpp>
#include "ElementsKernel/ProgramHeaders.h"
#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"
#include "LE3_2D_MASS_WL_CARTESIAN/CartesianAlgoKS.h"
#include "LE3_2D_MASS_WL_UTILITIES/DmInput.h"
#include "LE3_2D_MASS_WL_UTILITIES/DmOutput.h"

using boost::program_options::options_description;
using boost::program_options::variable_value;
using std::chrono::duration;
using std::chrono::high_resolution_clock;
using std::string;
using std::vector;

using namespace Euclid::WeakLensing::TwoDMass;
using namespace Euclid::WeakLensing::TwoDMass::CartesianKS;
using LE3_2D_MASS_WL_CARTESIAN::CartesianParam;

namespace po = boost::program_options;
namespace fs = boost::filesystem;

// Input namespace and classes
using namespace dpd::le3::wl::twodmass::inp::paramsconvergencepatch;

// Output namespace and classes
using namespace dpd::le3::wl::twodmass::out::convergencepatch;

static Elements::Logging logger = Elements::Logging::getLogger("LE3_2D_MASS_WL_CartesianMassMapping");

class LE3_2D_MASS_WL_CartesianMassMapping : public Elements::Program {

public:

  options_description defineSpecificProgramOptions() override {
  
    options_description options {};
   // input working directory
   options.add_options()
   ("workdir", po::value<string>()->default_value(""), "Working Directory");

   // input log directory: default is none
   options.add_options()
   ("logdir", po::value<string>()->default_value(""), "log Directory");

   // input shear map file
   options.add_options()
   ("inputShearMap", po::value<string>()->default_value(""), "input Shear Maps in json/fits format");
   // input convergence map file
   options.add_options()
   ("inputConvMap", po::value<string>()->default_value(""), "input convergence Map in fits format");

   // input parameter file
   options.add_options()
   ("paramFile", po::value<string>()->default_value(""), "Input Parameter File in XML");

   // output products (Convergence Map)
   options.add_options()
   ("outConvMap", po::value<string>()->default_value(""), "output convergence Map in fits format");

  // output products (Convergence Maps)
   options.add_options()
   ("outConvMaps", po::value<string>()->default_value(""), "output list of convergence Maps in json format");

   options.add_options()
   ("outConvMapXML", po::value<string>()->default_value("outConvMapXML.xml"), "output convergence Map in xml format");
   // output products (shear Map)
   options.add_options()
   ("outShearMap", po::value<string>()->default_value(""), "output shear Map in fits format");

    return options;
  }

  Elements::ExitCode mainMethod(std::map<std::string, variable_value>& args) override {

    //Elements::Logging logger = Elements::Logging::getLogger("LE3_2D_MASS_WL_CartesianMassMapping");
    //Start time and intro
    auto start_MappingTime = std::chrono::system_clock::now();
    logger.info("#");
    logger.info("# Entering mainMethod()");
    logger.info("#");

////////////////////////////////////////////////////////////////////////////////////////////////////////
  // get out time running ID for filename part
////////////////////////////////////////////////////////////////////////////////////////////////////////
  //std::string timeNamePart = getDateTimeString();
  //logger.info("Run ID: " + timeNamePart);

  logger.info() << "Using " << omp_get_max_threads() << " thread(s) max.";

////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Get the workdir, the data_dir and manage Input output
////////////////////////////////////////////////////////////////////////////////////////////////////////
   fs::path workdir {args["workdir"].as<string>()};
   fs::path datadir {workdir / "data"};
   fs::path InShearMap {args["inputShearMap"].as<string>()};
   fs::path InConvergenceMap {args["inputConvMap"].as<string>()};
   fs::path ParamFile {args["paramFile"].as<std::string>()};
   fs::path logPath (args["logdir"].as<string>());
   fs::path outShearMap {args["outShearMap"].as<std::string>()};
   fs::path outConvergenceMap {args["outConvMap"].as<std::string>()};
   fs::path outputConvMaps {args["outConvMaps"].as<std::string>()};
   auto out_xml_product_file_name =args["outConvMapXML"].as<std::string>();
////////////////////////////////////////////////////////////////////////////////////////////////////////
  // check parameter file exists
////////////////////////////////////////////////////////////////////////////////////////////////////////
  if (false == fs::is_regular(workdir/ParamFile)) {
   logger.info()<< "Parameter file not found or not in XML format ....";
   return Elements::ExitCode::NOINPUT;
  }
////////////////////////////////////////////////////////////////////////////////////////////////////////
  // check whether parameter file is of XML format or not
////////////////////////////////////////////////////////////////////////////////////////////////////////
  CartesianParam params;
  readParameterFile ((workdir/ParamFile), params);

////////////////////////////////////////////////////////////////////////////////////////////////////////
 // Object to perform Cartesian KS algorithm
////////////////////////////////////////////////////////////////////////////////////////////////////////
  Euclid::WeakLensing::TwoDMass::CartesianKS::CartesianAlgoKS CartesainAlgo(params);

////////////////////////////////////////////////////////////////////////////////////////////////////////
 // perform inverse KS Mass Mapping
////////////////////////////////////////////////////////////////////////////////////////////////////////
    if ((InConvergenceMap.string()).empty() == false) {
    logger.info("KS conversion from Convergence to shear Map");
     if ((outShearMap.string()).empty() == true) {
      outShearMap = fs::path("EUC_LE3_WL_ShearMapKS_" + getDateTimeString() + ".fits");
     }
     CartesainAlgo.performInverseKSMassMapping((datadir/InConvergenceMap).native(), (datadir/outShearMap).native());
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////
 // set FITS output filenames and perform KS Mass Mapping
////////////////////////////////////////////////////////////////////////////////////////////////////////
    if ((outputConvMaps.string()).empty() == true) {
       outputConvMaps = fs::path("ConvergenceMapsList_" + getDateTimeString() + ".json");
    }

    if (params.getNItReducedShear() != 0) {
     CartesainAlgo.performReducedShear(InShearMap, outConvergenceMap, workdir, out_xml_product_file_name);
    } else {
      // If input Shear map in fits format
      if (true == checkFileType((datadir/InShearMap).native(), Euclid::WeakLensing::TwoDMass::signFITS)) {
       if ((InShearMap.string()).empty() == false) {
        if ((InShearMap.string()).find("NReSample") != std::string::npos) {
         std::vector<fs::path> filenames;
         filenames.push_back(InShearMap);
         boost::filesystem::path NReConvMaps = CartesainAlgo.getResampledMaps(workdir, filenames);
        } else {
         CartesainAlgo.getNoisyConvergenceMap(workdir, InShearMap, outConvergenceMap);
         outConvergenceMap.clear();
         if (params.getSigmaGauss() != 0) {
            CartesainAlgo.getDenoisedConvergenceMap(workdir, InShearMap, outConvergenceMap);
         }
        }
       }// end of internal if-else
      } else { // input Shear map in json file
         std::ofstream outfile;
         outfile.open ((workdir / outputConvMaps).string(), std::ios_base::app);
         outfile << "[";
         std::vector<fs::path> filenames = read_filenames(workdir, InShearMap);
         bool NreSamples_Exists = false;
         for (size_t i = 0; i<filenames.size(); i++) {
           CartesainAlgo.getNoisyConvergenceMap(workdir, filenames[i], outConvergenceMap);
           if ((outConvergenceMap.string()).empty() == false) {
             outfile << outConvergenceMap.filename();
           }
           outConvergenceMap.clear();
           if (params.getSigmaGauss() != 0) {
             CartesainAlgo.getDenoisedConvergenceMap(workdir, filenames[i], outConvergenceMap);
             if ((outConvergenceMap.string()).empty() == false) {
               outfile << ",";
               outfile << outConvergenceMap.filename();
             }
             outConvergenceMap.clear();
           }
           if ((filenames[i].string()).find("NReSample") != std::string::npos) {
             NreSamples_Exists = true;
           }
         }// end of for loop
         if (NreSamples_Exists) {
          boost::filesystem::path NReConvMaps = CartesainAlgo.getResampledMaps (workdir, filenames);
          //outConvergenceMap.clear();
          CartesainAlgo.getSNRMap(workdir, NReConvMaps, outConvergenceMap);
          outfile << ",";
          outfile << outConvergenceMap.filename();
         }
         outfile << "]";
         outfile.close();
      } //end of else (input Shear map in json file)
    } //end of else
    logger.info("Done!");

////////////////////////////////////////////////////////////////////////////////////////////////////////
 //Write XML files
////////////////////////////////////////////////////////////////////////////////////////////////////////
  fs::path outputXML {workdir/out_xml_product_file_name};
  fs::path outputjson {workdir/outputConvMaps};
  CartesainAlgo.writeXMLfile (outputjson, outputXML);

////////////////////////////////////////////////////////////////////////////////////////////////////////
 //End of Mass Mapping
////////////////////////////////////////////////////////////////////////////////////////////////////////
   auto end_MappingTime = std::chrono::system_clock::now();
   std::chrono::duration<double> elapsed_seconds = end_MappingTime-start_MappingTime;
   std::time_t end_time = std::chrono::system_clock::to_time_t(end_MappingTime);
   logger.info() << "finished Program at " << std::ctime(&end_time);
   logger.info() << "elapsed time: " << elapsed_seconds.count() << "s";


    logger.info("#");
    logger.info("# Exiting mainMethod()");
    logger.info("#");

    return Elements::ExitCode::OK;
  }

};

MAIN_FOR(LE3_2D_MASS_WL_CartesianMassMapping)
